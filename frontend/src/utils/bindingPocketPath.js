/**
 * Smooth closed binding-pocket outline (Maestro LID–style) from residue positions.
 *
 * - `buildBindingPocketPath`: polar sort by angle → uniform radial expand → arc gaps → spline
 * - `buildWeightedPocketContourPath`: distance filter → polar → **variable radius**
 *   (`r + offset × interaction weight`) → angle sort → arc gaps → Catmull–Rom → Bézier
 * - `buildOrderedPocketContourPath`: **fixed** boundary order (no angle sort)
 */

const TAU = Math.PI * 2;

/**
 * Per-interaction multipliers for envelope radius:
 * `r_out = r_original + radialOffset * weight` (stronger contacts → slightly larger bump).
 */
export const DEFAULT_INTERACTION_WEIGHTS = {
  hbond: 1.25,
  salt_bridge: 1.3,
  pi_pi: 1.2,
  pi_cation: 1.2,
  hydrophobic: 1.05,
  default: 1.0,
};

function clamp(v, lo, hi) {
  return Math.max(lo, Math.min(hi, v));
}

function hypot2(ax, ay, bx, by) {
  const dx = bx - ax;
  const dy = by - ay;
  return Math.hypot(dx, dy);
}

/**
 * @param {number} x
 * @param {number} y
 * @param {number} cx
 * @param {number} cy
 * @returns {{ r: number, theta: number }} theta in (-π, π]
 */
export function toPolar(x, y, cx, cy) {
  const dx = x - cx;
  const dy = y - cy;
  return { r: Math.hypot(dx, dy), theta: Math.atan2(dy, dx) };
}

/**
 * @param {number} r
 * @param {number} theta
 * @param {number} cx
 * @param {number} cy
 */
export function toCartesian(r, theta, cx, cy) {
  return {
    x: cx + r * Math.cos(theta),
    y: cy + r * Math.sin(theta),
  };
}

/**
 * Sort residues by polar angle; stable for duplicate angles (by radius).
 * Theta normalized to [0, 2π).
 */
export function sortResiduesByAngle(residues, ligandCenter) {
  const cx = ligandCenter[0];
  const cy = ligandCenter[1];
  const tagged = residues.map((p, i) => {
    const { r, theta } = toPolar(p.x, p.y, cx, cy);
    const th = (theta + TAU) % TAU;
    return { x: p.x, y: p.y, r, theta: th, i, resid: p.resid, resname: p.resname };
  });
  tagged.sort((a, b) => {
    if (a.theta !== b.theta) return a.theta - b.theta;
    return a.r - b.r;
  });
  return tagged;
}

/** Match backend `interactions[].residue` e.g. "TYR156" */
export function matchResidueKey(residue) {
  const name = String(residue.resname ?? "").trim();
  const id = residue.resid != null ? String(residue.resid) : "";
  return `${name}${id}`;
}

/**
 * Max weight for a residue across matching interactions (strongest shapes envelope).
 * Unknown `type` uses `weights.default` or 1.
 */
export function maxInteractionWeightForResidue(residue, interactions, weights = DEFAULT_INTERACTION_WEIGHTS) {
  const wmap = { ...DEFAULT_INTERACTION_WEIGHTS, ...weights };
  const base = wmap.default ?? 1.0;
  if (!interactions?.length) return base;
  const key = matchResidueKey(residue);
  let maxW = base;
  for (const it of interactions) {
    if (String(it.residue) !== key) continue;
    const t = String(it.type ?? "").toLowerCase();
    const w = wmap[t] ?? base;
    if (w > maxW) maxW = w;
  }
  return maxW;
}

/** Keep residues whose 2D distance to ligand center is ≤ maxDistancePx (screen space). */
export function filterResiduesNearLigand(residues, ligandCenter, maxDistancePx) {
  const cx = Number(ligandCenter[0]);
  const cy = Number(ligandCenter[1]);
  const maxD =
    maxDistancePx == null || !Number.isFinite(maxDistancePx)
      ? Infinity
      : Math.max(0, Number(maxDistancePx));
  return residues.filter((r) => Math.hypot(Number(r.x) - cx, Number(r.y) - cy) <= maxD);
}

/**
 * Insert points along circular arcs when angular gaps are large (non-uniform distribution).
 *
 * @param {Array<{x:number,y:number}>} ordered - CCW ordered ring
 * @param {[number, number]} ligandCenter
 * @param {number} maxAngleStep - radians between inserted samples (e.g. π/16)
 */
export function interpolateArcGaps(ordered, ligandCenter, maxAngleStep = Math.PI / 16) {
  if (ordered.length < 2) return ordered.slice();
  const cx = ligandCenter[0];
  const cy = ligandCenter[1];
  const out = [];

  for (let i = 0; i < ordered.length; i++) {
    const a = ordered[i];
    const b = ordered[(i + 1) % ordered.length];
    out.push({ x: a.x, y: a.y });

    const ta = Math.atan2(a.y - cy, a.x - cx);
    const tb = Math.atan2(b.y - cy, b.x - cx);
    let d = tb - ta;
    while (d <= -Math.PI) d += TAU;
    while (d > Math.PI) d -= TAU;

    const steps = Math.ceil(Math.abs(d) / maxAngleStep) - 1;
    if (steps > 0) {
      const ra = Math.hypot(a.x - cx, a.y - cy);
      const rb = Math.hypot(b.x - cx, b.y - cy);
      for (let s = 1; s <= steps; s++) {
        const t = s / (steps + 1);
        const th = ta + d * t;
        const r = ra + (rb - ra) * t;
        out.push(toCartesian(r, th, cx, cy));
      }
    }
  }

  return out;
}

/**
 * Chord length raised to α (centripetal uses α = 0.5).
 */
function chordAlpha(p0, p1, alpha) {
  const d = hypot2(p0.x, p0.y, p1.x, p1.y);
  if (d < 1e-12) return 0;
  return Math.pow(d, alpha);
}

/**
 * Uniform Catmull–Rom segment (p1→p2) as cubic Bézier, with chord-based tension
 * scaling (centripetal-style: α = 0.5 reduces cusps vs uniform α = 1).
 * Handles are clamped to limit curvature (helps avoid self-intersection on dense data).
 */
function catmullRomToBezier(p0, p1, p2, p3, alpha, tension, curvatureClamp) {
  const t01 = chordAlpha(p0, p1, alpha);
  const t12 = chordAlpha(p1, p2, alpha);
  const t23 = chordAlpha(p2, p3, alpha);
  const denom = t01 + t23 + 1e-9;
  const scale = 1 + (t12 / denom) * 0.35;
  const s = 6 * tension * scale;

  let cp1x = p1.x + (p2.x - p0.x) / s;
  let cp1y = p1.y + (p2.y - p0.y) / s;
  let cp2x = p2.x - (p3.x - p1.x) / s;
  let cp2y = p2.y - (p3.y - p1.y) / s;

  const chord = hypot2(p1.x, p1.y, p2.x, p2.y);
  const maxH = chord * curvatureClamp;
  const clampHandle = (px, py, ax, ay) => {
    const vx = px - ax;
    const vy = py - ay;
    const L = Math.hypot(vx, vy);
    if (L <= maxH || L < 1e-12) return { x: px, y: py };
    const sc = maxH / L;
    return { x: ax + vx * sc, y: ay + vy * sc };
  };

  const c1 = clampHandle(cp1x, cp1y, p1.x, p1.y);
  const c2 = clampHandle(cp2x, cp2y, p2.x, p2.y);
  return { cp1: c1, cp2: c2, end: p2 };
}

function fmt(n) {
  return Number.isFinite(n) ? (Math.round(n * 1000) / 1000).toString() : "0";
}

/**
 * Shift each point radially from the ligand center (keeps boundary order).
 */
export function expandAnchorsOutward(points, ligandCenter, outwardPx) {
  const cx = ligandCenter[0];
  const cy = ligandCenter[1];
  const d = Number(outwardPx) || 0;
  if (d === 0) return points.map((p) => ({ x: p.x, y: p.y }));
  return points.map((p) => {
    const dx = p.x - cx;
    const dy = p.y - cy;
    const L = Math.hypot(dx, dy) || 1;
    return { x: p.x + (dx / L) * d, y: p.y + (dy / L) * d };
  });
}

/**
 * Remove consecutive duplicates (stabilizes splines when two anchors coincide).
 */
function dedupeConsecutivePoints(points, eps = 1e-3) {
  const out = [];
  for (const p of points) {
    const last = out[out.length - 1];
    if (last && hypot2(last.x, last.y, p.x, p.y) < eps) continue;
    out.push({ x: p.x, y: p.y });
  }
  return out;
}

/**
 * Insert points along each closed-ring edge so long chords are subdivided (smoother CR spline).
 * `maxSegmentLength` = target maximum Euclidean length between consecutive samples (px).
 */
export function densifyOrderedClosedRing(points, maxSegmentLength) {
  const n = points.length;
  if (n < 2 || !maxSegmentLength || maxSegmentLength <= 0) {
    return points.map((p) => ({ x: p.x, y: p.y }));
  }
  const out = [];
  for (let i = 0; i < n; i++) {
    const a = points[i];
    const b = points[(i + 1) % n];
    out.push({ x: a.x, y: a.y });
    const d = hypot2(a.x, a.y, b.x, b.y);
    const extra = Math.max(0, Math.ceil(d / maxSegmentLength) - 1);
    for (let k = 1; k <= extra; k++) {
      const t = k / (extra + 1);
      out.push({
        x: a.x + t * (b.x - a.x),
        y: a.y + t * (b.y - a.y),
      });
    }
  }
  return out;
}

/**
 * Closed Catmull–Rom chain → cubic Bézier SVG path (no implicit reordering of points).
 */
function svgPathClosedCatmullRom(pts, centripetalAlpha, tension, curvatureClamp) {
  const n = pts.length;
  if (n < 3) return null;
  const parts = [`M ${fmt(pts[0].x)} ${fmt(pts[0].y)}`];
  for (let i = 0; i < n; i++) {
    const p0 = pts[(i - 1 + n) % n];
    const p1 = pts[i];
    const p2 = pts[(i + 1) % n];
    const p3 = pts[(i + 2) % n];
    const { cp1, cp2, end } = catmullRomToBezier(
      p0,
      p1,
      p2,
      p3,
      centripetalAlpha,
      tension,
      curvatureClamp
    );
    parts.push(
      `C ${fmt(cp1.x)} ${fmt(cp1.y)} ${fmt(cp2.x)} ${fmt(cp2.y)} ${fmt(end.x)} ${fmt(end.y)}`
    );
  }
  parts.push("Z");
  return parts.join(" ");
}

function pathForTwoPoints(p0, p1, cx, cy, pad) {
  const mx = (p0.x + p1.x) / 2;
  const my = (p0.y + p1.y) / 2;
  const r = Math.max(hypot2(mx, my, cx, cy) * 0.08, pad);
  const d0 = Math.hypot(p0.x - cx, p0.y - cy);
  const d1 = Math.hypot(p1.x - cx, p1.y - cy);
  const rr = Math.max(d0, d1, r);
  return `M ${fmt(p0.x)} ${fmt(p0.y)} A ${fmt(rr)} ${fmt(rr)} 0 1 1 ${fmt(p1.x)} ${fmt(p1.y)} A ${fmt(rr)} ${fmt(rr)} 0 1 1 ${fmt(p0.x)} ${fmt(p0.y)}`;
}

/**
 * Dynamic pocket envelope: distance filter → polar → variable radius
 * `r_out = r + radialOffset * interactionWeight` → sort by angle → arc interpolation → spline.
 *
 * @param {Array<{x:number,y:number,resname?:string,resid?:number|string}>} residues
 * @param {[number, number]} ligandCenter
 * @param {Array<{ residue: string, type: string }>} [interactions]
 * @param {object} [options]
 * @param {number} [options.distanceThreshold] - max px from ligand center; default = no filter
 * @param {number} [options.radialOffset=14] - scales with interaction weight (loose hug)
 * @param {Record<string, number>} [options.interactionWeights] - overrides {@link DEFAULT_INTERACTION_WEIGHTS}
 * @param {number} [options.maxAngleStep] - arc gap fill (radians)
 * @param {number} [options.maxSegmentLength=0] - optional chord densification after arc step
 * @param {number} [options.minRadius=0] - floor on r_out (avoid degenerate center hits)
 */
export function buildWeightedPocketContourPath(
  residues,
  ligandCenter,
  interactions,
  options = {}
) {
  let inter = interactions;
  let opts = options;
  if (interactions != null && typeof interactions === "object" && !Array.isArray(interactions)) {
    opts = interactions;
    inter = [];
  }

  const cx = Number(ligandCenter[0]);
  const cy = Number(ligandCenter[1]);

  const radialOffset = Math.max(0, Number(opts.radialOffset ?? opts.offset ?? 14));
  const distanceThreshold = opts.distanceThreshold;
  const smoothing = clamp(opts.smoothing ?? 1, 0.35, 2.5);
  const centripetalAlpha = clamp(opts.centripetalAlpha ?? 0.5, 0.25, 1);
  const curvatureClamp = clamp(opts.curvatureClamp ?? 0.42, 0.12, 0.48);
  const maxAngleStep = opts.maxAngleStep ?? Math.PI / 14;
  const minRadius = Math.max(0, Number(opts.minRadius ?? 0));
  const maxSegmentLength = Math.max(
    0,
    Number(opts.maxSegmentLength ?? opts.interpolationDensity ?? 0)
  );
  const weights = { ...DEFAULT_INTERACTION_WEIGHTS, ...(opts.interactionWeights || {}) };

  let filtered = residues;
  if (distanceThreshold != null && Number.isFinite(Number(distanceThreshold))) {
    filtered = filterResiduesNearLigand(residues, [cx, cy], Number(distanceThreshold));
  }
  if (!filtered?.length) return "";

  const tagged = sortResiduesByAngle(filtered, [cx, cy]);
  const ring = tagged.map((p) => {
    const res = {
      resname: p.resname,
      resid: p.resid,
    };
    const w = maxInteractionWeightForResidue(res, inter || [], weights);
    const r0 = Math.max(p.r, 1e-6);
    const rOut = Math.max(minRadius, r0 + radialOffset * w);
    return toCartesian(rOut, p.theta, cx, cy);
  });

  if (ring.length === 1) {
    const r = Math.max(
      tagged[0].r +
        radialOffset *
          maxInteractionWeightForResidue(
            { resname: tagged[0].resname, resid: tagged[0].resid },
            inter || [],
            weights
          ),
      radialOffset * 2 || 12
    );
    const x = cx + r;
    return `M ${fmt(x)} ${fmt(cy)} A ${fmt(r)} ${fmt(r)} 0 1 1 ${fmt(cx - r)} ${fmt(cy)} A ${fmt(r)} ${fmt(r)} 0 1 1 ${fmt(x)} ${fmt(cy)}`;
  }

  if (ring.length === 2) {
    return pathForTwoPoints(ring[0], ring[1], cx, cy, radialOffset || 8);
  }

  let pts = interpolateArcGaps(ring, [cx, cy], maxAngleStep);
  if (pts.length < 3) pts = ring;

  if (maxSegmentLength > 0) {
    pts = densifyOrderedClosedRing(pts, maxSegmentLength);
    pts = dedupeConsecutivePoints(pts);
    if (pts.length < 3) {
      pts = interpolateArcGaps(ring, [cx, cy], maxAngleStep * 0.75);
    }
  }

  return (
    svgPathClosedCatmullRom(pts, centripetalAlpha, smoothing, curvatureClamp) || ""
  );
}

/**
 * Build a smooth closed SVG path around the ligand from residue screen positions.
 *
 * @param {Array<{x:number,y:number,resid?:number,resname?:string}>} residues
 * @param {[number, number]} ligandCenter - [cx, cy]
 * @param {object} [options]
 * @param {number} [options.radiusScale=1] - multiply polar radius after sorting
 * @param {number} [options.radialPadding=8] - extra outward shift (px) after scale
 * @param {number} [options.smoothing=1] - higher → gentler bends (less tight); typical 0.7–1.4
 * @param {number} [options.centripetalAlpha=0.5] - 0.5 = standard centripetal
 * @param {number} [options.maxAngleStep] - radians; insert arc samples if gap larger (default π/14)
 * @param {number} [options.curvatureClamp=0.42] - max control arm length as fraction of chord
 * @returns {string} SVG `d` attribute, or "" if not enough points
 */
export function buildBindingPocketPath(residues, ligandCenter, options = {}) {
  const cx = Number(ligandCenter[0]);
  const cy = Number(ligandCenter[1]);

  const radiusScale = options.radiusScale ?? 1;
  const radialPadding = options.radialPadding ?? 8;
  const smoothing = clamp(options.smoothing ?? 1, 0.35, 2.5);
  const centripetalAlpha = clamp(options.centripetalAlpha ?? 0.5, 0.25, 1);
  const maxAngleStep = options.maxAngleStep ?? Math.PI / 14;
  const tension = smoothing;
  const curvatureClamp = clamp(options.curvatureClamp ?? 0.42, 0.12, 0.48);

  if (!residues?.length) return "";

  const sorted = sortResiduesByAngle(residues, [cx, cy]).map((p) => {
    const rNew = p.r * radiusScale + radialPadding;
    const c = toCartesian(rNew, p.theta, cx, cy);
    return { x: c.x, y: c.y };
  });

  if (sorted.length === 1) {
    const r0 = Math.hypot(residues[0].x - cx, residues[0].y - cy);
    const r = Math.max(r0 * radiusScale + radialPadding, radialPadding * 2);
    const x = cx + r;
    return `M ${fmt(x)} ${fmt(cy)} A ${fmt(r)} ${fmt(r)} 0 1 1 ${fmt(cx - r)} ${fmt(cy)} A ${fmt(r)} ${fmt(r)} 0 1 1 ${fmt(x)} ${fmt(cy)}`;
  }

  if (sorted.length === 2) {
    return pathForTwoPoints(sorted[0], sorted[1], cx, cy, radialPadding);
  }

  let pts = interpolateArcGaps(sorted, [cx, cy], maxAngleStep);
  if (pts.length < 3) pts = sorted;

  return (
    svgPathClosedCatmullRom(pts, centripetalAlpha, tension, curvatureClamp) || ""
  );
}

/**
 * Pocket contour from an **explicit** ordered boundary list (e.g. Maestro-style walk around the site).
 * Does **not** re-sort by angle; order is preserved exactly.
 *
 * @example
 * const boundary = [
 *   { resname: "TYR", resid: 156, x: 120, y: 400, score: "+0.000" },
 *   { resname: "LEU", resid: 193, x: 140, y: 380, score: "+0.000" },
 * ];
 * const d = buildOrderedPocketContourPath(boundary, ligand_center, {
 *   outwardOffset: 12,
 *   smoothing: 1.05,
 *   maxSegmentLength: 28,
 * });
 * // <path d={d} fill="none" stroke="url(#pocketBandGrad)" strokeWidth={6} />
 *
 * @param {Array<{x:number,y:number,resname?:string,resid?:number|string,score?:string}>} orderedResidues
 * @param {[number, number]} ligandCenter - [cx, cy]
 * @param {object} [options]
 * @param {number} [options.outwardOffset=10] - px to push the contour outward from ligand center (per anchor)
 * @param {number} [options.smoothing=1] - spline tension (higher → gentler turns); alias: smoothingStrength
 * @param {number} [options.centripetalAlpha=0.5] - chord weighting (centripetal character)
 * @param {number} [options.curvatureClamp=0.42] - limit Bézier handle length vs chord (reduces loops / cusps)
 * @param {number} [options.maxSegmentLength=0] - if > 0, subdivide each boundary edge so chord length ≤ this (px); alias: interpolationDensity, maxChordLength
 * @returns {string} SVG `d` attribute
 */
export function buildOrderedPocketContourPath(orderedResidues, ligandCenter, options = {}) {
  const cx = Number(ligandCenter[0]);
  const cy = Number(ligandCenter[1]);

  const outwardOffset =
    options.outwardOffset ??
    options.outwardPx ??
    options.radialPadding ??
    10;
  const smoothing = clamp(
    options.smoothing ?? options.smoothingStrength ?? 1,
    0.35,
    2.5
  );
  const centripetalAlpha = clamp(options.centripetalAlpha ?? 0.5, 0.25, 1);
  const curvatureClamp = clamp(options.curvatureClamp ?? 0.42, 0.12, 0.48);
  const maxSegmentLength = Math.max(
    0,
    Number(
      options.maxSegmentLength ??
        options.interpolationDensity ??
        options.maxChordLength ??
        0
    )
  );

  if (!orderedResidues?.length) return "";

  let anchors = orderedResidues.map((r) => ({
    x: Number(r.x),
    y: Number(r.y),
  }));
  anchors = dedupeConsecutivePoints(anchors);
  if (anchors.length === 0) return "";

  let pts = expandAnchorsOutward(anchors, [cx, cy], outwardOffset);
  pts = dedupeConsecutivePoints(pts);

  if (maxSegmentLength > 0) {
    pts = densifyOrderedClosedRing(pts, maxSegmentLength);
    pts = dedupeConsecutivePoints(pts);
  }

  if (pts.length === 1) {
    const r = Math.max(
      Math.hypot(anchors[0].x - cx, anchors[0].y - cy) + outwardOffset,
      Math.abs(outwardOffset) * 2 || 12
    );
    const x = cx + r;
    return `M ${fmt(x)} ${fmt(cy)} A ${fmt(r)} ${fmt(r)} 0 1 1 ${fmt(cx - r)} ${fmt(cy)} A ${fmt(r)} ${fmt(r)} 0 1 1 ${fmt(x)} ${fmt(cy)}`;
  }

  if (pts.length === 2) {
    return pathForTwoPoints(pts[0], pts[1], cx, cy, Math.abs(outwardOffset) || 8);
  }

  return svgPathClosedCatmullRom(pts, centripetalAlpha, smoothing, curvatureClamp) || "";
}

export default buildBindingPocketPath;
