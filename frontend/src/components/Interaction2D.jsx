import React from "react";

/** Unit teardrop: tip points +local Y (down); rotate so tip aims at ligand center. */
const TEARDROP_D =
  "M 0,-11 C 7.2,-11 11.2,-6 11.2,-0.4 C 11.2,6.6 0,14.2 0,14.2 C 0,14.2 -11.2,6.6 -11.2,-0.4 C -11.2,-6 -7.2,-11 0,-11 Z";

const LegendItem = ({ color, label, line, dashed, arrow, borderColor, gradientColors }) => (
  <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
    {line ? (
      <div style={{ position: "relative", width: 30, height: 20, display: "flex", alignItems: "center" }}>
        {gradientColors && gradientColors.length >= 2 ? (
          <div
            style={{
              width: arrow ? 20 : 26,
              height: 3,
              background: `linear-gradient(to right, ${gradientColors[0]}, ${gradientColors[1]})`,
              borderRadius: 1,
            }}
          />
        ) : (
          <>
            <div
              style={{
                width: arrow ? 20 : 26,
                height: 0,
                borderTop: dashed ? `3px dashed ${color}` : `3px solid ${color}`,
              }}
            />
            {arrow && (
              <div
                style={{
                  width: 0,
                  height: 0,
                  borderLeft: `6px solid ${color}`,
                  borderTop: "4px solid transparent",
                  borderBottom: "4px solid transparent",
                  marginLeft: "-2px",
                }}
              />
            )}
          </>
        )}
      </div>
    ) : (
      <div 
        style={{ 
          width: 14, 
          height: 14, 
          background: color, 
          borderRadius: "50%", 
          border: `1px solid ${borderColor || "#cbd5e1"}` 
        }} 
      />
    )}
    <span style={{ color: "#0f172a" }}>{label}</span>
  </div>
);

function Interaction2D({ data }) {
  const safeData = data || {};
  const {
    svg,
    residues = [],
    pocket_outline_path,
    pocket_outline_paths,
    backbone_path,
    backbone_gap_dots = [],
    meta = {},
    interactions = [],
    ligand_atom_xy,
    ligand_atom_solvent_exposed,
    ligand_center,
    conformer_info,
  } = safeData;

  if (!data) return null;

  // Server layout + pocket share the same geometry; do not force a client-side circle.
  const positionedResidues = residues;

  const interactionResidueKeys = new Set(
    (interactions || []).map((it) => it.residue || `${it.resname}${it.resid}`)
  );
  const interactionTotal =
    meta.interaction_total != null ? Number(meta.interaction_total) : null;
  const interactionFiltered =
    meta.interaction_filtered != null
      ? Number(meta.interaction_filtered)
      : interactions.length;
  const stats = {
    residueNodes: positionedResidues.length,
    interactionTotal,
    interactionFiltered,
    uniqueInteractionResidues: interactionResidueKeys.size,
  };

  const residueIndex = (() => {
    const map = new Map();
    positionedResidues.forEach((r) => map.set(`${r.resname}${r.resid}`, r));
    return map;
  })();

  const residueGradientId = (cls) => {
    if (cls === "negative") return "lidGradNegative";
    if (cls === "positive") return "lidGradPositive";
    if (cls === "hydrophobic") return "lidGradHydrophobic";
    if (cls === "polar") return "lidGradPolar";
    if (cls === "cysteine") return "lidGradCys";
    if (cls === "glycine") return "lidGradGly";
    if (cls === "metal") return "lidGradMetal";
    if (cls === "nucleic") return "lidGradNucleic";
    return "lidGradOther";
  };

  const pocketPaths =
    Array.isArray(pocket_outline_paths) && pocket_outline_paths.length > 0
      ? pocket_outline_paths
      : pocket_outline_path
        ? [pocket_outline_path]
        : [];

  const styleFor = (it) => {
    switch (it.type) {
      case "hbond":
      {
        // Magenta arrow for H-bond (Maestro LID reference)
        const isBackbone = it.backbone;
        return {
          stroke: "#c026d3",
          strokeWidth: 2,
          dash: isBackbone ? "" : "4 2",
          markerEnd: "url(#arrowMagenta)",
          opacity: 0.95,
        };
      }
      case "pi_pi":
        // Green line for π–π stacking (Maestro LID)
        return {
          stroke: "#22c55e",
          strokeWidth: 2.5,
          dash: "",
          markerEnd: "",
          opacity: 0.9,
          showCircles: true,
          circleCount: 2
        };
      case "pi_cation":
        // Red line for π–cation (Maestro LID)
        return {
          stroke: "#ef4444",
          strokeWidth: 2.5,
          dash: "",
          markerEnd: "",
          opacity: 0.95,
          showCircles: true,
          circleCount: 1
        };
      case "halogen_bond":
        // Yellow line/arrow for halogen bond (Maestro LID)
        return { stroke: "#eab308", strokeWidth: 2, dash: "", markerEnd: "url(#arrowYellow)", opacity: 0.95 };
      case "metal_coordination":
        // Purple line for metal coordination (Maestro LID)
        return { stroke: "#7c3aed", strokeWidth: 2.5, dash: "", markerEnd: "", opacity: 0.95 };
      case "salt_bridge":
        // Red–blue gradient for salt bridge (Maestro LID); stroke uses gradient id
        return { stroke: "url(#saltBridgeGradient)", strokeWidth: 2.5, dash: "", markerEnd: "", opacity: 0.95 };
      case "hydrophobic":
        return { stroke: "#000000", strokeWidth: 1, dash: "", markerEnd: "", opacity: 0.5 };
      case "contact":
        return { stroke: "#64748b", strokeWidth: 1, dash: "3 2", markerEnd: "", opacity: 0.55 };
      case "distance":
        return { stroke: "#22c55e", strokeWidth: 1, dash: "2 2", markerEnd: "", opacity: 0.4 };
      default:
        return { stroke: "#9ca3af", strokeWidth: 1, dash: "2 2", markerEnd: "", opacity: 0.3 };
    }
  };

  const ligandAnchor = (idx) => {
    if (idx == null || idx < 0 || !ligand_atom_xy || !ligand_atom_xy[idx]) return null;
    return ligand_atom_xy[idx];
  };

  return (
    <div
      className="w-[95vw] mx-auto     p-5  "
      style={{ maxWidth: meta.svg_w + 100 }}
    >
      <div className="mb-4 flex flex-wrap items-center justify-between gap-3  pb-4">
        <div>
          <h2 className="text-base font-semibold text-slate-900">Ligand interaction diagram</h2>
          <p className="text-xs text-slate-500">
            Maestro-style LID: 2D ligand, teardrop residues toward the ligand, gradient pocket band with
            entrance gap, dotted backbone, magenta H-bond arrows, solvent-exposure markers.
          </p>
        </div>
        <dl className="flex flex-wrap gap-x-6 gap-y-2 text-sm">
          <div className="flex items-baseline gap-2">
            <dt className="text-slate-500">Residues shown</dt>
            <dd className="tabular-nums font-semibold text-slate-900">{stats.residueNodes}</dd>
          </div>
          <div className="flex items-baseline gap-2">
            <dt className="text-slate-500">Interactions (detected)</dt>
            <dd className="tabular-nums font-semibold text-slate-900">
              {stats.interactionTotal != null && !Number.isNaN(stats.interactionTotal)
                ? stats.interactionTotal
                : "—"}
            </dd>
          </div>
          <div className="flex items-baseline gap-2">
            <dt className="text-slate-500">Interactions (filtered)</dt>
            <dd className="tabular-nums font-semibold text-slate-900">{stats.interactionFiltered}</dd>
          </div>
          <div className="flex items-baseline gap-2">
            <dt className="text-slate-500">Unique contact residues</dt>
            <dd className="tabular-nums font-semibold text-slate-900">{stats.uniqueInteractionResidues}</dd>
          </div>
        </dl>
      </div>

      <div style={{ position: "relative", width: meta.svg_w, height: meta.svg_h, overflow: "visible" }}>
        {/* RDKit SVG - positioned at ligand center, BEHIND protein structures (zIndex: 0) */}
        {ligand_center && ligand_center.length === 2 ? (
          <div
            style={{
              position: "absolute",
              left: `${ligand_center[0]}px`,
              top: `${ligand_center[1]}px`,
              transform: "translate(-50%, -50%) scale(0.74)",
              pointerEvents: "none",
              zIndex: 0,
              opacity: 1,
            }}
            dangerouslySetInnerHTML={{ __html: svg }}
          />
        ) : null}

        {/* SVG with protein structures, interactions, and residues (zIndex: 1) */}
        <svg width={meta.svg_w} height={meta.svg_h} style={{ position: "absolute", left: 0, top: 0, overflow: "visible", zIndex: 1 }}>
          <defs>
            <marker id="arrowMagenta" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
              <path d="M 0 0 L 10 5 L 0 10 z" fill="#c026d3" />
            </marker>
            <marker id="arrowYellow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
              <path d="M 0 0 L 10 5 L 0 10 z" fill="#eab308" />
            </marker>
            <linearGradient id="saltBridgeGradient" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stopColor="#dc2626" />
              <stop offset="100%" stopColor="#2563eb" />
            </linearGradient>
            {/* Pocket band (Maestro: green → yellow → cyan glow) */}
            <linearGradient id="pocketBandGrad" x1="0%" y1="0%" x2="100%" y2="100%">
              <stop offset="0%" stopColor="#4ade80" />
              <stop offset="42%" stopColor="#facc15" />
              <stop offset="100%" stopColor="#38bdf8" />
            </linearGradient>
            {/* Residue “bubble” lighting (top highlight → deeper base) */}
            <linearGradient id="lidGradHydrophobic" x1="0" y1="0" x2="0" y2="1" gradientUnits="objectBoundingBox">
              <stop offset="0%" stopColor="#86efac" />
              <stop offset="55%" stopColor="#22c55e" />
              <stop offset="100%" stopColor="#15803d" />
            </linearGradient>
            <linearGradient id="lidGradPolar" x1="0" y1="0" x2="0" y2="1" gradientUnits="objectBoundingBox">
              <stop offset="0%" stopColor="#bae6fd" />
              <stop offset="55%" stopColor="#38bdf8" />
              <stop offset="100%" stopColor="#0284c7" />
            </linearGradient>
            <linearGradient id="lidGradPositive" x1="0" y1="0" x2="0" y2="1" gradientUnits="objectBoundingBox">
              <stop offset="0%" stopColor="#c4b5fd" />
              <stop offset="55%" stopColor="#7c3aed" />
              <stop offset="100%" stopColor="#5b21b6" />
            </linearGradient>
            <linearGradient id="lidGradNegative" x1="0" y1="0" x2="0" y2="1" gradientUnits="objectBoundingBox">
              <stop offset="0%" stopColor="#fda4af" />
              <stop offset="55%" stopColor="#f43f5e" />
              <stop offset="100%" stopColor="#be123c" />
            </linearGradient>
            <linearGradient id="lidGradCys" x1="0" y1="0" x2="0" y2="1" gradientUnits="objectBoundingBox">
              <stop offset="0%" stopColor="#d9f99d" />
              <stop offset="100%" stopColor="#65a30d" />
            </linearGradient>
            <linearGradient id="lidGradGly" x1="0" y1="0" x2="0" y2="1" gradientUnits="objectBoundingBox">
              <stop offset="0%" stopColor="#f9fafb" />
              <stop offset="100%" stopColor="#d1d5db" />
            </linearGradient>
            <linearGradient id="lidGradMetal" x1="0" y1="0" x2="0" y2="1" gradientUnits="objectBoundingBox">
              <stop offset="0%" stopColor="#e5e7eb" />
              <stop offset="100%" stopColor="#9ca3af" />
            </linearGradient>
            <linearGradient id="lidGradNucleic" x1="0" y1="0" x2="0" y2="1" gradientUnits="objectBoundingBox">
              <stop offset="0%" stopColor="#e5e7eb" />
              <stop offset="100%" stopColor="#9ca3af" />
            </linearGradient>
            <linearGradient id="lidGradOther" x1="0" y1="0" x2="0" y2="1" gradientUnits="objectBoundingBox">
              <stop offset="0%" stopColor="#f1f5f9" />
              <stop offset="100%" stopColor="#94a3b8" />
            </linearGradient>
          </defs>

          {/* Pocket boundary: soft glow + dashed gradient band (Maestro LID) */}
          {pocketPaths.length > 0 ? (
            <g style={{ filter: "drop-shadow(0 0 6px rgba(56, 189, 248, 0.45))" }}>
              {pocketPaths.map((d, pi) => (
                <g key={`pocket-${pi}`} opacity={0.95}>
                  <path
                    d={d}
                    fill="none"
                    stroke="url(#pocketBandGrad)"
                    strokeWidth={9}
                    strokeOpacity={0.22}
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                  <path
                    d={d}
                    fill="none"
                    stroke="#334155"
                    strokeWidth={2.4}
                    strokeOpacity={0.35}
                    strokeDasharray="4 9"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                  <path
                    d={d}
                    fill="none"
                    stroke="url(#pocketBandGrad)"
                    strokeWidth={4.2}
                    strokeDasharray="4 9"
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                </g>
              ))}
            </g>
          ) : null}

          {/* Smooth protein backbone + bead dashes (Maestro chain) */}
          {backbone_path && String(backbone_path).trim().length > 0 ? (
            <g>
              <path
                d={backbone_path}
                fill="none"
                stroke="#0f172a"
                strokeWidth={3.2}
                opacity={0.92}
                strokeLinecap="round"
                strokeLinejoin="round"
              />
              <path
                d={backbone_path}
                fill="none"
                stroke="#0f172a"
                strokeWidth={3.2}
                opacity={0.55}
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeDasharray="1 13"
              />
              {Array.isArray(backbone_gap_dots) &&
                backbone_gap_dots.map((p, gi) => (
                  <circle
                    key={`bb-gap-${gi}`}
                    cx={p.x}
                    cy={p.y}
                    r={2.8}
                    fill="#0f172a"
                    stroke="#0f172a"
                    strokeWidth={0.4}
                    opacity={0.92}
                  />
                ))}
            </g>
          ) : null}

          {/* Interaction lines: curved from ligand atom to residue, with height spacing */}
          {(() => {
            if (!ligand_center || ligand_center.length !== 2) return null;
            const lcx = ligand_center[0];
            const lcy = ligand_center[1];
            // Angle from ligand center to each residue node (for grouping nearby directions)
            const withAngle = interactions.map((it, idx) => {
              const node = residueIndex.get(it.residue);
              const anchor = ligandAnchor(it.ligand_atom_index);
              if (!node || !anchor) return null;
              const angle = Math.atan2(node.y - lcy, node.x - lcx);
              return { it, node, anchor, angle, idx };
            }).filter(Boolean);

            // Assign "lane" (height) index: interactions in similar direction get different lanes
            const angleBin = (Math.PI * 2) / 2; // ~15° bins
            const laneCount = {};
            const laneIndex = {};
            withAngle.forEach(({ angle, idx }) => {
              const key = Math.floor(angle / angleBin) * angleBin;
              const lane = laneCount[key] ?? 0;
              laneIndex[idx] = lane;
              laneCount[key] = lane + 1;
            });

            const heightSpacing = 100; // ~1 "cm" visual space between curved lanes
            const baseCurveOffset = 100;   // base bulge (curved on surface)

            return withAngle.map(({ it, node, anchor, idx }) => {
              const st = styleFor(it);
              const lane = laneIndex[idx] || 0;
              const offset = baseCurveOffset + lane * heightSpacing;

              const mx = (anchor.x + node.x) / 2;
              const my = (anchor.y + node.y) / 2;
              const dx = node.x - anchor.x;
              const dy = node.y - anchor.y;
              const len = Math.sqrt(dx * dx + dy * dy) || 1;
              const perpX = -dy / len;
              const perpY = dx / len;
              const fromLigandX = mx - lcx;
              const fromLigandY = my - lcy;
              const outward = perpX * fromLigandX + perpY * fromLigandY >= 0 ? 1 : -1;
              const cx = mx + perpX * offset * outward;
              const cy = my + perpY * offset * outward;

              const pathD = `M ${anchor.x} ${anchor.y} Q ${cx} ${cy} ${node.x} ${node.y}`;
              const midX = (anchor.x + cx + node.x) / 3;
              const midY = (anchor.y + cy + node.y) / 3;

              return (
                <g key={idx}>
                  <path
                    d={pathD}
                    fill="none"
                    stroke={st.stroke}
                    strokeWidth={st.strokeWidth}
                    strokeDasharray={st.dash}
                    opacity={st.opacity}
                    markerEnd={st.markerEnd}
                    style={{ filter: "drop-shadow(0 1px 1px rgba(0,0,0,0.1))" }}
                  />
                  {st.showCircles && st.circleCount && (() => {
                    const circles = [];
                    for (let i = 1; i <= st.circleCount; i++) {
                      const t = i / (st.circleCount + 1);
                      const t1 = 1 - t;
                      const x = t1 * t1 * anchor.x + 2 * t1 * t * cx + t * t * node.x;
                      const y = t1 * t1 * anchor.y + 2 * t1 * t * cy + t * t * node.y;
                      circles.push(
                        <circle key={i} cx={x} cy={y} r="3" fill={st.stroke} opacity={st.opacity} />
                      );
                    }
                    return circles;
                  })()}
                  {it.type === "distance" && (
                    <text
                      x={midX}
                      y={midY - 8}
                      textAnchor="middle"
                      fontSize="7"
                      fill="#22c55e"
                      fontWeight="500"
                      opacity={0.7}
                    >
                      {it.distance?.toFixed(1)}Å
                    </text>
                  )}
                </g>
              );
            });
          })()}

          {/* Fallback: per-segment curves only if API did not send backbone_path */}
          {(() => {
            if (backbone_path && String(backbone_path).trim().length > 0) return null;
            if (!ligand_center || ligand_center.length !== 2) return null;

            const connections = [];
            const lcx = ligand_center[0];
            const lcy = ligand_center[1];

            const residuesByChain = {};
            positionedResidues.forEach((r, idx) => {
              const chain = r.chain || "A";
              if (!residuesByChain[chain]) residuesByChain[chain] = [];
              residuesByChain[chain].push({ ...r, originalIndex: idx });
            });

            Object.keys(residuesByChain).forEach((chain) => {
              const chainResidues = residuesByChain[chain].sort(
                (a, b) => (a.resid || 0) - (b.resid || 0)
              );

              for (let i = 0; i < chainResidues.length - 1; i++) {
                const curr = chainResidues[i];
                const next = chainResidues[i + 1];
                const currResid = curr.resid || 0;
                const nextResid = next.resid || 0;
                const gap = nextResid - currResid;
                if (gap < 1) continue;
                // Only connect consecutive numbering, or a single skip (N → N+2) with one dot.
                if (gap !== 1 && gap !== 2) continue;

                const mx = (curr.x + next.x) / 2;
                const my = (curr.y + next.y) / 2;
                const dx = next.x - curr.x;
                const dy = next.y - curr.y;
                const len = Math.sqrt(dx * dx + dy * dy) || 1;
                const perpX = -dy / len;
                const perpY = dx / len;
                const fromLigandX = mx - lcx;
                const fromLigandY = my - lcy;
                const outward = perpX * fromLigandX + perpY * fromLigandY >= 0 ? 1 : -1;

                const curveOffset = 14;
                const cxCtrl = mx + perpX * curveOffset * outward;
                const cyCtrl = my + perpY * curveOffset * outward;

                const pathD = `M ${curr.x} ${curr.y} Q ${cxCtrl} ${cyCtrl} ${next.x} ${next.y}`;

                connections.push(
                  <path
                    key={`seq-${chain}-${currResid}-${nextResid}`}
                    d={pathD}
                    fill="none"
                    stroke="#0f172a"
                    strokeWidth="2"
                    opacity={0.88}
                    strokeLinecap="round"
                    strokeLinejoin="round"
                  />
                );

                if (gap === 2) {
                  const t = 0.5;
                  const omt = 1 - t;
                  const gx = omt * omt * curr.x + 2 * omt * t * cxCtrl + t * t * next.x;
                  const gy = omt * omt * curr.y + 2 * omt * t * cyCtrl + t * t * next.y;
                  connections.push(
                    <circle
                      key={`seq-gap-${chain}-${currResid}-${nextResid}`}
                      cx={gx}
                      cy={gy}
                      r={2.8}
                      fill="#0f172a"
                      opacity={0.92}
                    />
                  );
                }
              }
            });

            return connections;
          })()}

          {/* Residue nodes — Maestro teardrops (tip toward ligand) + horizontal labels */}
          {ligand_center && ligand_center.length === 2
            ? positionedResidues.map((r, i) => {
                const lcx = ligand_center[0];
                const lcy = ligand_center[1];
                const dx = lcx - r.x;
                const dy = lcy - r.y;
                const deg = (Math.atan2(dy, dx) * 180) / Math.PI - 90;

                const isGlycine = r.class === "glycine";
                const isNucleic = r.class === "nucleic";

                const strokeColor = isGlycine || isNucleic ? "#57534e" : "#0f172a";
                const strokeW = isGlycine || isNucleic ? 0.95 : 1.05;
                const gid = residueGradientId(r.class);

                const raw = (r.resname || "").toUpperCase();

                const resnameDisplay = isNucleic
                  ? raw.length <= 2
                    ? raw
                    : raw.slice(0, 2)
                  : raw.substring(0, 3);

                const residNum = r.resid ?? "";

                const nameYOffset = 20;
                const numberYOffset = 24;

                return (
                  <g key={i}>
                    <g
                      transform={`translate(${r.x},${r.y}) rotate(${deg})`}
                      style={{ filter: "drop-shadow(0 1px 2px rgba(0,0,0,0.14))" }}
                    >
                      <path
                        d={TEARDROP_D}
                        fill={`url(#${gid})`}
                        stroke={strokeColor}
                        strokeWidth={strokeW}
                        strokeLinejoin="round"
                      />
                    </g>

                    <text
                      x={r.x}
                      y={r.y - nameYOffset}
                      textAnchor="middle"
                      fontSize="9"
                      fill="#0f172a"
                      fontWeight="700"
                      style={{ textShadow: "0 1px 2px rgba(255,255,255,0.85)" }}
                    >
                      {resnameDisplay}
                    </text>

                    <text
                      x={r.x}
                      y={r.y + numberYOffset}
                      textAnchor="middle"
                      fontSize="8"
                      fill="#1e293b"
                      fontWeight="700"
                      style={{ textShadow: "0 1px 2px rgba(255,255,255,0.85)" }}
                    >
                      {residNum}
                    </text>
                  </g>
                );
              })
            : positionedResidues.map((r, i) => {
                const isGlycine = r.class === "glycine";
                const isNucleic = r.class === "nucleic";
                const gid = residueGradientId(r.class);
                const strokeColor = isGlycine || isNucleic ? "#57534e" : "#0f172a";
                const raw = (r.resname || "").toUpperCase();
                const resnameDisplay = isNucleic
                  ? raw.length <= 2
                    ? raw
                    : raw.slice(0, 2)
                  : raw.substring(0, 3);
                const residNum = r.resid ?? "";
                return (
                  <g key={i}>
                    <circle
                      cx={r.x}
                      cy={r.y}
                      r={11}
                      fill={`url(#${gid})`}
                      stroke={strokeColor}
                      strokeWidth={1.05}
                    />
                    <text x={r.x} y={r.y - 20} textAnchor="middle" fontSize="9" fill="#0f172a" fontWeight="700">
                      {resnameDisplay}
                    </text>
                    <text x={r.x} y={r.y + 24} textAnchor="middle" fontSize="8" fill="#1e293b" fontWeight="700">
                      {residNum}
                    </text>
                  </g>
                );
              })}

          {/* Solvent-exposed ligand atoms: small gray dots on top (Maestro / Nature LID) */}
          {Array.isArray(ligand_atom_xy) &&
          Array.isArray(ligand_atom_solvent_exposed) &&
          ligand_atom_solvent_exposed.length === ligand_atom_xy.length
            ? ligand_atom_xy.map((p, i) =>
                ligand_atom_solvent_exposed[i] ? (
                  <circle
                    key={`solvent-${i}`}
                    cx={p.x}
                    cy={p.y}
                    r={3.2}
                    fill="#94a3b8"
                    stroke="#64748b"
                    strokeWidth={0.35}
                    opacity={0.95}
                  />
                ) : null
              )
            : null}
        </svg>
      </div>

      {/* Conformer information */}
      {conformer_info && conformer_info.has_multiple && (
        <div className="mt-4 rounded-lg border border-sky-200 bg-sky-50 px-4 py-3 text-xs text-sky-900">
          <strong>Ligand conformers:</strong> {conformer_info.count} conformer(s) in structure
        </div>
      )}

      <div className="mt-4 rounded-lg border border-slate-200 bg-slate-50/80 px-3 py-2 text-xs text-slate-600">
        <span className="font-medium text-slate-700">Model</span>{" "}
        {data.selected_model || meta.selected_model || 1} / {data.total_models || meta.total_models || 1}
        <span className="mx-2 text-slate-300">·</span>
        <span className="font-medium text-slate-700">Ligand atoms</span>{" "}
        {data.ligand_atoms || meta.ligand_atoms || 0}
        <span className="mx-2 text-slate-300">·</span>
        <span className="font-medium text-slate-700">Protein atoms</span>{" "}
        {data.protein_atoms || meta.protein_atoms || 0}
        {meta.pocket_radius != null ? (
          <>
            <span className="mx-2 text-slate-300">·</span>
            <span className="font-medium text-slate-700">Pocket radius</span> {meta.pocket_radius} Å
          </>
        ) : null}
        {stats.interactionTotal != null && !Number.isNaN(stats.interactionTotal) ? (
          <>
            <span className="mx-2 text-slate-300">·</span>
            <span className="font-medium text-slate-700">Interactions</span>{" "}
            <span className="tabular-nums text-slate-600">
              {stats.interactionTotal} → {stats.interactionFiltered}
            </span>
            <span className="text-slate-400"> (detected → filtered)</span>
          </>
        ) : null}
      </div>

      <div
        className="mt-6 grid gap-6 border-t border-slate-200 pt-6 text-sm sm:grid-cols-2 lg:grid-cols-3"
        style={{ fontSize: 13 }}
      >
        {/* Column 1: Residue types (Maestro LID palette) */}
        <div>
          <div style={{ fontWeight: "600", marginBottom: "6px", color: "#374151" }}>Residue types</div>
          <LegendItem color="#f43f5e" label="Acidic (Asp, Glu)" />
          <LegendItem color="#6d28d9" label="Basic (Lys, Arg, His)" />
          <LegendItem color="#e5e7eb" label="Glycine" borderColor="#57534e" />
          <LegendItem color="#22c55e" label="Hydrophobic (incl. Tyr)" />
          <LegendItem color="#84cc16" label="Cysteine" />
          <LegendItem color="#9ca3af" label="Metal" />
        </div>

        {/* Column 2: Other */}
        <div>
          <div style={{ fontWeight: "600", marginBottom: "6px", color: "#374151" }}>Other</div>
          <LegendItem color="#38bdf8" label="Polar" />
          <LegendItem color="#9ca3af" label="DNA / RNA" />
          <LegendItem color="#94a3b8" label="Unspecified residue" />
          <LegendItem color="#ffffff" label="Water" borderColor="#cbd5e1" />
          <LegendItem color="#ffffff" label="Hydration site" borderColor="#9ca3af" />
          <LegendItem color="#ef4444" label="Hydration site (displaced)" />
        </div>

        {/* Column 3: Interaction types (Maestro LID–style: green π–π, red π–cation, purple metal, red–blue salt bridge, yellow halogen) */}
        <div>
          <div style={{ fontWeight: "600", marginBottom: "6px", color: "#374151" }}>Interaction types</div>
          <LegendItem
            gradientColors={["#4ade80", "#38bdf8"]}
            label="Pocket boundary (gradient + gap)"
            line
            dashed
          />
          <LegendItem color="#64748b" label="Contact (proximity)" line dashed />
          <LegendItem color="#22c55e" label="Distance" line dashed />
          <LegendItem color="#eab308" label="Halogen bond" line arrow />
          <LegendItem color="#7c3aed" label="Metal coordination" line />
          <LegendItem color="#22c55e" label="π–π stacking" line />
          <LegendItem color="#ef4444" label="π–cation" line />
          <LegendItem gradientColors={["#dc2626", "#2563eb"]} label="Salt bridge (red→blue)" line />
          <LegendItem color="#94a3b8" label="Solvent-exposed ligand atom (dot)" borderColor="#64748b" />
        </div>
      </div>
    </div>
  );
}

export default Interaction2D;