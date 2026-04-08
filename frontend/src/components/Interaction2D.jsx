import React from "react";

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
    backbone_path,
    meta = {},
    interactions = [],
    ligand_atom_xy,
    ligand_center,
    conformer_info,
  } = safeData;

  if (!data) return null;

  // Server layout + pocket share the same geometry; do not force a client-side circle.
  const positionedResidues = residues;

  const residueIndex = (() => {
    const map = new Map();
    positionedResidues.forEach((r) => map.set(`${r.resname}${r.resid}`, r));
    return map;
  })();

  // Maestro LID–style: orange acidic, purple basic, green hydrophobic, cyan polar,
  // grey gly/nucleic, yellow–green cysteine (matches Schrödinger palette)
  const residueFill = (cls) => {
    if (cls === "negative") return "#ea5806";     // orange – Acidic (Asp, Glu)
    if (cls === "positive") return "#9333ea";   // purple – Basic (Lys, Arg, His)
    if (cls === "hydrophobic") return "#16a34a"; // green – Hydrophobic
    if (cls === "polar") return "#06b6d4";       // cyan – Polar
    if (cls === "cysteine") return "#84cc16";   // yellow-green – Cys
    if (cls === "glycine") return "#d1d5db";    // light grey – Gly (Maestro)
    if (cls === "metal") return "#9ca3af";
    if (cls === "nucleic") return "#9ca3af";    // grey – DNA/RNA bases
    return "#94a3b8";
  };

  /** Teardrop points along +Y in local coords; rotate so tip aims toward ligand center. */
  const teardropRotationDeg = (rx, ry, lcx, lcy) => {
    const dx = lcx - rx;
    const dy = lcy - ry;
    return (Math.atan2(dy, dx) * 180) / Math.PI - 90;
  };

  // Local SVG path: tip toward +Y (down), bulb toward −Y — Maestro-like droplet
  const TEARDROP_PATH =
    "M 0 -13 Q 15 2 0 17 Q -15 2 0 -13 Z";

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
          label: isBackbone ? "backbone" : "sidechain"
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
      style={{
        margin: "30px auto",
        padding: 20,
        borderRadius: 14,
        width: meta.svg_w + 36,
      }}
    >
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
            {/* Pocket boundary: green → amber → cyan (Maestro LID “glow” band) */}
            <linearGradient id="pocketBandGrad" x1="0%" y1="0%" x2="100%" y2="100%">
              <stop offset="0%" stopColor="#4ade80" />
              <stop offset="45%" stopColor="#facc15" />
              <stop offset="100%" stopColor="#38bdf8" />
            </linearGradient>
          </defs>

          {/* Pocket: hugs ligand; gaps = solvent-exposed (no line) */}
          {pocket_outline_path ? (
            <g style={{ filter: "drop-shadow(0 0 4px rgba(56, 189, 248, 0.35))" }}>
              <path
                d={pocket_outline_path}
                fill="none"
                stroke="url(#pocketBandGrad)"
                strokeWidth={6}
                opacity={0.92}
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </g>
          ) : null}

          {/* Smooth protein backbone (B-spline through consecutive residues — matches pocket geometry from API) */}
          {backbone_path && String(backbone_path).trim().length > 0 ? (
            <path
              d={backbone_path}
              fill="none"
              stroke="#0f172a"
              strokeWidth={2.2}
              opacity={0.9}
              strokeLinecap="round"
              strokeLinejoin="round"
              style={{ filter: "drop-shadow(0 1px 1px rgba(0,0,0,0.08))" }}
            />
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
                  {it.type === "hbond" && st.label && (
                    <text
                      x={midX}
                      y={midY - 8}
                      textAnchor="middle"
                      fontSize="8"
                      fill="#c026d3"
                      fontWeight="600"
                      opacity={0.9}
                      style={{ textShadow: "0 1px 2px rgba(255,255,255,0.8)" }}
                    >
                      {st.label}
                    </text>
                  )}
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

                if (nextResid - currResid !== 1) continue;

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
              }
            });

            return connections;
          })()}

          {/* Residue nodes — Maestro-style teardrops (tip points toward ligand) */}
          {positionedResidues.map((r, i) => {


  const isGlycine = r.class === "glycine";
  const isNucleic = r.class === "nucleic";

  const fillColor = residueFill(r.class);
  const strokeColor = isGlycine || isNucleic ? "#57534e" : "#1e293b";
  const strokeW = isGlycine || isNucleic ? 0.2 : 0.1;

  const raw = (r.resname || "").toUpperCase();

  const resnameDisplay = isNucleic
    ? raw.length <= 2
      ? raw
      : raw.slice(0, 2)
    : raw.substring(0, 3);

  const residNum = r.resid ?? "";

  const scoreRaw = r.strain_score;
  const scoreStr =
    scoreRaw != null && !Number.isNaN(Number(scoreRaw))
      ? `${Number(scoreRaw) >= 0 ? "+" : ""}${Number(scoreRaw).toFixed(3)}`
      : "+0.000";

  return (
    <g key={i}>
      {/* 🔵 Circle Node */}
      <g
        transform={`translate(${r.x},${r.y}) scale(3)`}
        style={{ filter: "drop-shadow(0 1px 2px rgba(0,0,0,0.12))" }}
      >
        <circle
          cx={0}
          cy={0}
          r={8} // 👈 base radius (adjust this if needed)
          fill={fillColor}
          stroke={strokeColor}
          strokeWidth={strokeW}
        />
      </g>

      {/* 🔤 Residue Name */}
      <text
        x={r.x}
        y={r.y - 8}
        textAnchor="middle"
        fontSize="9"
        fill="#0f172a"
        fontWeight="700"
        style={{ textShadow: "0 1px 2px rgba(255,255,255,0.85)" }}
      >
        {resnameDisplay}
      </text>

      {/* 📊 Score */}
      <text
        x={r.x}
        y={r.y + 2}
        textAnchor="middle"
        fontSize="7.5"
        fill="#64748b"
        fontWeight="600"
        style={{ textShadow: "0 1px 2px rgba(255,255,255,0.85)" }}
      >
        {scoreStr}
      </text>

      {/* 🔢 Residue Number */}
      <text
        x={r.x}
        y={r.y + 12}
        textAnchor="middle"
        fontSize="8"
        fill="#334155"
        fontWeight="700"
        style={{ textShadow: "0 1px 2px rgba(255,255,255,0.85)" }}
      >
        {residNum}
      </text>
    </g>
  );
})}
        </svg>
      </div>

      {/* Conformer information */}
      {conformer_info && conformer_info.has_multiple && (
        <div
          style={{
            marginTop: 12,
            padding: "20px 20px",
            background: "#f0f9ff",
            border: "1px solid #bae6fd",
            borderRadius: 6,
            fontSize: 12,
            color: "#0369a1",
          }}
        >
          <strong>Ligand Conformers:</strong> {conformer_info.count} conformer(s) detected in structure
        </div>
      )}

      <div
        style={{
          marginTop: 150,
          padding: "8px 12px",
          background: "#f8fafc",
          border: "1px solid #e2e8f0",
          borderRadius: 6,
          fontSize: 12,
          color: "#334155",
        }}
      >
        {`Selected model: ${data.selected_model || meta.selected_model || 1} / ${data.total_models || meta.total_models || 1} | `}
        {`Ligand atoms: ${data.ligand_atoms || meta.ligand_atoms || 0} | `}
        {`Protein atoms: ${data.protein_atoms || meta.protein_atoms || 0}`}
      </div>

      {/* Legend matching reference image - 3 columns */}
      <div
        style={{
          marginTop: 140,
          paddingTop: 12,
          borderTop: "1px solid #e5e7eb",
          display: "grid",
          gridTemplateColumns: "repeat(3, minmax(200px, 1fr))",
          gap: "8px 14px",
          fontSize: 13,
        }}
      >
        {/* Column 1: Residue types (Maestro LID: Red=acidic, Purple=basic, Green=hydrophobic, Blue=polar) */}
        <div>
          <div style={{ fontWeight: "600", marginBottom: "6px", color: "#374151" }}>Residue types</div>
          <LegendItem color="#ea5806" label="Acidic (Asp, Glu)" />
          <LegendItem color="#9333ea" label="Basic (Lys, Arg, His)" />
          <LegendItem color="#d1d5db" label="Glycine" borderColor="#57534e" />
          <LegendItem color="#16a34a" label="Hydrophobic" />
          <LegendItem color="#84cc16" label="Cysteine" />
          <LegendItem color="#9ca3af" label="Metal" />
        </div>

        {/* Column 2: Other */}
        <div>
          <div style={{ fontWeight: "600", marginBottom: "6px", color: "#374151" }}>Other</div>
          <LegendItem color="#06b6d4" label="Polar" />
          <LegendItem color="#9ca3af" label="DNA / RNA" />
          <LegendItem color="#94a3b8" label="Unspecified residue" />
          <LegendItem color="#ffffff" label="Water" borderColor="#cbd5e1" />
          <LegendItem color="#ffffff" label="Hydration site" borderColor="#9ca3af" />
          <LegendItem color="#ef4444" label="Hydration site (displaced)" />
        </div>

        {/* Column 3: Interaction types (Maestro LID–style: green π–π, red π–cation, purple metal, red–blue salt bridge, yellow halogen) */}
        <div>
          <div style={{ fontWeight: "600", marginBottom: "6px", color: "#374151" }}>Interaction types</div>
          <LegendItem color="#22c55e" label="Distance" line dashed />
          <LegendItem color="#eab308" label="Halogen bond" line arrow />
          <LegendItem color="#7c3aed" label="Metal coordination" line />
          <LegendItem color="#22c55e" label="π–π stacking" line />
          <LegendItem color="#ef4444" label="π–cation" line />
          <LegendItem gradientColors={["#dc2626", "#2563eb"]} label="Salt bridge (red→blue)" line />
          <LegendItem color="#9ca3af" label="Solvent exposure" borderColor="#9ca3af" />
        </div>
      </div>
    </div>
  );
}

export default Interaction2D;