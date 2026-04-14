import React, { useEffect, useState } from "react";
import Plot from "react-plotly.js";
import api from "../api";

/**
 * 3D pocket demo: ligand (magenta), pocket residues (orange/red), chain (blue/green),
 * translucent pocket mesh — data from FastAPI `/api/viz3d/pocket-demo`.
 */
export default function Pocket3DPlotly() {
  const [figure, setFigure] = useState(null);
  const [err, setErr] = useState(null);

  useEffect(() => {
    let cancelled = false;
    (async () => {
      try {
        const { data } = await api.get("/api/viz3d/pocket-demo");
        if (!cancelled) setFigure(data);
      } catch (e) {
        if (!cancelled) setErr(e?.message ?? "Failed to load 3D plot");
      }
    })();
    return () => {
      cancelled = true;
    };
  }, []);

  if (err) {
    return (
      <div
        className="rounded-2xl border border-amber-200 bg-amber-50 px-4 py-3 text-sm text-amber-900"
        role="status"
      >
        <p className="font-medium">3D demo unavailable</p>
        <p className="mt-1 text-amber-800/90">{err}</p>
        <p className="mt-2 text-xs text-amber-800/70">
          Start the API (e.g. uvicorn on port 8008) so Plotly JSON can load.
        </p>
      </div>
    );
  }

  if (!figure) {
    return (
      <div className="flex min-h-[320px] flex-col items-center justify-center rounded-2xl border border-slate-200 bg-slate-50/80">
        <div
          className="mb-3 h-9 w-9 animate-spin rounded-full border-2 border-slate-200 border-t-violet-600"
          aria-hidden
        />
        <p className="text-sm font-medium text-slate-600">Loading 3D pocket…</p>
      </div>
    );
  }

  return (
    <div className="w-full overflow-hidden rounded-2xl border border-slate-200 bg-white shadow-sm">
      <Plot
        data={figure.data}
        layout={figure.layout}
        config={{
          responsive: true,
          displayModeBar: true,
          scrollZoom: true,
        }}
        style={{ width: "100%", minHeight: 520 }}
        useResizeHandler
      />
    </div>
  );
}
