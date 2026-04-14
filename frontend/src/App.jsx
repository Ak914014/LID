import React from "react";
import { useSelector } from "react-redux";
import UploadPanel from "./components/UploadPanel";
import Interaction2D from "./components/Interaction2D";
import Pocket3DPlotly from "./components/Pocket3DPlotly";

function App() {
  const { data, loading, error } = useSelector((state) => state.diagram);

  return (
    <div className="min-h-screen bg-white ">
  

      <main className="mx-auto  px-4 py-8">
        <div className="flex   flex-col gap-8 lg:flex-row lg:items-start">
          <aside className="w-full shrink-0 lg:sticky lg:top-24 lg:w-80">
            <UploadPanel />
          </aside>

          <div className="flex w-[95vw] flex-col gap-10">
            {error ? (
              <div
                className="mb-4 rounded-xl border border-red-200 bg-red-50 px-4 py-3 text-sm text-red-800 shadow-sm"
                role="alert"
              >
                {error}
              </div>
            ) : null}

            {loading ? (
              <div className="flex min-h-[320px] flex-col items-center justify-center rounded-2xl border border-slate-200 bg-white p-12 shadow-sm">
                <div
                  className="mb-3 h-9 w-9 animate-spin rounded-full border-2 border-slate-200 border-t-sky-600"
                  aria-hidden
                />
                <p className="text-sm font-medium text-slate-600">Generating diagram…</p>
                <p className="mt-1 text-xs text-slate-400">Detecting interactions and layout</p>
              </div>
            ) : null}

            {!loading && data ? <Interaction2D data={data} /> : null}

            <section className="w-full" aria-labelledby="pocket-3d-heading">
              <h2
                id="pocket-3d-heading"
                className="mb-3 text-lg font-semibold tracking-tight text-slate-800"
              >
                3D pocket (Plotly demo)
              </h2>
              <p className="mb-4 max-w-2xl text-sm text-slate-600">
                Parametric binding pocket: ligand at the center, concave pocket residues (orange/red),
                extended chain (blue/green), and a translucent cavity surface. Does not use your PDB
                coordinates yet — swap in real positions from your pipeline when ready.
              </p>
              <Pocket3DPlotly />
            </section>

            {!loading && !data && !error ? (
              <div className="rounded-2xl border border-dashed border-slate-300 bg-white/80 p-12 text-center shadow-sm">
                <p className="text-base font-medium text-slate-700">No diagram yet</p>
                <p className="mt-2 text-sm text-slate-500">
                  Upload a protein–ligand PDB and click Generate to preview the binding site.
                </p>
              </div>
            ) : null}
          </div>
        </div>
      </main>
    </div>
  );
}

export default App;
