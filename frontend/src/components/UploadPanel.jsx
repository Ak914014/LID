import React, { useState } from "react";
import { useDispatch, useSelector } from "react-redux";
import { generateDiagram } from "../features/diagramSlice";

function UploadPanel() {
  const dispatch = useDispatch();
  const { loading } = useSelector((state) => state.diagram);

  const [pdbFile, setPdbFile] = useState(null);
  const [pocketRadius, setPocketRadius] = useState(5);
  const [modelCount, setModelCount] = useState(1);
  const [modelIndex, setModelIndex] = useState(1);

  const handlePdbChange = async (file) => {
    setPdbFile(file);
    if (!file) {
      setModelCount(1);
      setModelIndex(1);
      return;
    }
    try {
      const text = await file.text();
      const matches = text.match(/^MODEL\b.*$/gm);
      const detected = matches ? matches.length : 0;
      const count = detected > 0 ? detected : 1;
      setModelCount(count);
      setModelIndex((prev) => Math.min(Math.max(prev, 1), count));
    } catch {
      setModelCount(1);
      setModelIndex(1);
    }
  };

  const handleSubmit = () => {
    if (!pdbFile) {
      alert("Please upload the complex PDB file");
      return;
    }

    dispatch(
      generateDiagram({
        pdbFile,
        sdfFile: null,
        ligandResname: "",
        pocketRadius,
        modelIndex,
      })
    );
  };

  return (
    <div className="rounded-2xl border border-slate-200/90 bg-white p-6 shadow-sm ring-1 ring-slate-900/5">
      <h2 className="text-sm font-semibold uppercase tracking-wide text-slate-500">
        Input
      </h2>
      <p className="mt-1 text-sm text-slate-600">
        Complex PDB (protein + ligand). Ligand is detected automatically.
      </p>

      <div className="mt-5 space-y-5">
        <div>
          <label className="text-sm font-medium text-slate-800" htmlFor="pdb-file">
            Complex PDB
          </label>
          <input
            id="pdb-file"
            type="file"
            accept=".pdb,.ent,.pdbqt"
            onChange={(e) => handlePdbChange(e.target.files[0])}
            className="mt-2 block w-full cursor-pointer rounded-lg border border-slate-200 bg-slate-50/80 px-3 py-2 text-sm file:mr-3 file:rounded-md file:border-0 file:bg-sky-600 file:px-3 file:py-1.5 file:text-sm file:font-medium file:text-white hover:file:bg-sky-700"
          />
          {pdbFile ? (
            <p className="mt-2 truncate text-xs text-slate-500" title={pdbFile.name}>
              Selected: <span className="font-medium text-slate-700">{pdbFile.name}</span>
            </p>
          ) : null}
        </div>

        <div>
          <label className="text-sm font-medium text-slate-800" htmlFor="pose">
            Pose / model
          </label>
          <select
            id="pose"
            value={modelIndex}
            onChange={(e) => setModelIndex(Number(e.target.value))}
            className="mt-2 w-full rounded-lg border border-slate-200 bg-white px-3 py-2 text-sm text-slate-800 shadow-sm focus:border-sky-500 focus:outline-none focus:ring-2 focus:ring-sky-500/20"
          >
            {Array.from({ length: modelCount }, (_, i) => i + 1).map((n) => (
              <option key={n} value={n}>
                Model {n}
              </option>
            ))}
          </select>
        </div>

        <div>
          <label className="text-sm font-medium text-slate-800" htmlFor="radius">
            Pocket radius (Å)
          </label>
          <div className="mt-2 flex items-center gap-3">
            <input
              id="radius"
              type="range"
              min={3}
              max={10}
              step={0.5}
              value={pocketRadius}
              onChange={(e) => setPocketRadius(Number(e.target.value))}
              className="h-2 flex-1 cursor-pointer accent-sky-600"
            />
            <span className="w-10 tabular-nums text-right text-sm font-medium text-slate-700">
              {pocketRadius}
            </span>
          </div>
          <p className="mt-1 text-xs text-slate-500">
            Controls how many protein/nucleic residues are included around the ligand (distance shell).
          </p>
        </div>

        <button
          type="button"
          onClick={handleSubmit}
          disabled={loading}
          className="w-full rounded-xl bg-sky-600 px-4 py-3 text-sm font-semibold text-white shadow-sm transition hover:bg-sky-700 focus:outline-none focus:ring-2 focus:ring-sky-500 focus:ring-offset-2 disabled:hover:bg-sky-600"
        >
          {loading ? "Generating…" : "Generate diagram"}
        </button>
      </div>
    </div>
  );
}

export default UploadPanel;
