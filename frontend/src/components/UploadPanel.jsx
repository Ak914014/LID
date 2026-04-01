import React, { useState } from "react";
import { useDispatch, useSelector } from "react-redux";
import { generateDiagram } from "../features/diagramSlice";

function UploadPanel() {
  const dispatch = useDispatch();
  const { loading } = useSelector((state) => state.diagram);

  const [pdbFile, setPdbFile] = useState(null);
  const [ligandResname, setLigandResname] = useState("");
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

    // Ligand resname is optional - will be auto-detected if not provided

    dispatch(
      generateDiagram({
        pdbFile,
        sdfFile: null,  // Always null - ligand extracted from PDB
        ligandResname,
        pocketRadius,
        modelIndex,
      })
    );
  };

  return (
    <div style={{ marginBottom: 20 }}>


      <div>
        <label><strong>Complex PDB File (Protein + Ligand)</strong></label>
        <br />
        <input
          type="file"
          accept=".pdb,.ent,.pdbqt"
          onChange={(e) => handlePdbChange(e.target.files[0])}
          style={{ marginTop: "5px" }}
        />
      </div>

      <br />

      {/* <div>
        <label><strong>Ligand Residue Name</strong> (optional - will auto-detect if not provided)</label>
        <br />
        <input
          value={ligandResname}
          onChange={(e) => setLigandResname(e.target.value)}
          placeholder="e.g. UNL, LIG, DRG (leave empty for auto-detection)"
          style={{ marginTop: "5px", padding: "5px", width: "300px" }}
        />
        <p style={{ fontSize: "12px", color: "#666", marginTop: "3px" }}>
          Leave empty to auto-detect, or specify: UNL, LIG, DRG, or your custom ligand name
        </p>
      </div> */}

      <br />

      <div>
        <label>Select Pose</label>
        <br />
        <select
          value={modelIndex}
          onChange={(e) => setModelIndex(Number(e.target.value))}
          style={{ marginTop: "5px", padding: "5px" }}
        >
          {Array.from({ length: modelCount }, (_, i) => i + 1).map((n) => (
            <option key={n} value={n}>
              {`Model ${n}`}
            </option>
          ))}
        </select>
      </div>

      <br />

      <div>
        <label>Pocket Radius (Å)</label>
        <br />
        <input
          type="number"
          value={pocketRadius}
          onChange={(e) => setPocketRadius(Number(e.target.value))}
          min="3"
          max="10"
        />
      </div>

      <br />

      <button onClick={handleSubmit} disabled={loading}>
        {loading ? "Generating..." : "Generate Diagram"}
      </button>
    </div>
  );
}

export default UploadPanel;