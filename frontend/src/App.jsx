import React from "react";
import { useSelector } from "react-redux";
import UploadPanel from "./components/UploadPanel";
import NglViewer from "./components/NglViewer";
import Interaction2D from "./components/Interaction2D";

function App() {
  const { data } = useSelector((state) => state.diagram);

  return (
    <div style={{ padding: 200, flexDirection: "row", display: "flex", gap: 20 }}>
      {/* <h2>Protein–Ligand Pocket Viewer</h2> */}

      <UploadPanel />

      {data && (
        <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 20 }}>
       {/* <NglViewer
  pdbUrl={data.pdb_url}
  ligandResname={data.meta.ligand_resname} // should be UNL
  pocketRadius={data.meta.pocket_radius}
/> */}
          <Interaction2D data={data} />
        </div>
      )}
    </div>
  );
}

export default App;