import React, { Suspense } from "react";
import { useSelector } from "react-redux";
import UploadPanel from "./components/UploadPanel";
import Interaction2D from "./components/Interaction2D";

// const MolstarViewer = React.lazy(() => import("./components/MolstarViewer"));

function App() {
  const { data } = useSelector((state) => state.diagram);

  return (
    <div
     className="p-4 flex  gap-"
    >
      <UploadPanel />

      
        <div
          
        >
          {data && (
          <Interaction2D data={data} />
          
        )}
        </div>
      
    </div>
  );
}

export default App;