import React, { useEffect, useRef } from "react";
import * as NGL from "ngl";

function NglViewer({ pdbUrl, ligandResname = "UNL", pocketRadius = 5.0 }) {
  const ref = useRef(null);

  useEffect(() => {
    if (!ref.current || !pdbUrl) return;

    const stage = new NGL.Stage(ref.current, { backgroundColor: "white" });

    stage.loadFile(pdbUrl).then((comp) => {
      // Protein cartoon
      comp.addRepresentation("cartoon", { sele: "protein", opacity: 1.0 });
      comp.addRepresentation("ball+stick", { sele: `resname ${ligandResname}` });

      // Pocket residues within radius
      const pocketSele = `protein and within ${pocketRadius} of (resname ${ligandResname})`;

      // Pocket sticks (so you see residues around ligand)
      comp.addRepresentation("ball+stick", {
        sele: pocketSele,
        opacity: 1.0,
      });

      // Pocket surface (the “real pocket” look)
      comp.addRepresentation("surface", {
        sele: pocketSele,
        opacity: 0.35,
      });

      comp.autoView(`resname ${ligandResname}`);
    });

    const handleResize = () => stage.handleResize();
    window.addEventListener("resize", handleResize);

    return () => {
      window.removeEventListener("resize", handleResize);
      stage.dispose();
    };
  }, [pdbUrl, ligandResname, pocketRadius]);

  return <div ref={ref} style={{ width: "100%", height: 520, borderRadius: 12, overflow: "hidden" }} />;
}

export default NglViewer;