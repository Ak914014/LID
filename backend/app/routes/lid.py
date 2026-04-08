import os
from typing import Optional

from fastapi import APIRouter, UploadFile, File, Query, HTTPException
from fastapi.responses import JSONResponse

from ..core.config import settings
from ..utils.files import save_upload_bytes
from ..services.diagram_service import build_diagram, extract_selected_model_from_text

router = APIRouter(prefix="/lid", tags=["lid"])


@router.post("/generate")
async def generate_lid(
    protein_ligand_pdb: UploadFile = File(...),
    ligand_sdf: Optional[UploadFile] = File(None),
    ligand_name: Optional[str] = Query(None, description="Ligand residue name (optional)"),
    model_index: int = Query(1, ge=1, description="Pose/model index from MODEL/ENDMDL"),
    pocket_radius: float = Query(5.0),
    svg_w: int = Query(600),
    svg_h: int = Query(600),
):
    try:
        pdb_bytes = await protein_ligand_pdb.read()
        sdf_bytes = await ligand_sdf.read() if ligand_sdf is not None else None

        cleaned_for_view, _ = extract_selected_model_from_text(
            pdb_bytes.decode("utf-8"), model_index=model_index
        )
        pdb_path = save_upload_bytes(
            settings.UPLOAD_DIR,
            protein_ligand_pdb.filename,
            cleaned_for_view.encode("utf-8"),
        )

        diagram = build_diagram(
            pdb_bytes=pdb_bytes,
            ligand_resname=ligand_name,
            pocket_radius=pocket_radius,
            svg_w=svg_w,
            svg_h=svg_h,
            sdf_bytes=sdf_bytes,
            model_index=model_index,
        )
        filename = os.path.basename(pdb_path)
        diagram["pdb_url"] = f"/static/{filename}"
        return diagram
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        return JSONResponse(
            {
                "error": "Internal error",
                "detail": str(e),
            },
            status_code=500,
        )
