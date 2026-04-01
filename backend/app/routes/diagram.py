import os
from typing import Optional

from fastapi import APIRouter, UploadFile, File, Query, HTTPException
from fastapi.responses import JSONResponse

from ..core.config import settings
from ..utils.files import save_upload_bytes
from ..services.diagram_service import build_diagram, extract_selected_model_from_text

router = APIRouter(prefix="/api/diagram", tags=["diagram"])


@router.post("/generate")
@router.post("/lid/generate")
async def generate_diagram(
    protein_ligand_pdb: UploadFile = File(...),
    ligand_sdf: Optional[UploadFile] = File(None),
    ligand_resname: Optional[str] = Query(None, description="Ligand residue name (auto-detected if not provided)"),
    ligand_name: Optional[str] = Query(None, description="Alias of ligand_resname"),
    model_index: int = Query(1, ge=1, description="Pose/model index from MODEL/ENDMDL"),
    pocket_radius: float = Query(5.0),
    svg_w: int = Query(820),
    svg_h: int = Query(520),
):
    try:
        # Read uploaded files
        pdb_bytes = await protein_ligand_pdb.read()

        sdf_bytes = None
        if ligand_sdf is not None:
            sdf_bytes = await ligand_sdf.read()

        # Save the same selected pose as the 2D diagram for NGL (receptor + ligand aligned).
        cleaned_for_view, _ = extract_selected_model_from_text(
            pdb_bytes.decode("utf-8"), model_index=model_index
        )
        pdb_path = save_upload_bytes(
            settings.UPLOAD_DIR,
            protein_ligand_pdb.filename,
            cleaned_for_view.encode("utf-8")
        )

        resolved_ligand_name = ligand_name or ligand_resname
        diagram = build_diagram(
            pdb_bytes=pdb_bytes,
            ligand_resname=resolved_ligand_name,
            pocket_radius=pocket_radius,
            svg_w=svg_w,
            svg_h=svg_h,
            sdf_bytes=sdf_bytes,
            model_index=model_index,
        )

        # Build static URL (no hardcoding port if possible)
        filename = os.path.basename(pdb_path)
        public_pdb_url = f"/static/{filename}"
        diagram["pdb_url"] = public_pdb_url

        return diagram

    except ValueError as e:
        import sys
        error_msg = str(e)
        print(f"ValueError in generate_diagram: {error_msg}", file=sys.stderr)
        raise HTTPException(status_code=400, detail=error_msg)

    except Exception as e:
        import traceback
        import sys
        exc_type, exc_value, exc_tb = sys.exc_info()
        error_detail = str(e)
        traceback_str = traceback.format_exc()
        
        # Log the error for debugging
        print(f"ERROR in generate_diagram: {error_detail}", file=sys.stderr)
        print(f"Traceback: {traceback_str}", file=sys.stderr)
        
        return JSONResponse(
            {
                "error": "Internal error",
                "detail": error_detail,
                "type": str(exc_type.__name__) if exc_type else "Unknown",
                "traceback": traceback_str,
            },
            status_code=500,
        )