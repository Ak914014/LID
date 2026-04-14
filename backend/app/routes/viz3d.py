from fastapi import APIRouter

from ..services.pocket_viz_3d import build_pocket_demo_figure, figure_to_plotly_json

router = APIRouter(prefix="/viz3d", tags=["viz3d"])


@router.get("/pocket-demo")
def pocket_demo_plotly():
    """Plotly figure JSON for react-plotly (data + layout)."""
    fig = build_pocket_demo_figure()
    return figure_to_plotly_json(fig)
