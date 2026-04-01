"""
SVG Renderer
Generates SVG output for interaction diagrams.
"""

from typing import List, Dict, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import numpy as np


def render_ligand_svg(mol, width: int = 820, height: int = 520) -> str:
    """
    Render ligand as clean 2D SVG.
    
    Args:
        mol: RDKit molecule
        width: SVG width
        height: SVG height
    
    Returns:
        SVG string
    """
    # Remove hydrogens
    mol_no_h = Chem.RemoveHs(mol)
    
    # Compute 2D coordinates
    if not mol_no_h.GetNumConformers():
        AllChem.Compute2DCoords(mol_no_h)
    
    # Create drawer
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    
    opts.addStereoAnnotation = True
    opts.addAtomIndices = False
    opts.addBondIndices = False
    opts.includeMetadata = False
    opts.kekulize = False  # Keep double bonds visible
    
    drawer.DrawMolecule(mol_no_h)
    drawer.FinishDrawing()
    
    return drawer.GetDrawingText()


def create_quadratic_bezier_path(
    x1: float, y1: float, 
    x2: float, y2: float, 
    curvature: float = 0.3,
    avoid_center: bool = True,
    center_threshold: float = 30.0
) -> str:
    """
    Create quadratic Bézier path between two points.
    
    Args:
        x1, y1: Start point
        x2, y2: End point
        curvature: Curvature factor (0.0 = straight, higher = more curved)
        avoid_center: If True, curve outward if line passes near origin
        center_threshold: Distance threshold for center avoidance
    """
    dx = x2 - x1
    dy = y2 - y1
    dist = np.sqrt(dx * dx + dy * dy)
    
    mid_x = (x1 + x2) / 2
    mid_y = (y1 + y2) / 2
    
    # Check if line passes too close to center
    if avoid_center:
        # Distance from origin to line segment
        # Using point-to-line distance formula
        a = dy
        b = -dx
        c = x1 * dy - y1 * dx
        dist_to_center = abs(c) / (np.sqrt(a * a + b * b) + 1e-9)
        
        if dist_to_center < center_threshold:
            # Increase curvature to avoid center
            curvature = max(curvature, 0.5)
    
    angle = np.arctan2(dy, dx)
    perp_angle = angle + np.pi / 2
    
    # Perpendicular offset direction
    offset = dist * curvature
    n_x = -dy / (dist + 1e-9)
    n_y = dx / (dist + 1e-9)
    
    cp_x = mid_x + offset * n_x
    cp_y = mid_y + offset * n_y
    
    return f"M {x1} {y1} Q {cp_x} {cp_y}, {x2} {y2}"


def create_interaction_edge_path(
    residue_pos: np.ndarray,
    ligand_pos: np.ndarray,
    interaction_type: str
) -> str:
    """
    Create edge path for interaction based on type.
    
    Args:
        residue_pos: 2D position of residue node [x, y]
        ligand_pos: 2D position of ligand anchor [x, y]
        interaction_type: Type of interaction (hbond, salt_bridge, pi_pi, etc.)
    
    Returns:
        SVG path string
    """
    x1, y1 = float(residue_pos[0]), float(residue_pos[1])
    x2, y2 = float(ligand_pos[0]), float(ligand_pos[1])
    
    if interaction_type == "hbond":
        # Curved quadratic Bézier for H-bonds
        return create_quadratic_bezier_path(x1, y1, x2, y2, curvature=0.4, avoid_center=True)
    elif interaction_type in ["pi_pi", "pi_cation"]:
        # Curved for pi interactions
        return create_quadratic_bezier_path(x1, y1, x2, y2, curvature=0.35, avoid_center=True)
    elif interaction_type == "salt_bridge":
        # Straight line for salt bridge
        return f"M {x1} {y1} L {x2} {y2}"
    elif interaction_type == "hydrophobic":
        # Slightly curved for hydrophobic
        return create_quadratic_bezier_path(x1, y1, x2, y2, curvature=0.2, avoid_center=True)
    else:
        # Default: slightly curved
        return create_quadratic_bezier_path(x1, y1, x2, y2, curvature=0.25, avoid_center=True)


def create_circular_path(
    residues: List[Dict],
    ligand_center: tuple,
    is_inner: bool = False
) -> Optional[str]:
    """Create circular path connecting consecutive residues."""
    if len(residues) < 2:
        return None
    
    path_points = []
    lig_cx, lig_cy = ligand_center
    
    for res in residues:
        angle = np.arctan2(res["y"] - lig_cy, res["x"] - lig_cx)
        path_points.append({
            "x": res["x"],
            "y": res["y"],
            "angle": angle,
        })
    
    path_points.sort(key=lambda p: p["angle"])
    
    curve_dir = -1 if is_inner else 1
    base_offset = 25 if is_inner else 35
    
    path_data = f"M {path_points[0]['x']} {path_points[0]['y']}"
    
    for i in range(len(path_points)):
        current = path_points[i]
        next_p = path_points[(i + 1) % len(path_points)]
        
        dx = next_p["x"] - current["x"]
        dy = next_p["y"] - current["y"]
        dist = np.sqrt(dx * dx + dy * dy)
        offset = min(dist * 0.4, base_offset) * curve_dir
        
        angle = np.arctan2(dy, dx)
        perp_angle = angle + np.pi / 2
        
        cp1x = current["x"] + np.cos(perp_angle) * offset
        cp1y = current["y"] + np.sin(perp_angle) * offset
        cp2x = next_p["x"] + np.cos(perp_angle) * offset
        cp2y = next_p["y"] + np.sin(perp_angle) * offset
        
        if i == len(path_points) - 1 and len(path_points) > 2:
            # Close circle
            first = path_points[0]
            dx_close = first["x"] - current["x"]
            dy_close = first["y"] - current["y"]
            angle_close = np.arctan2(dy_close, dx_close)
            perp_angle_close = angle_close + np.pi / 2
            offset_close = min(np.sqrt(dx_close * dx_close + dy_close * dy_close) * 0.4, base_offset) * curve_dir
            
            cp1_close_x = current["x"] + np.cos(perp_angle_close) * offset_close
            cp1_close_y = current["y"] + np.sin(perp_angle_close) * offset_close
            cp2_close_x = first["x"] + np.cos(perp_angle_close) * offset_close
            cp2_close_y = first["y"] + np.sin(perp_angle_close) * offset_close
            
            path_data += f" C {cp1_close_x} {cp1_close_y}, {cp2_close_x} {cp2_close_y}, {first['x']} {first['y']} Z"
        else:
            path_data += f" C {cp1x} {cp1y}, {cp2x} {cp2y}, {next_p['x']} {next_p['y']}"
    
    return path_data
