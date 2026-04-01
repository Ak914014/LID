"""
Interaction Anchor Mapper
Maps each residue interaction to its specific ligand atom anchor position.
"""

from typing import List, Dict, Tuple
import numpy as np


def map_interactions_to_ligand_atoms(
    interactions: List[Dict],
    ligand_atom_xy: List[Dict],
    ligand_center: Tuple[float, float]
) -> Dict[str, Dict]:
    """
    Map each residue interaction to its ligand atom anchor.
    
    Args:
        interactions: List of interaction dicts with "ligand_atom_index"
        ligand_atom_xy: List of ligand atom positions [{"x": float, "y": float}, ...]
        ligand_center: (x, y) coordinates of ligand center
    
    Returns:
        Dict mapping residue_id -> {
            "ligand_atom_index": int,
            "ligand_atom_pos": np.array([x, y]),
            "angle": float,  # Angle of anchor relative to center
            "interactions": [list of interaction dicts]
        }
    """
    residue_anchors = {}
    lig_cx, lig_cy = ligand_center
    
    for it in interactions:
        res_id = f"{it['resname']}{it['resid']}"
        lig_atom_idx = it.get("ligand_atom_index", -1)
        
        # Get ligand atom position
        if 0 <= lig_atom_idx < len(ligand_atom_xy):
            lig_atom = ligand_atom_xy[lig_atom_idx]
            lig_atom_pos = np.array([float(lig_atom["x"]), float(lig_atom["y"])])
        else:
            # Fallback to ligand center
            lig_atom_pos = np.array([lig_cx, lig_cy])
            lig_atom_idx = -1
        
        # Compute angle of anchor relative to ligand center
        dx = lig_atom_pos[0] - lig_cx
        dy = lig_atom_pos[1] - lig_cy
        angle = np.arctan2(dy, dx)
        
        # Convert numpy array to list for JSON serialization
        lig_atom_pos_list = [float(lig_atom_pos[0]), float(lig_atom_pos[1])]
        
        if res_id not in residue_anchors:
            residue_anchors[res_id] = {
                "ligand_atom_index": int(lig_atom_idx),
                "ligand_atom_pos": lig_atom_pos_list,
                "angle": float(angle),
                "interactions": [],
            }
        
        residue_anchors[res_id]["interactions"].append(it)
    
    return residue_anchors


# compute_sector_angles is now implemented in sector_layout.py with adaptive spans
# This function is kept for backward compatibility but delegates to the new implementation
def compute_sector_angles(
    residue_anchors: Dict[str, Dict],
    ring_radius_map: Dict[int, float] = None
) -> List[Dict]:
    """
    Compute adaptive angular sectors for each unique ligand atom anchor.
    
    This function delegates to the implementation in sector_layout.py.
    """
    from .sector_layout import compute_sector_angles as compute_adaptive_sectors
    
    if ring_radius_map is None:
        from .sector_layout import RING_RADIUS
        ring_radius_map = RING_RADIUS
    
    return compute_adaptive_sectors(residue_anchors, ring_radius_map)


def resolve_sector_overlaps(
    sectors: List[Dict],
    min_sector_width: float
) -> List[Dict]:
    """
    Resolve overlapping sectors by adjusting boundaries.
    """
    if len(sectors) <= 1:
        return sectors
    
    # Normalize angles to [0, 2π]
    for s in sectors:
        s["start_angle"] = ((s["start_angle"] % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
        s["end_angle"] = ((s["end_angle"] % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
        s["center_angle"] = ((s["center_angle"] % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
    
    sectors.sort(key=lambda s: s["center_angle"])
    
    # Check for overlaps and adjust
    for i in range(len(sectors) - 1):
        s1 = sectors[i]
        s2 = sectors[i + 1]
        
        # Check if sectors overlap
        if s1["end_angle"] > s2["start_angle"]:
            # Overlap detected - split the gap
            gap = s2["start_angle"] - s1["end_angle"]
            if gap < 0:
                # Wrap around case
                gap = (s2["start_angle"] + 2 * np.pi) - s1["end_angle"]
            
            # Adjust boundaries to create gap
            midpoint = (s1["end_angle"] + s2["start_angle"]) / 2
            if midpoint > s1["end_angle"]:
                s1["end_angle"] = midpoint - min_sector_width / 4
                s2["start_angle"] = midpoint + min_sector_width / 4
            else:
                # Wrap around case
                s1["end_angle"] = ((midpoint - min_sector_width / 4) % (2 * np.pi) + 2 * np.pi) % (2 * np.pi)
                s2["start_angle"] = ((midpoint + min_sector_width / 4) % (2 * np.pi) + 2 * np.pi) % (2 * np.pi)
    
    return sectors
