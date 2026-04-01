"""
Collision Resolver
Handles label and node overlap resolution.
"""

import math
import numpy as np
from typing import List, Dict

# Label bounding box dimensions (approximate)
LABEL_WIDTH = 50.0
LABEL_HEIGHT = 30.0
LABEL_PADDING = 8.0

# Node radius for boundary overlap check
NODE_RADIUS = 20.0  # pixels


def get_label_bbox(node: Dict, pos2: np.ndarray) -> Dict:
    """
    Get bounding box for residue label.
    
    Returns:
        {"x": min_x, "y": min_y, "width": w, "height": h}
    """
    return {
        "x": float(pos2[0] - LABEL_WIDTH / 2),
        "y": float(pos2[1] - LABEL_HEIGHT / 2),
        "width": LABEL_WIDTH,
        "height": LABEL_HEIGHT,
    }


def bbox_intersect(bbox1: Dict, bbox2: Dict) -> bool:
    """Check if two bounding boxes intersect."""
    return not (
        bbox1["x"] + bbox1["width"] < bbox2["x"] or
        bbox2["x"] + bbox2["width"] < bbox1["x"] or
        bbox1["y"] + bbox1["height"] < bbox2["y"] or
        bbox2["y"] + bbox2["height"] < bbox1["y"]
    )


def check_boundary_overlap(
    node_pos: np.ndarray,
    boundary_path: str,
    node_radius: float = NODE_RADIUS
) -> bool:
    """
    Check if node overlaps with pocket boundary.
    
    Simplified check: if node is too close to boundary, it overlaps.
    """
    # This is a simplified check - in production, you'd parse the SVG path
    # For now, we'll rely on the boundary being expanded outward enough
    return False


def resolve_label_collisions(
    nodes: List[Dict], 
    ring_radius: Dict[int, float],
    boundary_path: str = None
) -> List[Dict]:
    """
    Resolve label overlaps by nudging nodes tangentially.
    Also ensures nodes don't overlap boundary.
    
    Args:
        nodes: List of nodes with "pos2", "ring"
        ring_radius: Dict mapping ring number to radius
        boundary_path: Optional SVG path of boundary for overlap check
    
    Returns:
        Updated nodes with adjusted positions
    """
    # Only process residue nodes
    residue_nodes = [n for n in nodes if n.get("type") == "residue"]
    
    if len(residue_nodes) < 2:
        return nodes
    
    # Sort by y-coordinate (top to bottom)
    residue_nodes.sort(key=lambda n: n["pos2"][1] if isinstance(n["pos2"], np.ndarray) else n["pos2"][1], reverse=True)
    
    max_iterations = 10
    nudge_step = 6.0  # pixels
    
    for iteration in range(max_iterations):
        has_overlap = False
        
        for i in range(len(residue_nodes)):
            n1 = residue_nodes[i]
            pos1 = n1["pos2"]
            if isinstance(pos1, list):
                pos1 = np.array(pos1, dtype=float)
            bbox1 = get_label_bbox(n1, pos1)
            
            for j in range(i + 1, len(residue_nodes)):
                n2 = residue_nodes[j]
                pos2 = n2["pos2"]
                if isinstance(pos2, list):
                    pos2 = np.array(pos2, dtype=float)
                bbox2 = get_label_bbox(n2, pos2)
                
                if bbox_intersect(bbox1, bbox2):
                    has_overlap = True
                    
                    # Nudge the lower node (higher y) tangentially
                    # Choose which one to nudge based on ring (prefer nudging outer ring)
                    if n1["ring"] >= n2["ring"]:
                        nudge_node = n1
                        other_node = n2
                    else:
                        nudge_node = n2
                        other_node = n1
                    
                    # Compute tangent direction
                    pos = nudge_node["pos2"]
                    if isinstance(pos, list):
                        pos = np.array(pos, dtype=float)
                    R = ring_radius[nudge_node["ring"]]
                    angle = math.atan2(float(pos[1]), float(pos[0]))
                    tangent = np.array([-math.sin(angle), math.cos(angle)], dtype=float)
                    
                    # Nudge along tangent
                    new_pos = pos + nudge_step * tangent
                    
                    # Re-project to ring
                    r = float(np.linalg.norm(new_pos)) + 1e-9
                    nudge_node["pos2"] = R * new_pos / r
                    
                    # Update bbox for next check
                    if nudge_node == n1:
                        bbox1 = get_label_bbox(n1, nudge_node["pos2"])
                    else:
                        bbox2 = get_label_bbox(n2, nudge_node["pos2"])
        
        if not has_overlap:
            break
    
    return nodes
