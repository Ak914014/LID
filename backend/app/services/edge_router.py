"""
Edge Router
Routes interaction edges using Bezier curves to specific ligand atom anchors.
"""

import numpy as np
from typing import List, Dict, Tuple

EPS = 1e-9
CENTER_AVOIDANCE_THRESHOLD = 40.0  # pixels
CENTER_AVOIDANCE_CURVATURE = 0.5  # Increased curvature when near center
NORMALIZED_EDGE_LENGTH = 200.0  # Fixed visual length for normalization
LIGAND_EXCLUSION_PADDING = 15.0  # Additional padding around ligand


def normalize_edge_length(
    residue_pos: np.ndarray,
    ligand_anchor_pos: np.ndarray,
    ligand_center: Tuple[float, float],
    target_length: float = NORMALIZED_EDGE_LENGTH
) -> np.ndarray:
    """
    Normalize edge visual length by scaling the ligand anchor position.
    
    Args:
        residue_pos: 2D position of residue [x, y]
        ligand_anchor_pos: 2D position of ligand atom anchor [x, y]
        ligand_center: (x, y) coordinates of ligand center
        target_length: Target visual length in pixels
    
    Returns:
        Normalized ligand anchor position
    """
    lig_cx, lig_cy = ligand_center
    
    # Vector from residue to ligand anchor
    edge_vec = ligand_anchor_pos - residue_pos
    edge_length = np.linalg.norm(edge_vec) + EPS
    
    # If edge is too short, use original
    if edge_length < target_length * 0.5:
        return ligand_anchor_pos
    
    # Scale to target length
    scale_factor = target_length / edge_length
    
    # Compute normalized anchor position
    normalized_anchor = residue_pos + edge_vec * scale_factor
    
    return normalized_anchor


def compute_cubic_bezier_control_points(
    residue_pos: np.ndarray,
    ligand_anchor_pos: np.ndarray,
    ligand_center: Tuple[float, float],
    ligand_exclusion_radius: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute cubic Bezier control points that route outside ligand exclusion zone.
    
    Args:
        residue_pos: 2D position of residue [x, y]
        ligand_anchor_pos: 2D position of ligand atom anchor [x, y]
        ligand_center: (x, y) coordinates of ligand center
        ligand_exclusion_radius: Radius of ligand exclusion region
    
    Returns:
        (control_point_1, control_point_2) for cubic Bezier
    """
    lig_cx, lig_cy = ligand_center
    center = np.array([lig_cx, lig_cy], dtype=float)
    
    # Vector from residue to ligand anchor
    edge_vec = ligand_anchor_pos - residue_pos
    edge_length = np.linalg.norm(edge_vec) + EPS
    
    # Perpendicular direction
    perp = np.array([-edge_vec[1], edge_vec[0]], dtype=float)
    perp = perp / (np.linalg.norm(perp) + EPS)
    
    # Midpoint
    M = (residue_pos + ligand_anchor_pos) / 2
    
    # Check if midpoint is inside exclusion zone
    dist_mid_to_center = np.linalg.norm(M - center)
    
    if dist_mid_to_center < ligand_exclusion_radius:
        # Route via outer control points
        # Direction from center to midpoint
        dir_to_mid = M - center
        dir_to_mid = dir_to_mid / (np.linalg.norm(dir_to_mid) + EPS)
        
        # Place control points outside exclusion region
        outer_radius = ligand_exclusion_radius + 30.0
        cp1_pos = center + outer_radius * dir_to_mid + 20.0 * perp
        cp2_pos = center + outer_radius * dir_to_mid - 20.0 * perp
        
        # Adjust to be at 1/3 and 2/3 along the path
        t1, t2 = 1/3, 2/3
        cp1 = residue_pos + t1 * (ligand_anchor_pos - residue_pos) + (cp1_pos - M) * 0.5
        cp2 = residue_pos + t2 * (ligand_anchor_pos - residue_pos) + (cp2_pos - M) * 0.5
    else:
        # Standard cubic Bezier with curvature
        lambda_factor = min(edge_length * 0.3, 50.0)
        cp1 = residue_pos + edge_vec * 0.33 + lambda_factor * perp
        cp2 = residue_pos + edge_vec * 0.67 - lambda_factor * perp
    
    return cp1, cp2


def compute_bezier_control_point(
    residue_pos: np.ndarray,
    ligand_anchor_pos: np.ndarray,
    ligand_center: Tuple[float, float],
    interaction_type: str,
    ligand_exclusion_radius: float = None
) -> Tuple[np.ndarray, float]:
    """
    Compute Bezier control point for curved edge with clamped curvature.
    
    Args:
        residue_pos: 2D position of residue [x, y]
        ligand_anchor_pos: 2D position of ligand atom anchor [x, y]
        ligand_center: (x, y) coordinates of ligand center
        interaction_type: Type of interaction
    
    Returns:
        (control_point, curvature_factor)
    """
    # Midpoint
    M = (residue_pos + ligand_anchor_pos) / 2
    
    # Perpendicular direction (normalized)
    edge_vec = ligand_anchor_pos - residue_pos
    edge_length = np.linalg.norm(edge_vec) + EPS
    
    # Perpendicular vector: rotate 90 degrees
    perp = np.array([-edge_vec[1], edge_vec[0]], dtype=float)
    perp = perp / (np.linalg.norm(perp) + EPS)
    
    # Base curvature factor
    curvature_map = {
        "hbond": 0.4,
        "pi_pi": 0.35,
        "pi_cation": 0.35,
        "hydrophobic": 0.2,
        "salt_bridge": 0.0,  # Straight line
        "metal_coordination": 0.3,
        "halogen_bond": 0.3,
        "distance": 0.25,
    }
    
    base_curvature = curvature_map.get(interaction_type, 0.25)
    
    # Clamp curvature: lambda = clamp(edge_length * base_curvature, 30, 60)
    lambda_factor = edge_length * base_curvature
    lambda_factor = np.clip(lambda_factor, 30.0, 60.0)
    
    # Check if edge passes through ligand exclusion region
    lig_cx, lig_cy = ligand_center
    center_pos = np.array([lig_cx, lig_cy])
    
    if ligand_exclusion_radius is not None:
        # Check if midpoint is inside exclusion region
        dist_mid_to_center = np.linalg.norm(M - center_pos)
        
        if dist_mid_to_center < ligand_exclusion_radius:
            # Route via outer control point outside ligand boundary
            # Compute direction from center to midpoint
            dir_to_mid = M - center_pos
            dir_to_mid = dir_to_mid / (np.linalg.norm(dir_to_mid) + EPS)
            
            # Place control point outside exclusion region
            outer_point = center_pos + (ligand_exclusion_radius + 25.0) * dir_to_mid
            
            # Use outer point as control, but adjust lambda to reach it
            lambda_factor = np.linalg.norm(outer_point - M)
            lambda_factor = np.clip(lambda_factor, 45.0, 90.0)
            
            # Adjust perpendicular direction to point toward outer point
            perp = outer_point - M
            perp = perp / (np.linalg.norm(perp) + EPS)
        else:
            # Check if line segment intersects exclusion region
            # Distance from center to line segment
            a = edge_vec[1]
            b = -edge_vec[0]
            c = residue_pos[0] * edge_vec[1] - residue_pos[1] * edge_vec[0]
            dist_to_center = abs(
                a * lig_cx + b * lig_cy + c
            ) / (np.sqrt(a * a + b * b) + EPS)
            
            if dist_to_center < ligand_exclusion_radius + CENTER_AVOIDANCE_THRESHOLD:
                # Increase curvature to avoid ligand (but still clamp)
                lambda_factor = np.clip(
                    max(lambda_factor, edge_length * CENTER_AVOIDANCE_CURVATURE),
                    45.0, 90.0
                )
    else:
        # Fallback: check distance to center
        a = edge_vec[1]
        b = -edge_vec[0]
        c = residue_pos[0] * edge_vec[1] - residue_pos[1] * edge_vec[0]
        dist_to_center = abs(
            a * lig_cx + b * lig_cy + c
        ) / (np.sqrt(a * a + b * b) + EPS)
        
        if dist_to_center < CENTER_AVOIDANCE_THRESHOLD:
            lambda_factor = np.clip(
                max(lambda_factor, edge_length * CENTER_AVOIDANCE_CURVATURE),
                30.0, 60.0
            )
    
    # Control point
    control_point = M + lambda_factor * perp
    
    return control_point, base_curvature


def create_interaction_edge_path(
    residue_pos: np.ndarray,
    ligand_anchor_pos: np.ndarray,
    ligand_center: Tuple[float, float],
    interaction_type: str,
    normalize_length: bool = True,
    ligand_exclusion_radius: float = None
) -> str:
    """
    Create SVG path for interaction edge, avoiding ligand interior.
    Uses cubic Bezier if quadratic Bezier still intersects exclusion zone.
    
    Args:
        residue_pos: 2D position of residue [x, y]
        ligand_anchor_pos: 2D position of ligand atom anchor [x, y]
        ligand_center: (x, y) coordinates of ligand center
        interaction_type: Type of interaction
        normalize_length: If True, normalize edge visual length
        ligand_exclusion_radius: Radius of ligand exclusion region
    
    Returns:
        SVG path string
    """
    # Normalize edge length if requested
    if normalize_length:
        ligand_anchor_pos = normalize_edge_length(
            residue_pos, ligand_anchor_pos, ligand_center
        )
    
    x1, y1 = float(residue_pos[0]), float(residue_pos[1])
    x2, y2 = float(ligand_anchor_pos[0]), float(ligand_anchor_pos[1])
    
    # ENFORCE CURVED ROUTING: Check if straight segment intersects ligand exclusion
    if ligand_exclusion_radius is not None:
        lig_cx, lig_cy = ligand_center
        center = np.array([lig_cx, lig_cy], dtype=float)
        
        # Check if straight line segment intersects exclusion circle
        # Using point-to-line distance formula
        edge_vec = ligand_anchor_pos - residue_pos
        a = edge_vec[1]
        b = -edge_vec[0]
        c = residue_pos[0] * edge_vec[1] - residue_pos[1] * edge_vec[0]
        dist_to_center = abs(
            a * lig_cx + b * lig_cy + c
        ) / (np.sqrt(a * a + b * b) + EPS)
        
        # Check if midpoint is inside exclusion zone
        midpoint = (residue_pos + ligand_anchor_pos) / 2
        dist_mid_to_center = np.linalg.norm(midpoint - center)
        
        # If straight segment would intersect, force Bezier curvature
        if dist_to_center < ligand_exclusion_radius or dist_mid_to_center < ligand_exclusion_radius:
            # Force Bezier routing (even for salt bridge)
            pass  # Continue to Bezier routing below
        elif interaction_type == "salt_bridge":
            # Only use straight line if it doesn't intersect exclusion zone
            return f"M {x1} {y1} L {x2} {y2}"
    elif interaction_type == "salt_bridge":
        # No exclusion zone defined, use straight line
        return f"M {x1} {y1} L {x2} {y2}"
    
    # Try quadratic Bezier first
    control_point, _ = compute_bezier_control_point(
        residue_pos, ligand_anchor_pos, ligand_center, interaction_type,
        ligand_exclusion_radius=ligand_exclusion_radius
    )
    
    # Check if quadratic Bezier still intersects exclusion zone
    if ligand_exclusion_radius is not None:
        lig_cx, lig_cy = ligand_center
        center = np.array([lig_cx, lig_cy], dtype=float)
        
        # Sample points along quadratic Bezier curve
        intersects = False
        for t in np.linspace(0, 1, 10):
            # Quadratic Bezier: (1-t)^2*P0 + 2*(1-t)*t*P1 + t^2*P2
            p0 = residue_pos
            p1 = control_point
            p2 = ligand_anchor_pos
            point = (1-t)**2 * p0 + 2*(1-t)*t * p1 + t**2 * p2
            dist = np.linalg.norm(point - center)
            if dist < ligand_exclusion_radius:
                intersects = True
                break
        
        if intersects:
            # Use cubic Bezier with control points outside exclusion zone
            cp1, cp2 = compute_cubic_bezier_control_points(
                residue_pos, ligand_anchor_pos, ligand_center, ligand_exclusion_radius
            )
            cx1, cy1 = float(cp1[0]), float(cp1[1])
            cx2, cy2 = float(cp2[0]), float(cp2[1])
            return f"M {x1} {y1} C {cx1} {cy1}, {cx2} {cy2}, {x2} {y2}"
    
    # Use quadratic Bezier
    cx, cy = float(control_point[0]), float(control_point[1])
    return f"M {x1} {y1} Q {cx} {cy}, {x2} {y2}"


def bundle_hydrophobic_edges(
    edges: List[Dict],
    residue_nodes: Dict[str, Dict],
    ligand_center: Tuple[float, float],
    ligand_exclusion_radius: float = None
) -> List[Dict]:
    """
    Bundle multiple hydrophobic edges from same residue to same anchor.
    
    Args:
        edges: List of edge dicts
        residue_nodes: Dict mapping residue_id -> node
        ligand_center: (x, y) coordinates of ligand center
        ligand_exclusion_radius: Radius of ligand exclusion region
    
    Returns:
        List of edges with bundled hydrophobic interactions
    """
    # Group hydrophobic edges by (residue_id, ligand_atom_index)
    hydrophobic_groups = {}
    other_edges = []
    
    for edge in edges:
        if edge["type"] == "hydrophobic":
            key = (edge["residue_id"], edge["ligand_atom_index"])
            if key not in hydrophobic_groups:
                hydrophobic_groups[key] = []
            hydrophobic_groups[key].append(edge)
        else:
            other_edges.append(edge)
    
    # Bundle each group
    bundled_edges = []
    for key, group_edges in hydrophobic_groups.items():
        if len(group_edges) == 1:
            # Single edge, no bundling needed
            bundled_edges.append(group_edges[0])
        else:
            # Bundle multiple edges into one
            first_edge = group_edges[0]
            res_id = first_edge["residue_id"]
            
            if res_id not in residue_nodes:
                # Fallback: keep all edges
                bundled_edges.extend(group_edges)
                continue
            
            residue_node = residue_nodes[res_id]
            residue_pos = residue_node["pos2"]
            
            # Convert to numpy array if needed
            if isinstance(residue_pos, list):
                residue_pos = np.array(residue_pos, dtype=float)
            
            # Use first edge's anchor position
            ligand_anchor_pos = first_edge["ligand_atom_pos"]
            if isinstance(ligand_anchor_pos, list):
                ligand_anchor_pos = np.array(ligand_anchor_pos, dtype=float)
            
            # Create single bundled path (dashed arc)
            path = create_interaction_edge_path(
                residue_pos,
                ligand_anchor_pos,
                ligand_center,
                "hydrophobic",
                normalize_length=True,
                ligand_exclusion_radius=ligand_exclusion_radius
            )
            
            # Compute average distance
            avg_distance = np.mean([e["distance"] for e in group_edges])
            count = len(group_edges)
            
            bundled_edges.append({
                "residue_id": res_id,
                "path": path,
                "type": "hydrophobic",
                "backbone": False,
                "distance": float(avg_distance),
                "ligand_atom_index": first_edge["ligand_atom_index"],
                "ligand_atom_pos": first_edge["ligand_atom_pos"],
                "bundled_count": count,  # Number of bundled interactions
            })
    
    return other_edges + bundled_edges


def route_all_edges(
    interactions: List[Dict],
    residue_nodes: Dict[str, Dict],
    residue_anchors: Dict[str, Dict],
    ligand_center: Tuple[float, float],
    ligand_atom_xy: List[Dict] = None
) -> List[Dict]:
    """
    Route all interaction edges with hydrophobic bundling and ligand avoidance.
    
    Args:
        interactions: List of interaction dicts
        residue_nodes: Dict mapping residue_id -> node with "pos2"
        residue_anchors: Dict from map_interactions_to_ligand_atoms
        ligand_center: (x, y) coordinates of ligand center
        ligand_atom_xy: List of ligand atom positions for exclusion region
    
    Returns:
        List of edge dicts with "path", "type", "backbone", etc.
    """
    # Compute ligand exclusion region
    ligand_exclusion_radius = None
    if ligand_atom_xy:
        ligand_exclusion_radius, _ = compute_ligand_exclusion_region(
            ligand_atom_xy, ligand_center
        )
    
    # Include all interaction types (don't filter - show all meaningful interactions)
    edges = []
    
    for it in interactions:
        res_id = f"{it['resname']}{it['resid']}"
        
        if res_id not in residue_nodes:
            continue
        
        residue_node = residue_nodes[res_id]
        residue_pos = residue_node["pos2"]
        
        if res_id not in residue_anchors:
            continue
        
        anchor_info = residue_anchors[res_id]
        ligand_anchor_pos = anchor_info["ligand_atom_pos"]
        
        # Convert to numpy array if it's a list (for calculations)
        if isinstance(ligand_anchor_pos, list):
            ligand_anchor_pos_np = np.array(ligand_anchor_pos, dtype=float)
        else:
            ligand_anchor_pos_np = ligand_anchor_pos
        
        # Convert residue_pos to numpy array if needed
        if isinstance(residue_pos, list):
            residue_pos_np = np.array(residue_pos, dtype=float)
        else:
            residue_pos_np = residue_pos
        
        # Create edge path (with ligand exclusion)
        path = create_interaction_edge_path(
            residue_pos_np,
            ligand_anchor_pos_np,
            ligand_center,
            it["type"],
            normalize_length=True,
            ligand_exclusion_radius=ligand_exclusion_radius
        )
        
        # Ensure ligand_anchor_pos is a list for JSON serialization
        if isinstance(ligand_anchor_pos, np.ndarray):
            ligand_anchor_pos = [float(ligand_anchor_pos[0]), float(ligand_anchor_pos[1])]
        
        edges.append({
            "residue_id": res_id,
            "path": path,
            "type": it["type"],
            "backbone": it.get("backbone", False),
            "distance": float(it.get("distance", 0.0)),
            "ligand_atom_index": int(anchor_info["ligand_atom_index"]),
            "ligand_atom_pos": [float(ligand_anchor_pos[0]), float(ligand_anchor_pos[1])],  # Convert numpy array to list
        })
    
    # Bundle hydrophobic edges
    edges = bundle_hydrophobic_edges(edges, residue_nodes, ligand_center, ligand_exclusion_radius)
    
    return edges
