"""
Sector-Based Layout Engine
Places residues within angular sectors based on ligand atom anchors.
"""

import math
import numpy as np
from typing import List, Dict, Tuple

# Fixed ring radii
RING_RADIUS = {
    1: 170.0,  # Inner ring: H-bond, salt bridge
    2: 240.0,  # Middle ring: pi-stack, pi-cation
    3: 310.0,  # Outer ring: hydrophobic, others
}

# Node footprint parameters
NODE_DIAMETER = 40.0  # pixels (node radius * 2)
NODE_PADDING = 8.0  # pixels
MIN_NODE_SEPARATION = NODE_DIAMETER + NODE_PADDING  # d_min


def compute_minimum_angular_spacing(ring_radius: float) -> float:
    """
    Compute minimum angular spacing for nodes on a ring.
    
    Args:
        ring_radius: Radius of the ring in pixels
    
    Returns:
        Minimum angular spacing in radians
    """
    # delta = 2 * arcsin(d_min / (2 * R))
    d_min = MIN_NODE_SEPARATION
    R = ring_radius
    
    if R <= 0:
        return np.deg2rad(10.0)  # Fallback
    
    ratio = d_min / (2 * R)
    ratio = np.clip(ratio, -1.0, 1.0)  # Clamp for arcsin
    delta = 2 * np.arcsin(ratio)
    
    return max(delta, np.deg2rad(5.0))  # Minimum 5 degrees


def compute_sector_angles(
    residue_anchors: Dict[str, Dict],
    ring_radius_map: Dict[int, float] = None
) -> List[Dict]:
    """
    Compute adaptive angular sectors for each unique ligand atom anchor.
    
    Args:
        residue_anchors: Dict from map_interactions_to_ligand_atoms
        ring_radius_map: Dict mapping ring number to radius
    
    Returns:
        List of sector dicts with adaptive spans
    """
    if ring_radius_map is None:
        ring_radius_map = RING_RADIUS
    
    # Group residues by ligand atom anchor
    anchor_groups = {}
    
    for res_id, anchor_info in residue_anchors.items():
        lig_atom_idx = anchor_info["ligand_atom_index"]
        lig_atom_pos = anchor_info["ligand_atom_pos"]
        # Use tuple of list for key
        key = (lig_atom_idx, tuple(lig_atom_pos))
        
        if key not in anchor_groups:
            anchor_groups[key] = {
                "ligand_atom_index": lig_atom_idx,
                "ligand_atom_pos": lig_atom_pos,  # Already a list
                "center_angle": anchor_info["angle"],
                "residue_ids": [],
            }
        
        anchor_groups[key]["residue_ids"].append(res_id)
    
    # Compute required angular spans for each anchor
    anchor_spans = {}
    for key, group in anchor_groups.items():
        n_a = len(group["residue_ids"])
        
        # Use average ring radius for spacing calculation
        avg_ring_radius = np.mean(list(ring_radius_map.values()))
        delta_min = compute_minimum_angular_spacing(avg_ring_radius)
        
        # Required span: Theta_a = n_a * delta_min
        Theta_a = n_a * delta_min
        
        anchor_spans[key] = {
            "group": group,
            "required_span": Theta_a,
            "n_residues": n_a,
        }
    
    # Normalize spans so total <= 0.9 * 2π
    total_span = sum(a["required_span"] for a in anchor_spans.values())
    max_total_span = 0.9 * 2 * np.pi
    
    if total_span > max_total_span:
        # Scale down proportionally
        scale_factor = max_total_span / total_span
        for key in anchor_spans:
            anchor_spans[key]["required_span"] *= scale_factor
    
    # Sort anchors by center angle
    sorted_anchors = sorted(
        anchor_spans.items(),
        key=lambda x: x[1]["group"]["center_angle"]
    )
    
    # Allocate sectors sequentially around circle
    sectors = []
    current_angle = 0.0
    
    for key, span_info in sorted_anchors:
        group = span_info["group"]
        Theta_a = span_info["required_span"]
        
        sector_start = current_angle
        sector_end = current_angle + Theta_a
        
        # Convert numpy array to list for JSON serialization
        lig_atom_pos_list = group["ligand_atom_pos"]
        if isinstance(lig_atom_pos_list, np.ndarray):
            lig_atom_pos_list = [float(lig_atom_pos_list[0]), float(lig_atom_pos_list[1])]
        
        sectors.append({
            "center_angle": float(group["center_angle"]),
            "start_angle": float(sector_start),
            "end_angle": float(sector_end),
            "required_span": float(Theta_a),
            "residue_ids": group["residue_ids"],
            "ligand_atom_pos": lig_atom_pos_list,
            "ligand_atom_index": int(group["ligand_atom_index"]),
        })
        
        current_angle = sector_end
    
    return sectors


def assign_ring_by_interaction(interactions: List[Dict]) -> int:
    """
    Assign residue to ring based on strongest interaction type.
    
    Returns:
        Ring number (1, 2, or 3)
    """
    ring_map = {
        "hbond": 1,
        "salt_bridge": 1,
        "metal_coordination": 1,
        "halogen_bond": 1,
        "pi_pi": 2,
        "pi_cation": 2,
        "hydrophobic": 3,
        "distance": 3,
    }
    
    rings = [ring_map.get(it["type"], 3) for it in interactions]
    return min(rings)  # Use strongest (lowest ring number)


def place_residues_in_sectors(
    sectors: List[Dict],
    residue_anchors: Dict[str, Dict],
    interaction_graph: Dict
) -> List[Dict]:
    """
    Place residues within their assigned sectors on appropriate rings.
    
    Args:
        sectors: List of sector dicts from compute_sector_angles
        residue_anchors: Dict from map_interactions_to_ligand_atoms
        interaction_graph: Graph with "nodes" and "edges"
    
    Returns:
        List of residue nodes with "pos2" (2D position) and "sector" info
    """
    # Create residue node map
    residue_nodes = {}
    for node in interaction_graph.get("nodes", []):
        if node.get("type") == "residue":
            res_id = node.get("id")
            if res_id:
                residue_nodes[res_id] = node
    
    # Validate inputs
    if not sectors:
        # If no sectors, try to create a simple layout from residue nodes
        if not residue_nodes:
            return []
        # Fallback: place all residues evenly around circle
        placed_residues = []
        residue_list = list(residue_nodes.values())
        num_residues = len(residue_list)
        if num_residues == 0:
            return []
        
        for idx, node in enumerate(residue_list):
            angle = 2 * np.pi * idx / num_residues
            ring = node.get("ring", 3)
            R = RING_RADIUS.get(ring, RING_RADIUS[3])
            pos2 = np.array([
                R * np.cos(angle),
                R * np.sin(angle)
            ], dtype=float)
            node["pos2"] = pos2
            node["sector_angle"] = float(angle)
            placed_residues.append(node)
        return placed_residues
    
    # Place residues within sectors
    placed_residues = []
    
    # Debug: Check if residue IDs match
    anchor_res_ids = set(residue_anchors.keys())
    node_res_ids = set(residue_nodes.keys())
    missing_in_nodes = anchor_res_ids - node_res_ids
    missing_in_anchors = node_res_ids - anchor_res_ids
    
    # If there are mismatches, log them but try to continue
    if missing_in_nodes or missing_in_anchors:
        import logging
        logger = logging.getLogger(__name__)
        if missing_in_nodes:
            logger.warning(f"Residue IDs in anchors but not in nodes: {missing_in_nodes}")
        if missing_in_anchors:
            logger.warning(f"Residue IDs in nodes but not in anchors: {missing_in_anchors}")
        
        # Only use residues that exist in both
        common_res_ids = anchor_res_ids & node_res_ids
        if not common_res_ids:
            # No common residues - this is a critical error
            raise ValueError(
                f"No matching residue IDs between graph and anchors. "
                f"Graph has: {sorted(node_res_ids)[:10]}, "
                f"Anchors have: {sorted(anchor_res_ids)[:10]}"
            )
        
        # Filter sectors to only include common residues
        filtered_sectors = []
        for sector in sectors:
            sector_residue_ids = sector.get("residue_ids", [])
            filtered_ids = [rid for rid in sector_residue_ids if rid in common_res_ids]
            if filtered_ids:
                sector_copy = sector.copy()
                sector_copy["residue_ids"] = filtered_ids
                filtered_sectors.append(sector_copy)
        sectors = filtered_sectors
    
    for sector in sectors:
        sector_residue_ids = sector.get("residue_ids", [])
        
        if not sector_residue_ids:
            continue
        
        # Get residue nodes for this sector
        sector_residues = []
        for res_id in sector_residue_ids:
            if res_id in residue_nodes:
                node = residue_nodes[res_id]
                # Get interactions for this residue
                if res_id in residue_anchors:
                    interactions = residue_anchors[res_id].get("interactions", [])
                    if interactions:
                        ring = assign_ring_by_interaction(interactions)
                    else:
                        ring = node.get("ring", 3)  # Use ring from node if no interactions
                else:
                    ring = node.get("ring", 3)  # Default to outer ring
                
                sector_residues.append({
                    "node": node,
                    "res_id": res_id,
                    "ring": ring,
                })
            # If res_id not in residue_nodes, skip it (might be a mismatch)
        
        if not sector_residues:
            continue
        
        # Group by ring within sector
        residues_by_ring = {1: [], 2: [], 3: []}
        for res_info in sector_residues:
            ring = res_info["ring"]
            residues_by_ring[ring].append(res_info)
        
        # Place residues within sector for each ring
        for ring_num in [1, 2, 3]:
            ring_residues = residues_by_ring[ring_num]
            if not ring_residues:
                continue
            
            # Compute angular range for this ring within sector
            sector_start = sector["start_angle"]
            sector_end = sector["end_angle"]
            
            # Handle wrap-around
            if sector_end < sector_start:
                sector_end += 2 * np.pi
            
            sector_span = sector_end - sector_start
            num_residues = len(ring_residues)
            
            # Compute minimum angular spacing for this ring
            R = RING_RADIUS[ring_num]
            delta_min = compute_minimum_angular_spacing(R)
            
            # Evenly distribute within sector
            if num_residues == 1:
                assigned_angle = sector_start + sector_span / 2
            else:
                # Use minimum spacing
                min_span = (num_residues - 1) * delta_min
                if sector_span < min_span:
                    # Expand sector to accommodate minimum spacing
                    expansion = (min_span - sector_span) / 2
                    sector_start -= expansion
                    sector_end += expansion
                    sector_span = sector_end - sector_start
                
                # Evenly distribute with padding at edges
                angle_step = sector_span / (num_residues + 1)
                
                for idx, res_info in enumerate(ring_residues):
                    assigned_angle = sector_start + (idx + 1) * angle_step
                    
                    # Normalize angle to [0, 2π]
                    assigned_angle = ((assigned_angle % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
                    
                    # Compute position on ring
                    pos2 = np.array([
                        R * np.cos(assigned_angle),
                        R * np.sin(assigned_angle)
                    ], dtype=float)
                    
                    node = res_info["node"]
                    node["pos2"] = pos2
                    node["ring"] = ring_num
                    node["sector_angle"] = float(assigned_angle)
                    node["sector_center_angle"] = float(sector["center_angle"])
                    
                    placed_residues.append(node)
    
    return placed_residues


def hard_ring_projection(nodes: List[Dict], ring_radius_map: Dict[int, float]) -> List[Dict]:
    """
    Project all nodes back onto their assigned ring radius.
    
    Args:
        nodes: List of nodes with "pos2" and "ring"
        ring_radius_map: Dict mapping ring number to radius
    
    Returns:
        Updated nodes with positions projected onto rings
    """
    for node in nodes:
        if node.get("type") == "ligand":
            continue
        
        ring = node.get("ring", 3)
        R = ring_radius_map.get(ring, RING_RADIUS[3])
        
        pos2 = node.get("pos2")
        if pos2 is None:
            continue
        
        # Convert to numpy array if needed
        if isinstance(pos2, list):
            pos2 = np.array(pos2, dtype=float)
        
        # Project onto ring: p_i = R * p_i / ||p_i||
        norm = np.linalg.norm(pos2) + 1e-9
        pos2_projected = R * pos2 / norm
        
        node["pos2"] = pos2_projected
    
    return nodes


def count_edge_crossings(
    placed_residues: List[Dict],
    residue_anchors: Dict[str, Dict],
    ligand_center: Tuple[float, float]
) -> int:
    """
    Count number of edge crossings between interaction edges.
    
    Args:
        placed_residues: List of placed residue nodes
        residue_anchors: Dict mapping residue_id -> anchor info
        ligand_center: (x, y) coordinates of ligand center
    
    Returns:
        Number of edge crossings
    """
    crossings = 0
    lig_cx, lig_cy = ligand_center
    
    # Get residue positions
    residue_positions = {}
    for res in placed_residues:
        res_id = res.get("id")
        pos2 = res.get("pos2")
        if isinstance(pos2, list):
            pos2 = np.array(pos2, dtype=float)
        residue_positions[res_id] = pos2
    
    # Check all pairs of edges
    residue_ids = list(residue_positions.keys())
    for i in range(len(residue_ids)):
        for j in range(i + 1, len(residue_ids)):
            res_id1 = residue_ids[i]
            res_id2 = residue_ids[j]
            
            if res_id1 not in residue_anchors or res_id2 not in residue_anchors:
                continue
            
            pos1 = residue_positions[res_id1]
            pos2 = residue_positions[res_id2]
            
            anchor1 = residue_anchors[res_id1]["ligand_atom_pos"]
            anchor2 = residue_anchors[res_id2]["ligand_atom_pos"]
            
            if isinstance(anchor1, list):
                anchor1 = np.array(anchor1, dtype=float)
            if isinstance(anchor2, list):
                anchor2 = np.array(anchor2, dtype=float)
            
            # Check if edges cross (simplified: check if line segments intersect)
            # Edge 1: pos1 -> anchor1
            # Edge 2: pos2 -> anchor2
            # Using line segment intersection test
            def segments_intersect(p1, p2, p3, p4):
                """Check if line segments p1-p2 and p3-p4 intersect."""
                def ccw(A, B, C):
                    return (C[1] - A[1]) * (B[0] - A[0]) > (B[1] - A[1]) * (C[0] - A[0])
                return ccw(p1, p3, p4) != ccw(p2, p3, p4) and ccw(p1, p2, p3) != ccw(p1, p2, p4)
            
            if segments_intersect(pos1, anchor1, pos2, anchor2):
                crossings += 1
    
    return crossings


def optimize_sector_positions(
    placed_residues: List[Dict],
    sectors: List[Dict],
    residue_anchors: Dict[str, Dict] = None,
    ligand_center: Tuple[float, float] = None
) -> List[Dict]:
    """
    Global angular optimization: rotate sector blocks to minimize edge crossings and gaps.
    
    Preserves intra-sector ordering while allowing sector blocks to rotate.
    
    Args:
        placed_residues: List of placed residue nodes
        sectors: List of sector dicts
        residue_anchors: Optional dict for edge crossing calculation
        ligand_center: Optional ligand center for edge crossing calculation
    
    Returns:
        Updated placed_residues with optimized positions
    """
    if len(placed_residues) < 2:
        return placed_residues
    
    # Group residues by sector
    residues_by_sector = {}
    for res in placed_residues:
        sector_id = res.get("sector_center_angle", 0.0)
        if sector_id not in residues_by_sector:
            residues_by_sector[sector_id] = []
        residues_by_sector[sector_id].append(res)
    
    # Sort all residues by current angle
    all_residues = sorted(placed_residues, key=lambda r: r.get("sector_angle", 0.0))
    
    # Compute gaps between consecutive residues
    gaps = []
    for i in range(len(all_residues) - 1):
        angle1 = all_residues[i].get("sector_angle", 0.0)
        angle2 = all_residues[i + 1].get("sector_angle", 0.0)
        
        # Normalize angles
        angle1 = ((angle1 % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
        angle2 = ((angle2 % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
        
        if angle2 < angle1:
            angle2 += 2 * np.pi
        
        gap = angle2 - angle1
        gaps.append((i, gap))
    
    # Find maximum gap
    max_gap = max([g[1] for g in gaps], default=0.0)
    
    # Count edge crossings if anchors provided
    crossings = 0
    if residue_anchors and ligand_center:
        crossings = count_edge_crossings(placed_residues, residue_anchors, ligand_center)
    
    # Optimize: try rotating sectors to reduce gaps and crossings
    best_residues = [r.copy() for r in placed_residues]
    best_score = max_gap + crossings * 0.1  # Weight crossings less than gaps
    
    # Try small rotations of each sector block
    for sector_id, sector_residues in residues_by_sector.items():
        if len(sector_residues) == 0:
            continue
        
        # Try rotating this sector by small increments
        for rotation in np.linspace(-np.pi / 6, np.pi / 6, 7):  # ±30 degrees
            # Create test configuration
            test_residues = []
            for res in placed_residues:
                if res in sector_residues:
                    # Rotate this residue
                    new_res = res.copy()
                    old_angle = new_res.get("sector_angle", 0.0)
                    new_angle = old_angle + rotation
                    new_angle = ((new_angle % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
                    
                    ring = new_res.get("ring", 3)
                    R = RING_RADIUS[ring]
                    new_res["pos2"] = np.array([
                        R * np.cos(new_angle),
                        R * np.sin(new_angle)
                    ], dtype=float)
                    new_res["sector_angle"] = float(new_angle)
                    test_residues.append(new_res)
                else:
                    test_residues.append(res.copy())
            
            # Compute score
            test_gaps = []
            sorted_test = sorted(test_residues, key=lambda r: r.get("sector_angle", 0.0))
            for i in range(len(sorted_test) - 1):
                angle1 = sorted_test[i].get("sector_angle", 0.0)
                angle2 = sorted_test[i + 1].get("sector_angle", 0.0)
                angle1 = ((angle1 % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
                angle2 = ((angle2 % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
                if angle2 < angle1:
                    angle2 += 2 * np.pi
                test_gaps.append(angle2 - angle1)
            
            test_max_gap = max(test_gaps, default=0.0)
            test_crossings = 0
            if residue_anchors and ligand_center:
                test_crossings = count_edge_crossings(test_residues, residue_anchors, ligand_center)
            
            test_score = test_max_gap + test_crossings * 0.1
            
            if test_score < best_score:
                best_score = test_score
                best_residues = test_residues
    
    return best_residues


def compute_sector_layout(
    interaction_graph: Dict,
    residue_anchors: Dict[str, Dict],
    ligand_center: Tuple[float, float],
    sector_width: float = None  # Deprecated, now computed adaptively
) -> List[Dict]:
    """
    Main function: compute sector-based layout with adaptive spans and optimization.
    
    Args:
        interaction_graph: Graph with "nodes" and "edges"
        residue_anchors: Dict from map_interactions_to_ligand_atoms
        ligand_center: (x, y) coordinates of ligand center
        sector_width: Deprecated (kept for compatibility)
    
    Returns:
        List of placed residue nodes with "pos2" positions
    """
    # Validate inputs
    if not interaction_graph or "nodes" not in interaction_graph:
        return []
    
    if not residue_anchors:
        return []
    
    from .interaction_anchor import compute_sector_angles
    
    # Compute adaptive sectors
    sectors = compute_sector_angles(residue_anchors, RING_RADIUS)
    
    # Place residues within sectors
    placed_residues = place_residues_in_sectors(sectors, residue_anchors, interaction_graph)
    
    # If no residues placed, return empty list (will be caught by caller)
    if not placed_residues:
        return []
    
    # Global angular optimization: minimize edge crossings and gaps
    placed_residues = optimize_sector_positions(
        placed_residues, sectors, residue_anchors=residue_anchors, ligand_center=ligand_center
    )
    
    # Hard ring projection
    placed_residues = hard_ring_projection(placed_residues, RING_RADIUS)
    
    return placed_residues
