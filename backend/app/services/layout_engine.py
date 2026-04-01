"""
Layout Engine
Computes radial layout for residues around ligand center.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional


def compute_radial_layout(
    residues: List[Dict],
    interactions: List[Dict],
    ligand_center: Tuple[float, float],
    ligand_atoms: Optional[List[Dict]] = None,
    min_radius: float = 200.0,
    max_radius: float = 400.0,
    min_angular_spacing: float = 0.12,
    num_sectors: int = 8,
    inter_cluster_gap: float = 0.15
) -> List[Dict]:
    """
    Compute sector-based radial layout for residues.
    
    Args:
        residues: List of residue dicts with angle, dist, etc.
        interactions: List of interaction dicts
        ligand_center: (x, y) coordinates of ligand center
        ligand_atoms: List of ligand atom positions
        min_radius: Minimum radial distance
        max_radius: Maximum radial distance
        min_angular_spacing: Minimum angular spacing in radians
        num_sectors: Number of angular sectors for clustering
        inter_cluster_gap: Gap between clusters in radians
    
    Returns:
        List of residue dicts with x, y coordinates
    """
    if not residues:
        return []
    
    lig_center_x, lig_center_y = ligand_center
    
    # Compute interaction anchors for clustering
    residues_with_anchors = []
    for res in residues:
        # Find interaction anchor point
        anchor_angle = res.get("angle", 0.0)
        anchor_x = None
        anchor_y = None
        
        # Try to find interaction anchor from interactions
        res_interactions = [
            it for it in interactions 
            if it.get("residue") == f"{res.get('resname', '')}{res.get('resid', '')}"
        ]
        
        if res_interactions and ligand_atoms:
            first_it = res_interactions[0]
            lig_atom_idx = first_it.get("ligand_atom_index", -1)
            if 0 <= lig_atom_idx < len(ligand_atoms):
                anchor_x = ligand_atoms[lig_atom_idx].get("x")
                anchor_y = ligand_atoms[lig_atom_idx].get("y")
                if anchor_x is not None and anchor_y is not None:
                    anchor_angle = np.arctan2(
                        anchor_y - lig_center_y,
                        anchor_x - lig_center_x
                    )
        
        residues_with_anchors.append({
            **res,
            "anchor_angle": anchor_angle,
            "anchor_x": anchor_x,
            "anchor_y": anchor_y,
        })
    
    # Cluster by sector
    clusters = _cluster_by_sector(residues_with_anchors, num_sectors)
    
    # Allocate arc ranges
    allocated_clusters = _allocate_arc_ranges(clusters, inter_cluster_gap)
    
    # Place residues within arcs
    layout_residues = []
    for cluster in allocated_clusters:
        cluster_residues = cluster["residues"]
        arc_start = cluster["arcStart"]
        arc_size = cluster["arcSize"]
        
        if not cluster_residues:
            continue
        
        spacing = max(arc_size / len(cluster_residues), min_angular_spacing)
        
        for idx, res in enumerate(cluster_residues):
            dist_3d = res.get("dist", 0.0)
            normalized_dist = np.clip(dist_3d, 0, 20)
            radius = min_radius + (normalized_dist / 20.0) * (max_radius - min_radius)
            
            assigned_angle = arc_start + idx * spacing
            
            x = lig_center_x + radius * np.cos(assigned_angle)
            y = lig_center_y + radius * np.sin(assigned_angle)
            
            layout_residues.append({
                **res,
                "angle": float(assigned_angle),
                "radius": float(radius),
                "x": float(x),
                "y": float(y),
            })
    
    # Resolve collisions
    layout_residues = _resolve_collisions(
        layout_residues, 
        ligand_center, 
        min_radius, 
        max_radius
    )
    
    return layout_residues


def _cluster_by_sector(residues: List[Dict], num_sectors: int) -> List[Dict]:
    """Cluster residues by angular sector."""
    # Normalize angles
    for res in residues:
        angle = res.get("anchor_angle", 0.0)
        res["anchor_angle"] = ((angle % (2 * np.pi)) + 2 * np.pi) % (2 * np.pi)
    
    # Sort by anchor angle
    residues_sorted = sorted(residues, key=lambda r: r.get("anchor_angle", 0.0))
    
    # Cluster into sectors
    sector_size = (2 * np.pi) / num_sectors
    clusters = []
    current_cluster = []
    current_sector = 0
    
    for res in residues_sorted:
        sector = int(res["anchor_angle"] / sector_size)
        
        if sector == current_sector:
            current_cluster.append(res)
        else:
            if current_cluster:
                clusters.append({
                    "sector": current_sector,
                    "residues": current_cluster,
                })
            current_cluster = [res]
            current_sector = sector
    
    if current_cluster:
        clusters.append({
            "sector": current_sector,
            "residues": current_cluster,
        })
    
    return clusters


def _allocate_arc_ranges(clusters: List[Dict], min_gap: float) -> List[Dict]:
    """Allocate arc ranges per cluster with gaps."""
    if not clusters:
        return []
    
    total_angle = 2 * np.pi
    total_gaps = (len(clusters) - 1) * min_gap
    available_angle = total_angle - total_gaps
    
    total_residues = sum(len(c["residues"]) for c in clusters)
    
    current_angle = 0.0
    allocated = []
    
    for idx, cluster in enumerate(clusters):
        cluster_weight = len(cluster["residues"]) / max(total_residues, 1)
        arc_size = available_angle * cluster_weight
        
        allocated.append({
            **cluster,
            "arcStart": current_angle,
            "arcEnd": current_angle + arc_size,
            "arcSize": arc_size,
        })
        
        current_angle += arc_size + (min_gap if idx < len(clusters) - 1 else 0)
    
    return allocated


def _resolve_collisions(
    residues: List[Dict],
    ligand_center: Tuple[float, float],
    min_radius: float,
    max_radius: float
) -> List[Dict]:
    """Resolve label overlaps with iterative relaxation."""
    resolved = [dict(r) for r in residues]
    max_iterations = 50
    push_factor = 1.1
    min_angular_sep = 0.12
    
    for _ in range(max_iterations):
        has_overlap = False
        resolved.sort(key=lambda r: r.get("angle", 0.0))
        
        for i in range(len(resolved)):
            for j in range(i + 1, len(resolved)):
                r1 = resolved[i]
                r2 = resolved[j]
                
                angle_diff = abs(r2.get("angle", 0.0) - r1.get("angle", 0.0))
                normalized_diff = min(angle_diff, 2 * np.pi - angle_diff)
                
                if normalized_diff < min_angular_sep:
                    radius1 = r1.get("radius", min_radius)
                    radius2 = r2.get("radius", min_radius)
                    
                    if radius1 < radius2:
                        r1["radius"] = min(radius1 * push_factor, max_radius)
                        r1["x"] = ligand_center[0] + r1["radius"] * np.cos(r1.get("angle", 0.0))
                        r1["y"] = ligand_center[1] + r1["radius"] * np.sin(r1.get("angle", 0.0))
                    else:
                        r2["radius"] = min(radius2 * push_factor, max_radius)
                        r2["x"] = ligand_center[0] + r2["radius"] * np.cos(r2.get("angle", 0.0))
                        r2["y"] = ligand_center[1] + r2["radius"] * np.sin(r2.get("angle", 0.0))
                    
                    has_overlap = True
        
        if not has_overlap:
            break
    
    return resolved
