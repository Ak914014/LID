"""
Radial Layout Engine
Implements hybrid deterministic rings + force relaxation layout.
"""

import math
import numpy as np
from typing import List, Dict, Tuple

EPS = 1e-9

# Fixed ring radii (in pixels)
RING_RADIUS = {
    1: 170.0,  # Inner ring: H-bond, salt bridge, metal coordination
    2: 240.0,  # Middle ring: pi-stack, pi-cation
    3: 310.0,  # Outer ring: hydrophobic, others
}

# Node footprint radius for collision detection
NODE_FOOTPRINT = 20.0
MIN_ANGULAR_SPACING_PADDING = 8.0


def wrap_angle(a: float) -> float:
    """Wrap angle to [-π, π]."""
    while a <= -math.pi:
        a += 2 * math.pi
    while a > math.pi:
        a -= 2 * math.pi
    return a


def compute_pca_basis(ligand_atoms_3d: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute PCA basis from ligand heavy atom positions.
    Returns two orthonormal vectors (u, w) for 2D projection.
    """
    if ligand_atoms_3d.shape[0] < 2:
        # Fallback to standard basis
        return np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
    
    # Center the atoms
    centered = ligand_atoms_3d - ligand_atoms_3d.mean(axis=0)
    
    # Compute covariance matrix
    cov = np.cov(centered.T)
    
    # Get eigenvectors (principal components)
    eigenvals, eigenvecs = np.linalg.eigh(cov)
    
    # Sort by eigenvalue (descending)
    idx = np.argsort(eigenvals)[::-1]
    eigenvecs = eigenvecs[:, idx]
    
    # First two principal components
    u = eigenvecs[:, 0]
    w = eigenvecs[:, 1] if eigenvecs.shape[1] > 1 else np.array([0.0, 1.0, 0.0])
    
    # Ensure orthonormal
    u = u / (np.linalg.norm(u) + EPS)
    w = w - np.dot(w, u) * u  # Gram-Schmidt
    w = w / (np.linalg.norm(w) + EPS)
    
    return u, w


def deterministic_ring_layout(
    nodes: List[Dict],
    ligand_centroid: np.ndarray,
    basis_u: np.ndarray,
    basis_w: np.ndarray,
    ring_radius: Dict[int, float] = RING_RADIUS
) -> List[Dict]:
    """
    Deterministic initial placement on rings.
    
    Args:
        nodes: List of node dicts with "centroid3", "ring", "id"
        ligand_centroid: 3D position of ligand center
        basis_u, basis_w: Orthonormal basis vectors for 2D projection
        ring_radius: Dict mapping ring number to radius
    
    Returns:
        Updated nodes with "theta", "alpha", "pos2" fields
    """
    # Compute initial theta from 3D ordering (for sorting only)
    for n in nodes:
        if n["type"] == "ligand":
            n["theta"] = 0.0
            n["pos2"] = np.array([0.0, 0.0], dtype=float)
            continue
        
        v = n["centroid3"] - ligand_centroid
        x = float(np.dot(v, basis_u))
        y = float(np.dot(v, basis_w))
        n["theta"] = math.atan2(y, x)
    
    # Group by ring and assign evenly spaced angles
    for k in [1, 2, 3]:
        ring_nodes = [n for n in nodes if n.get("ring") == k]
        if not ring_nodes:
            continue
        
        # Sort by theta to preserve biological orientation
        ring_nodes.sort(key=lambda t: t["theta"])
        
        m = len(ring_nodes)
        if m == 0:
            continue
        
        # Stable phase anchored to smallest theta
        phi = ring_nodes[0]["theta"]
        
        # Evenly distribute around circle
        for j, n in enumerate(ring_nodes):
            n["alpha"] = phi + 2 * math.pi * j / m
    
    # Compute 2D positions on rings
    for n in nodes:
        if n["type"] == "ligand":
            continue
        
        R = ring_radius[n["ring"]]
        a = n["alpha"]
        n["pos2"] = np.array([R * math.cos(a), R * math.sin(a)], dtype=float)
    
    return nodes


def hybrid_relax(
    nodes: List[Dict],
    ring_radius: Dict[int, float] = RING_RADIUS,
    k_rep: float = 500.0,
    k_spr: float = 1.2,
    k_tan: float = 0.2,
    dt: float = 0.08,
    damp: float = 0.9,
    steps: int = 120
) -> List[Dict]:
    """
    Relax positions with constrained forces while maintaining ring structure.
    
    Args:
        nodes: List of node dicts with "pos2", "ring", "theta"
        ring_radius: Dict mapping ring number to radius
        k_rep: Repulsion force constant
        k_spr: Spring force constant (keeps nodes on ring)
        k_tan: Tangential force constant (preserves order)
        dt: Time step
        damp: Damping factor
        steps: Number of iterations
    
    Returns:
        Updated nodes with relaxed "pos2"
    """
    # Initialize velocities
    for n in nodes:
        if n["type"] == "ligand":
            continue
        n["vel"] = np.zeros(2, dtype=float)
    
    for _ in range(steps):
        # Compute forces for each node
        for i, ni in enumerate(nodes):
            if ni["type"] == "ligand":
                continue
            
            pi = ni["pos2"]
            Fi = np.zeros(2, dtype=float)
            
            # Repulsion from other nodes
            for j, nj in enumerate(nodes):
                if i == j or nj["type"] == "ligand":
                    continue
                pj = nj["pos2"]
                d = pi - pj
                r2 = float(np.dot(d, d)) + EPS
                Fi += k_rep * d / r2
            
            # Ring spring (radial constraint)
            R = ring_radius[ni["ring"]]
            r = float(np.linalg.norm(pi)) + EPS
            radial = pi / r
            Fi += -k_spr * (r - R) * radial
            
            # Tangential angle preservation (weak)
            theta_target = ni.get("theta", 0.0)
            theta_now = math.atan2(pi[1], pi[0])
            delta = wrap_angle(theta_target - theta_now)
            tang = np.array([-math.sin(theta_now), math.cos(theta_now)], dtype=float)
            Fi += k_tan * delta * tang
            
            # Integrate
            ni["vel"] = damp * ni["vel"] + dt * Fi
            ni["pos2"] = ni["pos2"] + dt * ni["vel"]
            
            # Hard project back to ring (stabilizing)
            r = float(np.linalg.norm(ni["pos2"])) + EPS
            ni["pos2"] = ring_radius[ni["ring"]] * ni["pos2"] / r
    
    return nodes


def compute_radial_layout(
    interaction_graph: Dict,
    ligand_atoms_3d: np.ndarray = None
) -> List[Dict]:
    """
    Main layout function: deterministic init + hybrid relaxation.
    
    Args:
        interaction_graph: Graph dict with "nodes", "edges", "ligand_centroid"
        ligand_atoms_3d: Optional 3D positions of ligand atoms for PCA
    
    Returns:
        List of nodes with "pos2" (2D positions) and updated angles
    """
    nodes = interaction_graph["nodes"]
    ligand_centroid = interaction_graph["ligand_centroid"]
    
    # Compute PCA basis from ligand geometry
    if ligand_atoms_3d is not None and ligand_atoms_3d.shape[0] > 1:
        basis_u, basis_w = compute_pca_basis(ligand_atoms_3d)
    else:
        # Fallback to standard basis
        basis_u = np.array([1.0, 0.0, 0.0])
        basis_w = np.array([0.0, 1.0, 0.0])
    
    # Step 1: Deterministic initial placement
    nodes = deterministic_ring_layout(nodes, ligand_centroid, basis_u, basis_w)
    
    # Step 2: Hybrid relaxation
    nodes = hybrid_relax(nodes, steps=120)
    
    # Step 3: Final collision resolution (angular spacing enforcement)
    nodes = resolve_angular_collisions(nodes)
    
    return nodes


def resolve_angular_collisions(nodes: List[Dict]) -> List[Dict]:
    """
    Enforce minimum angular spacing within each ring.
    """
    # Group by ring
    for k in [1, 2, 3]:
        ring_nodes = [n for n in nodes if n.get("ring") == k]
        if len(ring_nodes) < 2:
            continue
        
        # Sort by current angle
        ring_nodes.sort(key=lambda n: math.atan2(n["pos2"][1], n["pos2"][0]))
        
        R = RING_RADIUS[k]
        d_min = 2 * NODE_FOOTPRINT + MIN_ANGULAR_SPACING_PADDING
        delta_alpha_min = 2 * math.asin(d_min / (2 * R + EPS))
        
        # Walk around ring and enforce spacing
        for i in range(len(ring_nodes) - 1):
            n1 = ring_nodes[i]
            n2 = ring_nodes[i + 1]
            
            a1 = math.atan2(n1["pos2"][1], n1["pos2"][0])
            a2 = math.atan2(n2["pos2"][1], n2["pos2"][0])
            
            # Normalize angles
            if a2 < a1:
                a2 += 2 * math.pi
            
            delta = a2 - a1
            if delta < delta_alpha_min:
                # Push n2 away
                a2_new = a1 + delta_alpha_min
                n2["pos2"] = np.array([
                    R * math.cos(a2_new),
                    R * math.sin(a2_new)
                ], dtype=float)
        
        # Normalize drift (center the ring)
        angles = [math.atan2(n["pos2"][1], n["pos2"][0]) for n in ring_nodes]
        mean_angle = sum(angles) / len(angles)
        for n in ring_nodes:
            a = math.atan2(n["pos2"][1], n["pos2"][0])
            a_shifted = wrap_angle(a - mean_angle)
            n["pos2"] = np.array([
                R * math.cos(a_shifted),
                R * math.sin(a_shifted)
            ], dtype=float)
    
    return nodes
