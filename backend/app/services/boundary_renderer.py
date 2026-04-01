"""
Boundary Renderer
Computes smooth irregular circular pocket boundary around origin (0,0).
"""

import numpy as np
from typing import List, Dict, Tuple
from scipy.spatial import ConvexHull
import re

# Boundary expansion margin
BOUNDARY_MARGIN = 25.0  # pixels


def compute_offset_normals(points: np.ndarray) -> np.ndarray:
    """
    Compute outward-pointing normals for each edge of the polygon.
    
    Args:
        points: Array of 2D points forming a closed polygon
    
    Returns:
        Array of normalized outward normals
    """
    n = len(points)
    normals = []
    
    for i in range(n):
        p1 = points[i]
        p2 = points[(i + 1) % n]
        
        # Edge vector
        edge = p2 - p1
        
        # Perpendicular (rotated 90 degrees counterclockwise for outward)
        normal = np.array([-edge[1], edge[0]], dtype=float)
        norm = np.linalg.norm(normal) + 1e-9
        normal = normal / norm
        
        normals.append(normal)
    
    return np.array(normals, dtype=float)


def expand_hull_outward(points: np.ndarray, margin: float = BOUNDARY_MARGIN) -> np.ndarray:
    """
    Expand convex hull outward by margin using offset normals.
    
    Args:
        points: Array of 2D points from convex hull
        margin: Expansion margin in pixels
    
    Returns:
        Expanded points array
    """
    if len(points) < 3:
        return points
    
    normals = compute_offset_normals(points)
    
    # Expand each point outward
    expanded = []
    for i in range(len(points)):
        # Use average of adjacent normals for smoother expansion
        prev_normal = normals[(i - 1) % len(normals)]
        curr_normal = normals[i]
        avg_normal = (prev_normal + curr_normal) / 2
        avg_normal = avg_normal / (np.linalg.norm(avg_normal) + 1e-9)
        
        expanded_point = points[i] + margin * avg_normal
        expanded.append(expanded_point)
    
    return np.array(expanded, dtype=float)


def chaikin_smooth(points: np.ndarray, iterations: int = 4) -> np.ndarray:
    """
    Apply Chaikin's corner-cutting algorithm for smooth curves.
    
    Args:
        points: Array of 2D points [[x, y], ...]
        iterations: Number of smoothing iterations (3-5 recommended)
    
    Returns:
        Smoothed points array
    """
    if len(points) < 3:
        return points
    
    smoothed = points.copy()
    
    for _ in range(iterations):
        new_points = []
        
        for i in range(len(smoothed)):
            p1 = smoothed[i]
            p2 = smoothed[(i + 1) % len(smoothed)]
            
            # Chaikin: 1/4 and 3/4 points
            q1 = 0.75 * p1 + 0.25 * p2
            q2 = 0.25 * p1 + 0.75 * p2
            
            new_points.append(q1)
            new_points.append(q2)
        
        smoothed = np.array(new_points, dtype=float)
    
    return smoothed


def compute_pocket_boundary_from_ring_radius(
    max_ring_radius: float,
    margin: float = BOUNDARY_MARGIN
) -> str:
    """
    Compute smooth irregular circular pocket boundary around origin (0,0).
    
    Uses polar function: r(theta) = R + A1*sin(3θ) + A2*cos(5θ)
    This ensures pocket surrounds ligand at origin.
    
    Args:
        max_ring_radius: Maximum ring radius from sector layout
        margin: Expansion margin in pixels
    
    Returns:
        SVG path string for closed smooth boundary around origin
    """
    # Base radius = max ring radius + margin
    R = max_ring_radius + margin
    
    # Amplitude parameters for irregular shape
    A1 = R * 0.15  # 15% variation
    A2 = R * 0.10  # 10% variation
    
    # Generate smooth irregular circular contour
    num_points = 64  # Number of points for smooth curve
    theta_values = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
    
    # Compute radius for each angle: r(theta) = R + A1*sin(3θ) + A2*cos(5θ)
    r_values = R + A1 * np.sin(3 * theta_values) + A2 * np.cos(5 * theta_values)
    
    # Convert to Cartesian coordinates around origin (0,0)
    points = []
    for i, theta in enumerate(theta_values):
        r = r_values[i]
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        points.append([x, y])
    
    points = np.array(points, dtype=float)
    
    # Apply Chaikin smoothing for extra smoothness
    smoothed_points = chaikin_smooth(points, iterations=2)
    
    # Create SVG path
    path_parts = [f"M {smoothed_points[0][0]:.1f},{smoothed_points[0][1]:.1f}"]
    
    for i in range(1, len(smoothed_points)):
        path_parts.append(f"L {smoothed_points[i][0]:.1f},{smoothed_points[i][1]:.1f}")
    
    # Close path
    path_parts.append("Z")
    
    return " ".join(path_parts)


def translate_svg_path(path: str, dx: float, dy: float) -> str:
    """
    Translate SVG path by (dx, dy).
    
    Args:
        path: SVG path string
        dx: X translation
        dy: Y translation
    
    Returns:
        Translated SVG path string
    """
    if not path:
        return path
    
    # Parse and translate coordinates in path
    # Match numbers (including negative) in path commands
    def translate_coords(match):
        num = float(match.group(0))
        # Determine if this is x or y coordinate based on position
        # Simple approach: translate all numbers
        # More sophisticated: parse path commands
        return str(num)
    
    # Simple translation: find all coordinate pairs and translate
    # Pattern: M x,y or L x,y or C x1,y1 x2,y2 x3,y3 or Q x1,y1 x2,y2
    def translate_pair(match):
        x = float(match.group(1))
        y = float(match.group(2))
        return f"{x + dx:.1f},{y + dy:.1f}"
    
    # Replace coordinate pairs
    pattern = r'(-?\d+\.?\d*),(-?\d+\.?\d*)'
    translated = re.sub(pattern, translate_pair, path)
    
    return translated


def compute_pocket_boundary(
    residue_positions: List[np.ndarray],
    iterations: int = 4,
    margin: float = BOUNDARY_MARGIN
) -> str:
    """
    Legacy function: Compute boundary from residue positions.
    For unified coordinate system, use compute_pocket_boundary_from_ring_radius instead.
    """
    if not residue_positions:
        return ""
    
    # Compute max radius from residue positions
    max_radius = max([np.linalg.norm(pos) for pos in residue_positions], default=350.0)
    return compute_pocket_boundary_from_ring_radius(max_radius, margin)
