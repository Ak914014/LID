"""
Interaction Detection
Implements geometric detection rules for protein-ligand interactions.
"""

import numpy as np
import warnings
from typing import List, Dict, Optional
import MDAnalysis as mda

from .residue_classifier import (
    classify_residue,
    is_backbone_atom,
    is_aromatic_residue,
    is_aromatic_atom,
    POSITIVE,
    NEGATIVE,
    HYDROPHOBIC,
)


def compute_angle_degrees(v1: np.ndarray, v2: np.ndarray) -> float:
    """Compute angle between two vectors in degrees."""
    v1_norm = v1 / (np.linalg.norm(v1) + 1e-10)
    v2_norm = v2 / (np.linalg.norm(v2) + 1e-10)
    dot = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
    angle_rad = np.arccos(dot)
    return np.degrees(angle_rad)


def detect_hbond(
    donor_pos: np.ndarray,
    acceptor_pos: np.ndarray,
    donor_h_pos: Optional[np.ndarray] = None,
    distance_cutoff: float = 3.5,
    angle_cutoff: float = 120.0
) -> bool:
    """
    Detect hydrogen bond using geometric rules.
    
    Args:
        donor_pos: Donor atom position (N or O)
        acceptor_pos: Acceptor atom position (O or N)
        donor_h_pos: Hydrogen position (if available)
        distance_cutoff: Maximum distance in Å
        angle_cutoff: Minimum angle in degrees
    
    Returns:
        True if H-bond detected
    """
    dist = np.linalg.norm(acceptor_pos - donor_pos)
    
    if dist > distance_cutoff:
        return False
    
    # If hydrogen position available, check angle
    if donor_h_pos is not None:
        vec_donor_acceptor = acceptor_pos - donor_pos
        vec_donor_h = donor_h_pos - donor_pos
        angle = compute_angle_degrees(vec_donor_acceptor, vec_donor_h)
        return angle >= angle_cutoff
    
    # Without H position, use distance only (simplified)
    return dist <= distance_cutoff


def detect_salt_bridge(
    pos_res_pos: np.ndarray,
    neg_res_pos: np.ndarray,
    distance_cutoff: float = 4.0
) -> bool:
    """Detect salt bridge between opposite charges."""
    dist = np.linalg.norm(neg_res_pos - pos_res_pos)
    return dist <= distance_cutoff


def detect_pi_pi_stacking(
    aromatic_centroid1: np.ndarray,
    aromatic_centroid2: np.ndarray,
    min_dist: float = 4.5,
    max_dist: float = 6.0
) -> bool:
    """Detect pi-pi stacking between aromatic systems."""
    dist = np.linalg.norm(aromatic_centroid2 - aromatic_centroid1)
    return min_dist <= dist <= max_dist


def detect_pi_cation(
    aromatic_centroid: np.ndarray,
    cation_pos: np.ndarray,
    distance_cutoff: float = 10.0
) -> bool:
    """Detect pi-cation interaction."""
    dist = np.linalg.norm(cation_pos - aromatic_centroid)
    return dist <= distance_cutoff


def get_aromatic_centroid(res_atoms, resname: str) -> Optional[np.ndarray]:
    """Compute centroid of aromatic ring atoms."""
    if not is_aromatic_residue(resname):
        return None
    
    aromatic_atom_names = {
        "CG", "CD1", "CD2", "CE1", "CE2", "CZ", 
        "CE3", "CZ2", "CZ3", "CH2"
    }
    
    aromatic_positions = []
    for atom in res_atoms:
        if atom.name.strip().upper() in aromatic_atom_names:
            aromatic_positions.append(atom.position)
    
    if len(aromatic_positions) >= 3:
        return np.mean(aromatic_positions, axis=0)
    return None


def detect_interactions(
    u: mda.Universe,
    mol,  # RDKit molecule
    ligand_resname: str,
    map_pdb_to_rd: Dict[int, int],
    cutoff_contact: float = 4.0
) -> List[Dict]:
    """
    Detect all protein-ligand interactions using geometric rules.
    
    Returns list of interaction dictionaries with type, distance, atoms, etc.
    """
    lig = u.select_atoms(f"resname {ligand_resname}")
    if lig.n_atoms == 0:
        return []
    
    # Only check residues within cutoff_contact distance from ligand
    # This makes the pocket_radius parameter actually affect the results
    protein = u.select_atoms(f"(protein or nucleic) and around {cutoff_contact} (resname {ligand_resname})")
    interactions = []
    
    for res in protein.residues:
        res_atoms = res.atoms
        if res_atoms.n_atoms == 0:
            continue
        
        # Compute distance matrix
        dmat = mda.lib.distances.distance_array(lig.positions, res_atoms.positions)
        if dmat.size == 0:
            continue
        
        min_idx = np.unravel_index(np.argmin(dmat), dmat.shape)
        lig_i = int(min_idx[0])
        res_j = int(min_idx[1])
        
        # Bounds check to prevent index out of bounds
        if lig_i >= lig.n_atoms or res_j >= res_atoms.n_atoms:
            continue
        
        min_dist = float(dmat[min_idx])
        
        # Initial distance check (will be refined by specific interaction types)
        if min_dist > cutoff_contact:
            continue
        
        # Additional bounds check before accessing atoms
        if lig_i < 0 or lig_i >= len(lig.atoms) or res_j < 0 or res_j >= len(res_atoms):
            continue
        
        lig_atom = lig.atoms[lig_i]
        prot_atom = res_atoms[res_j]
        prot_atom_name = prot_atom.name.strip()
        
        rd_idx = map_pdb_to_rd.get(lig_i, -1)
        rname = res.resname
        lig_elem = (lig_atom.element or lig_atom.name[0]).strip().upper()
        prot_elem = (prot_atom.element or prot_atom.name[0]).strip().upper()
        
        itype = None  # Will be set by interaction detection, None means no meaningful interaction
        
        # H-bond detection (geometric: distance ≤ 3.5 Å, angle > 120°)
        if min_dist <= 3.5 and lig_elem in {"N", "O"} and prot_elem in {"N", "O"}:
            # Try to find H atom for angle calculation
            donor_h_pos = None
            if lig_elem == "N":
                # Look for attached H (use atom position, not index in selection)
                try:
                    lig_atom_pos = lig_atom.position
                    h_atoms = lig.select_atoms(f"around 1.2 (point {lig_atom_pos[0]} {lig_atom_pos[1]} {lig_atom_pos[2]}) and name H*")
                    if h_atoms.n_atoms > 0:
                        donor_h_pos = h_atoms[0].position
                except Exception:
                    donor_h_pos = None
            elif prot_elem == "N":
                # Look for H attached to protein N
                try:
                    prot_atom_pos = prot_atom.position
                    h_atoms = res_atoms.select_atoms(f"around 1.2 (point {prot_atom_pos[0]} {prot_atom_pos[1]} {prot_atom_pos[2]}) and name H*")
                    if h_atoms.n_atoms > 0:
                        donor_h_pos = h_atoms[0].position
                except Exception:
                    donor_h_pos = None
            
            # Check both directions
            if detect_hbond(
                lig_atom.position,
                prot_atom.position,
                donor_h_pos,
                distance_cutoff=3.5,
                angle_cutoff=120.0
            ) or detect_hbond(
                prot_atom.position,
                lig_atom.position,
                donor_h_pos,
                distance_cutoff=3.5,
                angle_cutoff=120.0
            ):
                itype = "hbond"
        
        # Salt bridge (opposite charges ≤ 4.0 Å)
        elif min_dist <= 4.0:
            if (rname in POSITIVE and lig_elem in {"O", "N"}) or \
               (rname in NEGATIVE and lig_elem in {"N"}):
                itype = "salt_bridge"
        
        # Pi-Pi stacking (aromatic-aromatic, 4.5-6.0 Å)
        elif min_dist <= 6.0 and is_aromatic_residue(rname) and is_aromatic_atom(prot_atom_name, rname):
            if lig_elem == "C":
                # Check if ligand has aromatic system (simplified)
                aromatic_centroid1 = lig_atom.position  # Simplified
                aromatic_centroid2 = prot_atom.position
                if detect_pi_pi_stacking(aromatic_centroid1, aromatic_centroid2, 4.5, 6.0):
                    itype = "pi_pi"
        
        # Pi-cation (aromatic to cation ≤ 6 Å)
        elif min_dist <= 6.0 and rname in POSITIVE and lig_elem == "C":
            aromatic_centroid = lig_atom.position  # Simplified
            cation_pos = prot_atom.position
            if detect_pi_cation(aromatic_centroid, cation_pos, 6.0):
                itype = "pi_cation"
        
        # Metal coordination (metal ions ≤ 3.0 Å)
        elif min_dist <= 3.0:
            metal_elements = {"ZN", "FE", "MG", "CA", "MN", "CU", "NA", "K", "CO", "NI"}
            if lig_elem in metal_elements or prot_elem in metal_elements:
                itype = "metal_coordination"
        
        # Halogen bond (halogen to electron-rich ≤ 3.8 Å)
        elif min_dist <= 3.8 and lig_elem in {"F", "CL", "BR", "I"}:
            if prot_elem in {"N", "O", "S"} or is_aromatic_atom(prot_atom_name, rname):
                itype = "halogen_bond"
        
        # Hydrophobic contacts (heavy atoms ≤ 4.0 Å, strict cutoff)
        elif min_dist <= 4.0 and rname in HYDROPHOBIC and lig_elem in {"C", "S", "F", "CL", "BR", "I"}:
            itype = "hydrophobic"
        
        # Only add interactions of meaningful types (exclude generic "distance" contacts)
        meaningful_types = {
            "hbond", "salt_bridge", "pi_pi", "pi_cation", 
            "hydrophobic", "metal_coordination", "halogen_bond"
        }
        
        if itype is None or itype not in meaningful_types:
            continue  # Skip residues without meaningful interactions
        
        interactions.append({
            "residue": f"{rname}{int(res.resid)}",
            "resname": rname,
            "resid": int(res.resid),
            "chain": getattr(res, "segid", "") or "A",
            "res_class": classify_residue(rname),
            "type": itype,
            "distance": round(min_dist, 2),
            "ligand_atom_index": int(rd_idx),
            "protein_atom_name": prot_atom_name,
            "backbone": bool(is_backbone_atom(prot_atom_name)),
        })
    
    return interactions
