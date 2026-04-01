"""
Interaction Graph Builder
Constructs a graph representation of protein-ligand interactions.
"""

from typing import List, Dict, Set, Tuple
import numpy as np
import MDAnalysis as mda


def build_interaction_graph(
    u: mda.Universe,
    ligand_resname: str,
    interactions: List[Dict]
) -> Dict:
    """
    Build interaction graph from detected interactions.
    
    Returns:
        {
            "nodes": [
                {"id": "ligand", "type": "ligand", "centroid3": np.array([x,y,z])},
                {"id": "RES123", "type": "residue", "resname": "ARG", "resid": 123, 
                 "chain": "A", "centroid3": np.array([x,y,z]), "ring": 1, "class": "positive"}
            ],
            "edges": [
                {"source": "ligand", "target": "RES123", "type": "hbond", "backbone": False}
            ]
        }
    """
    lig = u.select_atoms(f"resname {ligand_resname}")
    if lig.n_atoms == 0:
        raise ValueError(f"No ligand found with resname '{ligand_resname}'")
    
    ligand_centroid = lig.center_of_mass()
    
    # Create ligand node
    nodes = [{
        "id": "ligand",
        "type": "ligand",
        "centroid3": ligand_centroid.copy(),
    }]
    
    # Map interaction types to rings
    # Ring 1 (inner): H-bond + salt bridge (strong interactions)
    # Ring 2 (middle): pi-stack + pi-cation (medium interactions)
    # Ring 3 (outer): hydrophobic + others (weak interactions)
    ring_map = {
        "hbond": 1,
        "salt_bridge": 1,
        "pi_pi": 2,
        "pi_cation": 2,
        "hydrophobic": 3,
        "metal_coordination": 1,  # Strong interaction
        "halogen_bond": 1,  # Strong interaction
        "distance": 3,  # Generic/weak
    }
    
    # Group interactions by residue
    residue_interactions = {}
    for it in interactions:
        res_id = f"{it['resname']}{it['resid']}"
        if res_id not in residue_interactions:
            residue_interactions[res_id] = []
        residue_interactions[res_id].append(it)
    
    # Create residue nodes
    residue_nodes = {}
    edges = []
    
    for res_id, res_interactions in residue_interactions.items():
        # Get residue from universe
        first_it = res_interactions[0]
        resname = first_it['resname']
        resid = first_it['resid']
        chain = first_it.get('chain', '')
        
        # Try to select residue with chain information if available
        if chain and chain.strip():
            # Try with chain ID first
            res_atoms = u.select_atoms(
                f"resname {resname} and resid {resid} and segid {chain}"
            )
            if res_atoms.n_atoms == 0:
                # Fallback: try without chain (in case chain is stored differently)
                res_atoms = u.select_atoms(
                    f"resname {resname} and resid {resid}"
                )
        else:
            # No chain info, try without
            res_atoms = u.select_atoms(
                f"resname {resname} and resid {resid}"
            )
        
        if res_atoms.n_atoms == 0:
            # Last resort: try to find by resname and resid only (ignore chain)
            # This handles cases where chain info might be inconsistent
            all_residues = u.residues
            matching_res = None
            for res in all_residues:
                if res.resname == resname and int(res.resid) == int(resid):
                    matching_res = res
                    break
            
            if matching_res is None:
                # Still can't find it - skip this residue
                continue
            else:
                # Found it, get its atoms
                res_atoms = matching_res.atoms
        
        # Compute residue centroid (mean of contacting atoms or CA)
        ca_atoms = res_atoms.select_atoms("name CA")
        if ca_atoms.n_atoms > 0:
            centroid = ca_atoms.center_of_mass()
        else:
            centroid = res_atoms.center_of_mass()
        
        # Determine ring based on strongest interaction
        rings = [ring_map.get(it["type"], 3) for it in res_interactions]
        assigned_ring = min(rings)  # Use strongest (lowest ring number)
        
        # Create residue node
        node = {
            "id": res_id,
            "type": "residue",
            "resname": first_it["resname"],
            "resid": first_it["resid"],
            "chain": first_it.get("chain", "A"),
            "centroid3": centroid.copy(),
            "ring": assigned_ring,
            "class": first_it.get("res_class", "other"),
        }
        nodes.append(node)
        residue_nodes[res_id] = node
        
        # Create edges
        for it in res_interactions:
            edges.append({
                "source": "ligand",
                "target": res_id,
                "type": it["type"],
                "backbone": it.get("backbone", False),
                "distance": it.get("distance", 0.0),
                "ligand_atom_index": it.get("ligand_atom_index", -1),
            })
    
    return {
        "nodes": nodes,
        "edges": edges,
        "ligand_centroid": ligand_centroid,
    }
