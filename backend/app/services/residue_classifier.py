"""
Residue Classifier
Classifies amino acid residues into biochemical classes for color coding.
"""

# Residue classification dictionaries
HYDROPHOBIC = {
    "ALA", "VAL", "LEU", "ILE", "PHE", "MET", "PRO", "TRP", "TYR"
}

POLAR = {
    "SER", "THR", "ASN", "GLN", "CYS"
}

POSITIVE = {
    "ARG", "LYS", "HIS"
}

NEGATIVE = {
    "ASP", "GLU"
}

NUCLEIC_ACID = {
    "A", "T", "G", "C", "U",  # RNA/DNA bases
    "DA", "DT", "DG", "DC", "DU",  # DNA
    "A5", "A3", "T5", "T3", "G5", "G3", "C5", "C3", "U5", "U3",  # Terminal
    "RA", "RU", "RG", "RC"  # RNA
}

METAL = {
    "ZN", "FE", "MG", "CA", "MN", "CU", "NA", "K", "CO", "NI"
}

WATER = {
    "HOH", "WAT", "TIP", "TIP3", "TIP4"
}


def classify_residue(resname: str) -> str:
    """
    Classify residue into biochemical class.
    
    Returns:
        "hydrophobic", "polar", "positive", "negative", "glycine", 
        "metal", "water", "nucleic", or "other"
    """
    resname_upper = resname.strip().upper()
    
    if resname_upper == "GLY":
        return "glycine"
    elif resname_upper in HYDROPHOBIC:
        return "hydrophobic"
    elif resname_upper in POLAR:
        return "polar"
    elif resname_upper in POSITIVE:
        return "positive"
    elif resname_upper in NEGATIVE:
        return "negative"
    elif resname_upper in METAL:
        return "metal"
    elif resname_upper in WATER:
        return "water"
    elif resname_upper in NUCLEIC_ACID:
        return "nucleic"
    else:
        return "other"


def get_residue_color(res_class: str) -> str:
    """
    Get color code for residue class.
    
    Returns hex color code.
    """
    color_map = {
        "hydrophobic": "#86efac",  # green
        "polar": "#7dd3fc",        # cyan
        "positive": "#3b82f6",     # blue
        "negative": "#fb923c",     # orange
        "glycine": "#fef08a",      # beige
        "metal": "#4b5563",        # dark grey
        "water": "#9ca3af",        # light grey
        "nucleic": "#4b5563",      # dark grey
        "other": "#9ca3af",        # light grey
    }
    return color_map.get(res_class, "#9ca3af")


def is_backbone_atom(atom_name: str) -> bool:
    """Check if atom is part of protein backbone."""
    backbone_atoms = {"N", "CA", "C", "O", "OXT"}
    return atom_name.strip().upper() in backbone_atoms


def is_aromatic_residue(resname: str) -> bool:
    """Check if residue has aromatic ring."""
    aromatic = {"PHE", "TYR", "TRP", "HIS"}
    return resname.strip().upper() in aromatic


def is_aromatic_atom(atom_name: str, resname: str) -> bool:
    """Check if atom is part of aromatic ring."""
    if not is_aromatic_residue(resname):
        return False
    aromatic_patterns = {
        "CG", "CD1", "CD2", "CE1", "CE2", "CZ", 
        "CE3", "CZ2", "CZ3", "CH2"
    }
    return atom_name.strip().upper() in aromatic_patterns
