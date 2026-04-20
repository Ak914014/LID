import logging
import os
import tempfile
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.optimize import linear_sum_assignment
from scipy.spatial import ConvexHull

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from app.services.ai_interaction_service import rank_and_filter_interactions

import MDAnalysis as mda
try:
    import prolif as plf
except ImportError:  # pragma: no cover - optional runtime dependency
    plf = None

logger = logging.getLogger(__name__)


# =========================
# Residue classification
# =========================

NEG_RES = {"ASP", "GLU"}
POS_RES = {"ARG", "LYS", "HIS"}
# Maestro LID commonly colors TYR with hydrophobic / aromatic greens (not polar blue).
HYDROPHOBIC_RES = {"VAL", "LEU", "ILE", "ALA", "MET", "PHE", "TRP", "PRO", "TYR"}
POLAR_RES = {"ASN", "GLN", "SER", "THR"}  # CYS classified separately (Maestro-style)

NUCLEIC_RES = {
    "A", "T", "G", "C", "U",
    "DA", "DT", "DG", "DC", "DU",
    "RA", "RU", "RG", "RC",
    "ADE", "THY", "GUA", "CYT", "URA",
}

# Standard protein residue names (for ligand auto-detection: ligand = first non-protein)
STANDARD_PROTEIN_RES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HOH", "WAT", "H2O", "NAG", "NDG", "MAN", "BMA", "CL", "NA", "CA", "ZN",
}


def _pdb_line_resname(line: str) -> str:
    """Parse residue name from ATOM/HETATM line; supports 3- and 4-char (columns 17-20)."""
    if len(line) < 20:
        return ""
    return (line[17:21] if len(line) >= 21 else line[17:20]).strip()


def _pdb_line_serial(line: str) -> int | None:
    """Parse atom serial from ATOM/HETATM line (columns 7-11)."""
    if len(line) < 11:
        return None
    try:
        return int(line[6:11].strip())
    except ValueError:
        return None


def _is_pdb_atomish_record(line: str) -> bool:
    rec = line[0:6].strip() if len(line) >= 6 else ""
    return rec in {"ATOM", "HETATM", "TER", "CONECT"}


def _parse_model_number(line: str, fallback: int) -> int:
    parts = line.split()
    if len(parts) >= 2 and parts[1].isdigit():
        return int(parts[1])
    return fallback


def _atom_identity_key(line: str) -> tuple[str, str, str, str, str]:
    rec = line[0:6].strip() if len(line) >= 6 else ""
    atom_name = line[12:16].strip() if len(line) >= 16 else ""
    resname = _pdb_line_resname(line)
    chain = line[21:22].strip() if len(line) >= 22 else ""
    resid = line[22:26].strip() if len(line) >= 26 else ""
    return rec, atom_name, resname, chain, resid


def _dedupe_outside_lines_keep_one_protein(outside_lines: list[str]) -> list[str]:
    """
    Some docking exports repeat the same protein between each MODEL block.
    Keep only one copy of repeated outside-model atoms.
    """
    seen_atom_keys: set[tuple[str, str, str, str, str]] = set()
    seen_misc_lines: set[str] = set()
    out: list[str] = []

    for ln in outside_lines:
        rec = ln[0:6].strip() if len(ln) >= 6 else ""
        if rec in {"ATOM", "HETATM"}:
            key = _atom_identity_key(ln)
            if key in seen_atom_keys:
                continue
            seen_atom_keys.add(key)
            out.append(ln)
            continue

        # TER/CONECT repeated with each copied protein block: keep first instance.
        if ln in seen_misc_lines:
            continue
        seen_misc_lines.add(ln)
        out.append(ln)
    return out


def extract_selected_model(pdb_path: str, model_index: int) -> str:
    """
    Returns a clean PDB string containing only one selected model plus
    relevant non-model protein content (if present outside MODEL blocks).
    """
    with open(pdb_path, "r", encoding="utf-8") as f:
        pdb_text = f.read()
    selected_pdb, _ = extract_selected_model_from_text(pdb_text, model_index=model_index)
    return selected_pdb


def extract_selected_model_from_text(pdb_text: str, model_index: int) -> tuple[str, int]:
    lines = pdb_text.splitlines()
    outside_lines: list[str] = []
    model_blocks: dict[int, list[str]] = {}

    in_model = False
    current_model = None
    fallback_model_num = 1

    for line in lines:
        if line.startswith("MODEL"):
            in_model = True
            current_model = _parse_model_number(line, fallback_model_num)
            fallback_model_num += 1
            model_blocks.setdefault(current_model, [])
            continue
        if line.startswith("ENDMDL"):
            in_model = False
            current_model = None
            continue

        if in_model and current_model is not None:
            if _is_pdb_atomish_record(line):
                model_blocks[current_model].append(line)
            continue

        if _is_pdb_atomish_record(line):
            outside_lines.append(line)

    if not model_blocks:
        # No MODEL records: treat as single-model input.
        cleaned = "\n".join(outside_lines) + "\nEND\n"
        return cleaned, 1

    if model_index not in model_blocks:
        available = ", ".join(str(k) for k in sorted(model_blocks.keys()))
        raise ValueError(f"Requested model_index={model_index} not found. Available models: {available}")

    selected_lines = model_blocks[model_index]
    first_model_index = sorted(model_blocks.keys())[0]
    first_model_lines = model_blocks[first_model_index]

    def _is_protein_atom_line(ln: str) -> bool:
        rec = ln[0:6].strip() if len(ln) >= 6 else ""
        if rec != "ATOM":
            return False
        return _pdb_line_resname(ln) in STANDARD_PROTEIN_RES

    # Split selected model lines into protein and non-protein parts
    selected_protein = [ln for ln in selected_lines if _is_protein_atom_line(ln)]
    selected_non_protein = [ln for ln in selected_lines if not _is_protein_atom_line(ln)]
    first_model_protein = [ln for ln in first_model_lines if _is_protein_atom_line(ln)]

    out: list[str] = []
    if selected_protein:
        # Full complex in this MODEL (NMR/MD ensembles, flexible receptor, etc.):
        # receptor and ligand coordinates must both come from the selected pose.
        out.extend(selected_lines)
    elif first_model_protein and selected_non_protein:
        # Rigid-receptor docking: MODEL k may list only the ligand pose; keep protein from model 1.
        out.extend(first_model_protein)
        out.extend(selected_non_protein)
    else:
        # Protein outside MODEL records, or ligand-only file layout
        out.extend(_dedupe_outside_lines_keep_one_protein(outside_lines))
        out.extend(selected_lines)

    out.append("END")
    return "\n".join(out) + "\n", len(model_blocks)


def extract_single_protein_all_models_from_text(pdb_text: str) -> tuple[str, int]:
    """
    Build a visualization-friendly PDB with:
    - one protein scaffold (from first model, when duplicated per MODEL)
    - all model-specific non-protein atoms kept in MODEL/ENDMDL blocks
    """
    lines = pdb_text.splitlines()
    model_blocks: dict[int, list[str]] = {}
    outside_lines: list[str] = []
    in_model = False
    current_model = None
    fallback_model_num = 1

    for line in lines:
        if line.startswith("MODEL"):
            in_model = True
            current_model = _parse_model_number(line, fallback_model_num)
            fallback_model_num += 1
            model_blocks.setdefault(current_model, [])
            continue
        if line.startswith("ENDMDL"):
            in_model = False
            current_model = None
            continue

        if in_model and current_model is not None:
            if _is_pdb_atomish_record(line):
                model_blocks[current_model].append(line)
            continue

        if _is_pdb_atomish_record(line):
            outside_lines.append(line)

    if not model_blocks:
        return (pdb_text if pdb_text.endswith("\n") else pdb_text + "\n"), 1

    def _is_protein_atom_line(ln: str) -> bool:
        rec = ln[0:6].strip() if len(ln) >= 6 else ""
        return rec == "ATOM" and _pdb_line_resname(ln) in STANDARD_PROTEIN_RES

    model_ids = sorted(model_blocks.keys())
    first_model_lines = model_blocks[model_ids[0]]
    first_model_protein = [ln for ln in first_model_lines if _is_protein_atom_line(ln)]

    out: list[str] = []
    if first_model_protein:
        out.extend(first_model_protein)
        for mid in model_ids:
            lig_lines = [ln for ln in model_blocks[mid] if not _is_protein_atom_line(ln)]
            out.append(f"MODEL {mid}")
            out.extend(lig_lines)
            out.append("ENDMDL")
        out.append("END")
        return "\n".join(out) + "\n", len(model_ids)

    # Fallback for "protein outside MODEL + ligand in MODEL" layout.
    out.extend(_dedupe_outside_lines_keep_one_protein(outside_lines))
    for mid in model_ids:
        out.append(f"MODEL {mid}")
        out.extend(model_blocks[mid])
        out.append("ENDMDL")
    out.append("END")
    return "\n".join(out) + "\n", len(model_ids)


def classify_residue(resname: str) -> str:
    r = (resname or "").strip().upper()
    if r in NUCLEIC_RES:
        return "nucleic"
    if r == "GLY":
        return "glycine"
    if r == "CYS":
        return "cysteine"
    if r in NEG_RES:
        return "negative"
    if r in POS_RES:
        return "positive"
    if r in HYDROPHOBIC_RES:
        return "hydrophobic"
    if r in POLAR_RES:
        return "polar"
    return "other"


# =========================
# RDKit SVG + 2D coordinates
# =========================

def rdkit_svg(mol: Chem.Mol, w: int = 820, h: int = 300) -> str:
    if not mol.GetNumConformers():
        AllChem.Compute2DCoords(mol)
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def ligand_2d_coords(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    pts = []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        pts.append((p.x, p.y))
    return np.array(pts, dtype=float)


def rdkit_full_to_heavy_index_map(mol: Chem.Mol) -> dict[int, int]:
    """Map RDKit atom index (molecule possibly with H) to 0..N-1 heavy-atom-only indices."""
    m: dict[int, int] = {}
    j = 0
    for i in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(i).GetAtomicNum() != 1:
            m[i] = j
            j += 1
    return m


def compute_ligand_atom_solvent_exposed(
    u: mda.Universe,
    ligand_resname: str,
    mol: Chem.Mol,
    mol_draw: Chem.Mol,
    min_protein_dist_angstrom: float = 3.45,
) -> list[bool]:
    """
    Per mol_draw atom (same order as ligand_atom_xy): True if that ligand heavy atom is
    not in tight contact with protein (Maestro-style solvent-exposed marker dots).

    Uses minimum heavy-atom distance to protein; atoms with min_dist > threshold are "exposed".
    """
    n_draw = mol_draw.GetNumAtoms()
    if n_draw == 0:
        return []

    lig = u.select_atoms(f"resname {ligand_resname}")
    prot = u.select_atoms("protein")
    if lig.n_atoms == 0 or prot.n_atoms == 0:
        return [False] * n_draw

    dmat = mda.lib.distances.distance_array(lig.positions, prot.positions)
    min_per_pdb = np.min(dmat, axis=1)

    try:
        pdb_to_rd = map_pdb_ligand_atoms_to_rdkit(u, mol, ligand_resname)
    except Exception:
        return [False] * n_draw

    rd_to_pdb = {int(v): int(k) for k, v in pdb_to_rd.items()}
    full_to_heavy = rdkit_full_to_heavy_index_map(mol)
    heavy_to_full = {int(hi): int(fi) for fi, hi in full_to_heavy.items()}

    out: list[bool] = []
    for h in range(n_draw):
        full_i = heavy_to_full.get(h)
        if full_i is None:
            out.append(False)
            continue
        pdb_i = rd_to_pdb.get(full_i)
        if pdb_i is None or pdb_i < 0 or pdb_i >= len(min_per_pdb):
            out.append(False)
            continue
        md = float(min_per_pdb[pdb_i])
        out.append(md > float(min_protein_dist_angstrom))
    return out


def normalize_to_svg_coords(xy: np.ndarray, svg_w: int, svg_h: int, pad: int = 160) -> np.ndarray:
    """
    Increased pad shrinks ligand so it fits more nicely inside pocket.
    """
    xy = np.array(xy, dtype=float)
    minx, miny = xy.min(axis=0)
    maxx, maxy = xy.max(axis=0)
    spanx = max(maxx - minx, 1e-6)
    spany = max(maxy - miny, 1e-6)

    scale = min((svg_w - 2 * pad) / spanx, (svg_h - 2 * pad) / spany)

    out = np.zeros_like(xy)
    out[:, 0] = pad + (xy[:, 0] - minx) * scale
    out[:, 1] = svg_h - (pad + (xy[:, 1] - miny) * scale)  # invert y for SVG
    return out


# =========================
# Universe helper
# =========================

def universe_from_pdb_bytes(pdb_bytes: bytes) -> mda.Universe:
    normalized = _normalize_pdb_for_universe(pdb_bytes)
    fd, path = tempfile.mkstemp(suffix=".pdb")
    os.write(fd, normalized)
    os.close(fd)
    try:
        u = mda.Universe(path)
        return u
    finally:
        try:
            os.remove(path)
        except Exception:
            pass


def _normalize_pdb_for_universe(pdb_bytes: bytes) -> bytes:
    """
    Produce a single-frame PDB that MDAnalysis can load reliably.
    When file has MODEL 1 (e.g. ligand) then ENDMDL then ATOM lines (protein),
    merge into one block so we get one structure with both. Otherwise return as-is.
    """
    text = pdb_bytes.decode("utf-8")
    lines = text.splitlines()
    model_start = None
    endmdl_at = None
    for i, line in enumerate(lines):
        if line.startswith("MODEL "):
            model_start = i
            break
    if model_start is None:
        return pdb_bytes
    for i in range(model_start + 1, len(lines)):
        if lines[i].startswith("ENDMDL"):
            endmdl_at = i
            break
    if endmdl_at is None:
        return pdb_bytes
    has_atom_after = False
    for j in range(endmdl_at + 1, len(lines)):
        if lines[j].startswith("MODEL "):
            break
        if len(lines[j]) >= 22 and lines[j][0:6].strip() in ("ATOM", "HETATM"):
            has_atom_after = True
            break
    if not has_atom_after:
        return pdb_bytes
    out = []
    for i in range(model_start + 1, endmdl_at):
        line = lines[i]
        if len(line) >= 22 and line[0:6].strip() in ("ATOM", "HETATM"):
            out.append(line)
    n_ligand = len(out)
    first_serials = set()
    for i in range(model_start + 1, endmdl_at):
        line = lines[i]
        if len(line) >= 11 and line[0:6].strip() in ("ATOM", "HETATM"):
            try:
                first_serials.add(int(line[6:11].strip()))
            except ValueError:
                pass
    for line in lines:
        if line.startswith("CONECT") and len(line) >= 16:
            parts = line.split()
            if len(parts) >= 2 and parts[1].isdigit() and int(parts[1]) in first_serials:
                out.append(line)
    # Protein block: renumber serials so they continue after ligand (no duplicate serials)
    next_serial = n_ligand + 1
    for i in range(endmdl_at + 1, len(lines)):
        if lines[i].startswith("MODEL "):
            break
        line = lines[i]
        if len(line) >= 22 and line[0:6].strip() in ("ATOM", "HETATM"):
            # PDB serial is columns 7-11 (1-indexed), right-justified
            serial_str = str(next_serial).rjust(5)
            new_line = line[:6] + serial_str + line[11:]
            out.append(new_line)
            next_serial += 1
    if len(out) <= n_ligand:
        return pdb_bytes
    return ("\n".join(out) + "\n").encode("utf-8")


# =========================
# Ligand from PDB (when no SDF: PDB has ATOM/HETATM + CONECT = full ligand data)
# =========================

def _parse_pdb_ligand_block(pdb_text: str, ligand_resname: str) -> tuple[list[str], list[str], set[int]]:
    """
    Parse PDB text and return (atom_lines, conect_lines, ligand_serials).
    Scans the entire file so any PDB layout works (single model, multi-model,
    ligand first or last). Only includes ATOM/HETATM with resname == ligand_resname
    and CONECT that reference only those serials.
    """
    lines = pdb_text.splitlines()
    atom_lines: list[str] = []
    ligand_serials: set[int] = set()

    for line in lines:
        if len(line) < 22:
            continue
        rec = line[0:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        resname = _pdb_line_resname(line)
        if resname != ligand_resname:
            continue
        serial = _pdb_line_serial(line)
        if serial is None:
            continue
        ligand_serials.add(serial)
        atom_lines.append(line)

    # CONECT: include if central serial and all bonded serials are in ligand_serials
    conect_lines: list[str] = []
    for line in lines:
        if line.startswith("CONECT") and len(line) >= 16:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                ser = int(parts[1])
            except ValueError:
                continue
            if ser not in ligand_serials:
                continue
            others = [int(x) for x in parts[2:] if x.isdigit()]
            if all(o in ligand_serials for o in others):
                conect_lines.append(line)

    return atom_lines, conect_lines, ligand_serials


def extract_ligand_pdb_block(pdb_bytes: bytes, ligand_resname: str) -> str:
    """
    Extract a PDB block containing only the ligand (ATOM/HETATM + CONECT)
    so RDKit can build a molecule from it (no separate SDF needed).
    """
    pdb_text = pdb_bytes.decode("utf-8")
    atom_lines, conect_lines, _ = _parse_pdb_ligand_block(pdb_text, ligand_resname)
    if not atom_lines:
        raise ValueError(f"No ligand atoms found with resname '{ligand_resname}' in PDB.")
    block = "\n".join(atom_lines) + "\n"
    if conect_lines:
        block += "\n".join(conect_lines) + "\n"
    block += "END\n"
    return block


def detect_ligand_resname_from_pdb(pdb_bytes: bytes) -> str:
    """
    Detect ligand residue name: first non-protein residue found in the PDB
    (whole file, any MODEL). Works for all layouts: single model, multi-model,
    ligand before or after protein.
    """
    pdb_text = pdb_bytes.decode("utf-8")
    for line in pdb_text.splitlines():
        if len(line) < 20:
            continue
        rec = line[0:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        resname = _pdb_line_resname(line)
        if resname and resname not in STANDARD_PROTEIN_RES:
            return resname
    raise ValueError(
        "Could not detect ligand residue name in PDB. "
        "Provide ligand_resname (e.g. UNL, LIG) or use a PDB with a non-protein residue."
    )


def _resolve_single_ligand_resname(u: mda.Universe, ligand_resname: str | None) -> str:
    water_resnames = {"HOH", "WAT", "H2O"}
    if ligand_resname:
        residues = [r for r in u.select_atoms(f"resname {ligand_resname}").residues if r.resname not in water_resnames]
        if len(residues) == 0:
            raise ValueError(f"No ligand residue found for ligand_name '{ligand_resname}' in selected model.")
        if len(residues) > 1:
            raise ValueError(
                f"Expected exactly one ligand residue for ligand_name '{ligand_resname}', found {len(residues)}."
            )
        return ligand_resname

    candidates = [
        r for r in u.select_atoms("not protein and not (resname HOH WAT H2O)").residues
        if r.resname not in STANDARD_PROTEIN_RES
    ]
    if len(candidates) == 0:
        raise ValueError("No ligand residue detected in selected model.")
    if len(candidates) > 1:
        names = ", ".join(f"{r.resname}{int(r.resid)}" for r in candidates[:8])
        raise ValueError(f"Expected one ligand residue in selected model, found {len(candidates)}: {names}")
    return candidates[0].resname


def ligand_mol_from_pdb(pdb_bytes: bytes, ligand_resname: str) -> Chem.Mol:
    """
    Build an RDKit molecule from the ligand in the PDB (ATOM + CONECT).
    Use this when no separate SDF file is available.
    Uses sanitize=False because PDB/CONECT often produce valence warnings;
    2D depiction and coordinates still work.
    """
    block = extract_ligand_pdb_block(pdb_bytes, ligand_resname)
    mol = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError(
            f"Failed to build ligand molecule from PDB for resname '{ligand_resname}'. "
            "Check that the PDB has CONECT records for the ligand."
        )
    return mol


# =========================
# Ligand atom mapping (PDB ligand atoms -> RDKit atom indices)
# =========================

def _mol_get_3d_coords(mol: Chem.Mol) -> np.ndarray | None:
    if mol.GetNumConformers() == 0:
        return None
    conf = mol.GetConformer()
    pts = []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        pts.append((p.x, p.y, p.z))
    pts = np.array(pts, dtype=float)

    # if Z range ~0, treat as 2D-only
    if np.ptp(pts[:, 2]) < 1e-3:
        return None
    return pts


def _kabsch(P: np.ndarray, Q: np.ndarray):
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = Q.mean(axis=0) - (P.mean(axis=0) @ R)
    return R, t


def map_pdb_ligand_atoms_to_rdkit(u: mda.Universe, mol: Chem.Mol, ligand_resname: str) -> dict[int, int]:
    lig = u.select_atoms(f"resname {ligand_resname}")
    if lig.n_atoms == 0:
        raise ValueError(f"No ligand found with resname '{ligand_resname}' in PDB.")

    pdb_xyz = lig.positions.astype(float)
    pdb_elems = np.array([a.element.strip().upper() if a.element else a.name[0].upper() for a in lig.atoms])

    rd_xyz = _mol_get_3d_coords(mol)
    rd_elems = np.array([mol.GetAtomWithIdx(i).GetSymbol().upper() for i in range(mol.GetNumAtoms())])

    # If SDF has 3D coords with same atom count, try alignment
    if rd_xyz is not None and rd_xyz.shape[0] == pdb_xyz.shape[0]:
        R, t = _kabsch(rd_xyz, pdb_xyz)
        rd_aligned = rd_xyz @ R + t

        rmsd = float(np.sqrt(np.mean(np.sum((rd_aligned - pdb_xyz) ** 2, axis=1))))
        if rmsd < 2.0:
            return {i: i for i in range(lig.n_atoms)}

        cost = np.linalg.norm(rd_aligned[:, None, :] - pdb_xyz[None, :, :], axis=2)
        mismatch = (rd_elems[:, None] != pdb_elems[None, :])
        cost = cost + mismatch.astype(float) * 1e6

        row_ind, col_ind = linear_sum_assignment(cost)
        return {int(col): int(row) for row, col in zip(row_ind, col_ind)}

    # Fallback: element-based assignment using centered/scaled 2D coords
    AllChem.Compute2DCoords(mol)
    rd_xy = ligand_2d_coords(mol)
    rd_xyz2 = np.c_[rd_xy, np.zeros((rd_xy.shape[0],), dtype=float)]

    rd_xyz2c = rd_xyz2 - rd_xyz2.mean(axis=0)
    pdb_centered = pdb_xyz - pdb_xyz.mean(axis=0)

    s = (np.linalg.norm(pdb_centered) / max(np.linalg.norm(rd_xyz2c), 1e-6))
    rd_scaled = rd_xyz2c * s

    cost = np.linalg.norm(rd_scaled[:, None, :] - pdb_centered[None, :, :], axis=2)
    mismatch = (rd_elems[:, None] != pdb_elems[None, :])
    cost = cost + mismatch.astype(float) * 1e6

    row_ind, col_ind = linear_sum_assignment(cost)
    return {int(col): int(row) for row, col in zip(row_ind, col_ind)}


# =========================
# Pocket residue ring placement
# =========================

def place_residues_ring(u: mda.Universe, ligand_resname: str, pocket_radius: float):
    lig = u.select_atoms(f"resname {ligand_resname}")
    if lig.n_atoms == 0:
        raise ValueError(f"No ligand found with resname '{ligand_resname}' in PDB.")
    lig_c = lig.center_of_mass()

    pocket_sel = f"(protein or nucleic) and around {pocket_radius} (resname {ligand_resname})"
    pocket = u.select_atoms(pocket_sel)

    res_nodes = []
    for res in pocket.residues:
        rc = res.atoms.center_of_mass()
        v = rc - lig_c
        chain = getattr(res, "segid", "") or "A"
        res_nodes.append({
            "resname": res.resname,
            "resid": int(res.resid),
            "chain": chain,
            "dist": float(np.linalg.norm(v)),
            "class": classify_residue(res.resname),
        })

    if not res_nodes:
        return []

    # Sort by sequence order: chain first, then resid (small number first)
    res_nodes.sort(key=lambda n: (n["chain"], n["resid"]))

    # One node per residue (avoid duplicate MDAnalysis segments / alternate records)
    seen_keys: set[tuple[str, int, str]] = set()
    deduped: list[dict] = []
    for n in res_nodes:
        key = (str(n["chain"]), int(n["resid"]), str(n["resname"]))
        if key in seen_keys:
            continue
        seen_keys.add(key)
        deduped.append(n)
    res_nodes = deduped

    # Assign angles so that order around the ring = sequence order (210, 211, 212...)
    n = len(res_nodes)
    if n == 0:
        return []
    for i, node in enumerate(res_nodes):
        # Evenly spaced angles: 0, 2π/n, 4π/n, ... so first residue in sequence is at angle 0
        node["angle"] = 2.0 * np.pi * i / n

    return res_nodes


def _point_at_arc_fraction_closed(poly: np.ndarray, t: float) -> np.ndarray:
    """Point at fraction t ∈ [0,1) along closed polyline (arc length)."""
    M = int(poly.shape[0])
    if M < 2:
        return poly[0].astype(float).copy()
    seg_lens: list[float] = []
    for i in range(M):
        a = poly[i]
        b = poly[(i + 1) % M]
        seg_lens.append(float(np.hypot(b[0] - a[0], b[1] - a[1])))
    total = float(sum(seg_lens))
    if total < 1e-9:
        return poly[0].astype(float).copy()
    dist = (float(t) % 1.0) * total
    acc = 0.0
    for i in range(M):
        sl = seg_lens[i]
        if acc + sl >= dist - 1e-9:
            alpha = (dist - acc) / sl if sl > 1e-9 else 0.0
            a = poly[i]
            b = poly[(i + 1) % M]
            return a + alpha * (b - a)
        acc += sl
    return poly[0].astype(float).copy()


def _offset_closed_curve_radially(
    curve: np.ndarray,
    lig_center: np.ndarray,
    outward_px: float,
) -> np.ndarray:
    """Push each sample outward from ligand center (room for labels outside pocket band)."""
    L = np.asarray(lig_center, dtype=float).reshape(2)
    out = np.zeros_like(curve, dtype=float)
    d = float(outward_px)
    for i in range(curve.shape[0]):
        P = curve[i]
        v = P - L
        nv = float(np.linalg.norm(v))
        if nv < 1e-9:
            out[i] = P
        else:
            out[i] = P + (v / nv) * d
    return out


def _resolve_xy_overlaps(
    xy: np.ndarray,
    lig_center: np.ndarray,
    min_dist: float,
    push_px: float = 14.0,
    max_iter: int = 36,
) -> np.ndarray:
    """Push nodes apart radially when labels overlap."""
    out = np.array(xy, dtype=float, copy=True)
    L = np.asarray(lig_center, dtype=float).reshape(2)
    n = out.shape[0]
    for _ in range(max_iter):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                d = float(np.linalg.norm(out[i] - out[j]))
                if d >= min_dist or d < 1e-9:
                    continue
                for k in (i, j):
                    v = out[k] - L
                    nv = float(np.linalg.norm(v)) or 1.0
                    out[k] = out[k] + (v / nv) * push_px
                moved = True
        if not moved:
            break
    return out


def compute_residue_positions_on_pocket_offset(
    ligand_atom_xy: list[dict],
    lig_center: np.ndarray,
    n_residues: int,
    margin_px: float = 14.0,
    atom_radius_px: float = 9.0,
    residue_outward_px: float = 54.0,
    curve_eval_points: int = 520,
    closed_spline_s: float = 0.1,
) -> np.ndarray | None:
    """
    Place residue centers on a loop parallel to the ligand pocket (union-of-disks hull),
    just outside the coloured pocket band — fills the gap vs a circular residue ring.
    Returns (n, 2) float array or None if geometry cannot be built.
    """
    if n_residues < 1 or not ligand_atom_xy:
        return None
    hp = _pocket_disk_sample_hull_points(
        ligand_atom_xy,
        margin_px,
        atom_radius_px,
        samples_per_atom=26,
    )
    if hp.shape[0] < 3:
        return None
    curve = _closed_curve_from_hull_polyline(
        hp, n_eval=max(curve_eval_points, n_residues * 24), s_scale=closed_spline_s
    )
    offset_loop = _offset_closed_curve_radially(curve, lig_center, residue_outward_px)
    xy = np.zeros((n_residues, 2), dtype=float)
    for i in range(n_residues):
        t = (i + 0.5) / float(n_residues)
        xy[i] = _point_at_arc_fraction_closed(offset_loop, t)
    return xy


def _smooth_backbone_subpath_from_points(pts_array: np.ndarray) -> str:
    """Single open stroke: smooth spline through ordered points (one contiguous sequence run)."""
    if len(pts_array) < 2:
        return ""
    if len(pts_array) == 2:
        return (
            f"M {pts_array[0][0]:.1f} {pts_array[0][1]:.1f} "
            f"L {pts_array[1][0]:.1f} {pts_array[1][1]:.1f}"
        )
    # splprep requires m > k (point count > spline degree)
    n = len(pts_array)
    k = min(3, max(1, n - 1))
    smoothing = max(n * 1.15, 12.0)
    tck, _ = splprep(
        [pts_array[:, 0], pts_array[:, 1]], s=float(smoothing), per=False, k=k
    )
    uu = np.linspace(0, 1, max(220, len(pts_array) * 18))
    xs, ys = splev(uu, tck)
    path_d = f"M {xs[0]:.1f} {ys[0]:.1f} "
    path_d += " ".join([f"L {x:.1f} {y:.1f}" for x, y in zip(xs[1:], ys[1:])])
    return path_d


def _backbone_quad_control_xy(
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    lig_center: np.ndarray,
    curve_offset: float = 22.0,
) -> tuple[float, float]:
    """Same bulge heuristic as the 2D client: midpoint + perpendicular offset from ligand center."""
    lcx, lcy = float(lig_center[0]), float(lig_center[1])
    mx = (x0 + x1) * 0.5
    my = (y0 + y1) * 0.5
    dx = x1 - x0
    dy = y1 - y0
    length = float(np.hypot(dx, dy)) or 1.0
    perp_x = -dy / length
    perp_y = dx / length
    from_lig_x = mx - lcx
    from_lig_y = my - lcy
    outward = 1.0 if (perp_x * from_lig_x + perp_y * from_lig_y) >= 0 else -1.0
    cx = mx + perp_x * float(curve_offset) * outward
    cy = my + perp_y * float(curve_offset) * outward
    return cx, cy


def _quad_bezier_point(t: float, p0: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> tuple[float, float]:
    omt = 1.0 - float(t)
    t = float(t)
    x = omt * omt * p0[0] + 2.0 * omt * t * p1[0] + t * t * p2[0]
    y = omt * omt * p0[1] + 2.0 * omt * t * p1[1] + t * t * p2[1]
    return float(x), float(y)


def generate_backbone_path(residues, lig_center):
    """
    Smooth backbone for sequence runs (N→N+1).

    If two displayed residues differ by **exactly 2** in sequence number (one missing index,
    e.g. …291, 293…), a short quadratic bridge plus **one** gap dot is drawn. If the skip is
    **larger than 2**, the backbone is **not** bridged (line breaks).
    """
    _path, _dots = generate_backbone_path_and_gap_dots(residues, lig_center)
    return _path


def generate_backbone_path_and_gap_dots(
    residues: list[dict],
    lig_center: np.ndarray,
) -> tuple[str, list[dict[str, float]]]:
    """
    Returns (svg_path, gap_dot_positions).

    Gap dots appear only when ``next_resid - prev_resid == 2`` (single missing residue number).
    """
    if not residues or len(residues) < 2:
        return "", []

    L = np.asarray(lig_center, dtype=float).reshape(2)

    residues_by_chain: dict = {}
    for r in residues:
        chain = r.get("chain", "A")
        if not chain or str(chain).strip() == "" or chain == ".000":
            chain = "A"
        residues_by_chain.setdefault(chain, []).append(r)

    subpaths: list[str] = []
    gap_dots: list[dict[str, float]] = []

    for chain in sorted(residues_by_chain.keys()):
        chain_residues = sorted(
            residues_by_chain[chain],
            key=lambda r: int(r.get("resid", 0)),
        )
        if len(chain_residues) < 2:
            continue

        run_start = 0
        for i in range(1, len(chain_residues)):
            prev_id = int(chain_residues[i - 1].get("resid", 0))
            curr_id = int(chain_residues[i].get("resid", 0))
            if curr_id - prev_id == 1:
                continue

            segment = chain_residues[run_start:i]
            if len(segment) >= 2:
                pts = [
                    (float(r["x"]), float(r["y"]))
                    for r in segment
                    if "x" in r and "y" in r
                ]
                if len(pts) >= 2:
                    subpaths.append(_smooth_backbone_subpath_from_points(np.array(pts, dtype=float)))

            # Single-number skip only: N and N+2 present (difference == 2) → bridge + one dot.
            # Larger skips break the backbone (no connector).
            step = int(curr_id - prev_id)
            if step == 2:
                ra = chain_residues[i - 1]
                rb = chain_residues[i]
                if "x" in ra and "y" in ra and "x" in rb and "y" in rb:
                    x0, y0 = float(ra["x"]), float(ra["y"])
                    x1, y1 = float(rb["x"]), float(rb["y"])
                    cx, cy = _backbone_quad_control_xy(x0, y0, x1, y1, L, curve_offset=22.0)
                    subpaths.append(
                        f"M {x0:.1f} {y0:.1f} Q {cx:.1f} {cy:.1f} {x1:.1f} {y1:.1f}"
                    )
                    p0 = np.array([x0, y0], dtype=float)
                    p1 = np.array([cx, cy], dtype=float)
                    p2 = np.array([x1, y1], dtype=float)
                    gx, gy = _quad_bezier_point(0.5, p0, p1, p2)
                    gap_dots.append({"x": gx, "y": gy})

            run_start = i

        segment = chain_residues[run_start:]
        if len(segment) >= 2:
            pts = [
                (float(r["x"]), float(r["y"]))
                for r in segment
                if "x" in r and "y" in r
            ]
            if len(pts) >= 2:
                subpaths.append(_smooth_backbone_subpath_from_points(np.array(pts, dtype=float)))

    return " ".join(p for p in subpaths if p), gap_dots


def pocket_outline_path(res_nodes, center_xy, base_r=190, thickness=55):
    if not res_nodes:
        return ""
    cx, cy = center_xy

    dists = np.array([n["dist"] for n in res_nodes], dtype=float)
    dmin, dmax = float(dists.min()), float(dists.max())
    denom = max(dmax - dmin, 1e-6)

    pts = []
    for n in res_nodes:
        a = n["angle"]
        r = base_r + thickness * ((n["dist"] - dmin) / denom)
        x = cx + r * np.cos(a)
        y = cy + r * np.sin(a)
        pts.append((x, y))

    pts = np.array(pts, dtype=float)
    pts2 = np.vstack([pts, pts[0]])

    tck, _ = splprep([pts2[:, 0], pts2[:, 1]], s=500, per=True)
    uu = np.linspace(0, 1, 220)
    xs, ys = splev(uu, tck)

    d = f"M {xs[0]:.1f} {ys[0]:.1f} " + " ".join(
        [f"L {x:.1f} {y:.1f}" for x, y in zip(xs[1:], ys[1:])]
    ) + " Z"
    return d


def _pocket_disk_union_boundary_points(
    ligand_atom_xy: list[dict],
    margin_px: float,
    atom_radius_px: float,
    buffer_quad_segs: int = 28,
) -> np.ndarray | None:
    """
    True boundary of ∪ disks around atoms (preserves concavities / “bays” between groups).

    This matches Maestro-style shrink-wrap around the ligand silhouette; unlike ConvexHull,
    indentations between functional groups are kept.
    """
    try:
        from shapely.geometry import Point
        from shapely.ops import unary_union
    except ImportError:
        return None

    if not ligand_atom_xy:
        return np.zeros((0, 2), dtype=float)

    R = float(margin_px + atom_radius_px)
    polys = [
        Point(float(p["x"]), float(p["y"])).buffer(R, quad_segs=max(8, buffer_quad_segs))
        for p in ligand_atom_xy
    ]
    try:
        u = unary_union(polys)
    except Exception:
        return None

    if u.is_empty:
        return None

    if u.geom_type == "Polygon":
        geom = u
    elif u.geom_type == "MultiPolygon":
        geom = max(u.geoms, key=lambda g: float(g.area))
    else:
        return None

    coords = np.asarray(geom.exterior.coords, dtype=float)
    if coords.shape[0] > 1 and np.allclose(coords[0], coords[-1]):
        coords = coords[:-1]
    return coords if coords.shape[0] >= 3 else None


def _pocket_disk_sample_hull_points(
    ligand_atom_xy: list[dict],
    margin_px: float,
    atom_radius_px: float,
    samples_per_atom: int,
) -> np.ndarray:
    """
    Outer envelope around ligand atoms in 2D.

    Prefer the **union-of-disks** boundary (concave, ligand-following). If Shapely is
    unavailable, fall back to convex hull of dense circle samples (smoother blob).
    """
    union_pts = _pocket_disk_union_boundary_points(
        ligand_atom_xy, margin_px, atom_radius_px
    )
    if union_pts is not None and union_pts.shape[0] >= 3:
        return union_pts

    atoms = np.array([[float(p["x"]), float(p["y"])] for p in ligand_atom_xy], dtype=float)
    n_atom = atoms.shape[0]
    if n_atom == 0:
        return np.zeros((0, 2), dtype=float)
    R = float(margin_px + atom_radius_px)
    if n_atom == 1:
        t = np.linspace(0.0, 2.0 * np.pi, max(samples_per_atom, 12), endpoint=False)
        return np.stack(
            [atoms[0, 0] + R * np.cos(t), atoms[0, 1] + R * np.sin(t)],
            axis=1,
        )
    k = max(int(samples_per_atom), 12)
    thetas = np.linspace(0.0, 2.0 * np.pi, k, endpoint=False)
    blocks: list[np.ndarray] = []
    for i in range(n_atom):
        cx, cy = float(atoms[i, 0]), float(atoms[i, 1])
        blocks.append(
            np.stack([cx + R * np.cos(thetas), cy + R * np.sin(thetas)], axis=1)
        )
    all_pts = np.vstack(blocks)
    try:
        hull = ConvexHull(all_pts)
    except Exception:
        return all_pts
    return np.asarray(all_pts[np.asarray(hull.vertices, dtype=int)], dtype=float)


def _closed_curve_from_hull_polyline(hp: np.ndarray, n_eval: int, s_scale: float) -> np.ndarray:
    """Periodic B-spline through ordered hull vertices → dense smooth closed polyline."""
    if hp.shape[0] < 3:
        return hp
    # small periodic smoothing: follows wavy outline without collapsing to a blob
    s = max(float(hp.shape[0]) * float(s_scale), 0.5)
    try:
        tck, _ = splprep([hp[:, 0], hp[:, 1]], s=s, per=True)
    except Exception:
        return hp
    uu = np.linspace(0.0, 1.0, int(n_eval), endpoint=False)
    xs, ys = splev(uu, tck)
    return np.stack([xs, ys], axis=1)


def _buried_mask_for_pocket_points(
    curve_pts: np.ndarray,
    L: np.ndarray,
    res_pts: np.ndarray,
    sector_half_deg: float,
    max_hull_res_dist_px: float,
) -> np.ndarray:
    n = curve_pts.shape[0]
    buried = np.zeros(n, dtype=bool)
    half_w = np.deg2rad(sector_half_deg)
    for i in range(n):
        P = curve_pts[i]
        th = float(np.arctan2(P[1] - L[1], P[0] - L[0]))
        ok = False
        for j in range(res_pts.shape[0]):
            d = res_pts[j] - L
            if np.linalg.norm(d) < 1e-6:
                continue
            ang = float(np.arctan2(d[1], d[0]))
            da = (ang - th + np.pi) % (2.0 * np.pi) - np.pi
            if abs(da) < half_w:
                if np.linalg.norm(res_pts[j] - P) <= max_hull_res_dist_px:
                    ok = True
                    break
        buried[i] = ok
    return buried


def _circular_runs(mask: np.ndarray, link_wrapped_ends: bool = True) -> list[list[int]]:
    n = int(mask.size)
    idx_b = np.where(mask)[0]
    if len(idx_b) == 0:
        return []
    runs: list[list[int]] = []
    cur = [int(idx_b[0])]
    for k in range(1, len(idx_b)):
        prev_i = int(idx_b[k - 1])
        this_i = int(idx_b[k])
        if this_i == prev_i + 1 or (prev_i == n - 1 and this_i == 0):
            cur.append(this_i)
        else:
            runs.append(cur)
            cur = [this_i]
    runs.append(cur)
    if (
        link_wrapped_ends
        and len(runs) >= 2
    ):
        first, last = runs[0], runs[-1]
        if first[0] == 0 and last[-1] == n - 1:
            runs[-1] = last + first
            runs.pop(0)
    return runs


def _hull_segment_to_svg_path(pts: np.ndarray) -> str:
    """Smooth open path through ordered 2D points (tight fit + dense eval = curved, not faceted)."""
    if pts.shape[0] < 2:
        return ""
    if pts.shape[0] == 2:
        return (
            f"M {pts[0, 0]:.1f} {pts[0, 1]:.1f} "
            f"L {pts[1, 0]:.1f} {pts[1, 1]:.1f}"
        )
    if pts.shape[0] == 3:
        return (
            f"M {pts[0, 0]:.1f} {pts[0, 1]:.1f} "
            f"L {pts[1, 0]:.1f} {pts[1, 1]:.1f} "
            f"L {pts[2, 0]:.1f} {pts[2, 1]:.1f}"
        )
    # Low s: follow the wavy outline; many samples: visually smooth, not straight chords
    smoothing = max(float(pts.shape[0]) * 0.06, 0.8)
    try:
        tck, _ = splprep([pts[:, 0], pts[:, 1]], s=smoothing, per=False)
        uu = np.linspace(0, 1, max(48, min(220, pts.shape[0] * 6)))
        xs, ys = splev(uu, tck)
    except Exception:
        parts = [f"M {pts[0, 0]:.1f} {pts[0, 1]:.1f}"]
        for i in range(1, pts.shape[0]):
            parts.append(f"L {pts[i, 0]:.1f} {pts[i, 1]:.1f}")
        return " ".join(parts)

    out = f"M {xs[0]:.1f} {ys[0]:.1f} "
    out += " ".join([f"L {x:.1f} {y:.1f}" for x, y in zip(xs[1:], ys[1:])])
    return out


def _closed_hull_curve_to_svg_path(curve: np.ndarray) -> str:
    """Closed SVG path: dense smooth loop hugging the ligand (Maestro-style pocket band)."""
    if curve.shape[0] < 2:
        return ""
    xs = curve[:, 0]
    ys = curve[:, 1]
    parts = [f"M {xs[0]:.1f} {ys[0]:.1f}"]
    for i in range(1, xs.shape[0]):
        parts.append(f"L {xs[i]:.1f} {ys[i]:.1f}")
    parts.append("Z")
    return " ".join(parts)


def pocket_outline_paths_ligand_hug(
    ligand_atom_xy: list[dict],
    lig_center: np.ndarray,
    residue_xy: list[tuple[float, float]],
    margin_px: float = 12.0,
    atom_radius_px: float = 9.0,
    samples_per_atom: int = 26,
    curve_eval_points: int = 320,
    closed_spline_s: float = 0.12,
    sector_half_deg: float = 28.0,
    max_hull_res_dist_px: float = 280.0,
    continuous_band: bool = True,
) -> list[str]:
    """
    Pocket boundary as one or more SVG path strings (Maestro-style: with continuous_band=False,
    segments skip solvent-facing arcs so the dashed outline shows a visible pocket entrance gap).
    """
    if not ligand_atom_xy:
        return []
    L = np.array([float(lig_center[0]), float(lig_center[1])], dtype=float)
    res_pts = np.array(residue_xy, dtype=float) if residue_xy else np.zeros((0, 2), dtype=float)

    hp = _pocket_disk_sample_hull_points(ligand_atom_xy, margin_px, atom_radius_px, samples_per_atom)
    if hp.shape[0] < 3:
        return []

    curve = _closed_curve_from_hull_polyline(hp, n_eval=curve_eval_points, s_scale=closed_spline_s)

    if continuous_band:
        return [_closed_hull_curve_to_svg_path(curve)]

    if res_pts.shape[0] == 0:
        buried = np.ones(curve.shape[0], dtype=bool)
    else:
        buried = _buried_mask_for_pocket_points(curve, L, res_pts, sector_half_deg, max_hull_res_dist_px)

    if not np.any(buried):
        return [_closed_hull_curve_to_svg_path(curve)]

    runs = _circular_runs(buried)
    path_chunks: list[str] = []
    for run in runs:
        if len(run) < 2:
            continue
        pts = curve[np.array(run, dtype=int)]
        chunk = _hull_segment_to_svg_path(pts)
        if chunk:
            path_chunks.append(chunk)

    if path_chunks:
        return path_chunks
    return [_closed_hull_curve_to_svg_path(curve)]


def pocket_outline_paths_maestro_entrance_gap(
    ligand_atom_xy: list[dict],
    lig_center: np.ndarray,
    residue_xy: list[tuple[float, float]],
    margin_px: float = 14.0,
    atom_radius_px: float = 9.0,
    samples_per_atom: int = 26,
    curve_eval_points: int = 320,
    closed_spline_s: float = 0.12,
    gap_half_deg: float = 26.0,
) -> list[str]:
    """
    Full ligand-hugging pocket contour split into one or two arcs with a **gap** centered on
    the direction opposite the pocket residue centroid (Maestro / Nature LID: pocket entrance).
    """
    if not ligand_atom_xy:
        return []
    L = np.array([float(lig_center[0]), float(lig_center[1])], dtype=float)
    res_pts = np.array(residue_xy, dtype=float) if residue_xy else np.zeros((0, 2), dtype=float)

    hp = _pocket_disk_sample_hull_points(ligand_atom_xy, margin_px, atom_radius_px, samples_per_atom)
    if hp.shape[0] < 3:
        return []

    curve = _closed_curve_from_hull_polyline(hp, n_eval=curve_eval_points, s_scale=closed_spline_s)

    if res_pts.shape[0] == 0:
        return [_closed_hull_curve_to_svg_path(curve)]

    centroid = res_pts.mean(axis=0)
    v = centroid - L
    norm = float(np.linalg.norm(v))
    if norm < 1e-9:
        return [_closed_hull_curve_to_svg_path(curve)]

    # Opening toward bulk solvent: opposite the in-pocket residue cluster.
    # Carve indices (not just angles) so the gap always appears in the polyline samples.
    gap_center = float(np.arctan2(v[1], v[0]) + np.pi)
    ang = np.arctan2(curve[:, 1] - L[1], curve[:, 0] - L[0])
    da = (ang - gap_center + np.pi) % (2.0 * np.pi) - np.pi
    i_c = int(np.argmin(np.abs(da)))
    n = int(curve.shape[0])
    # Total removed arc ~2*gap_half_deg around the pocket entrance
    width = max(12, int(n * (2.0 * float(gap_half_deg) / 360.0)))
    width = min(width, max(8, n // 2 - 2))
    keep = np.ones(n, dtype=bool)
    for k in range(-width, width + 1):
        keep[(i_c + k) % n] = False

    if np.all(keep) or not np.any(keep):
        return [_closed_hull_curve_to_svg_path(curve)]

    # Do not merge head/tail: we need two separate arcs with a real gap between them.
    runs = _circular_runs(keep, link_wrapped_ends=False)
    path_chunks: list[str] = []
    for run in runs:
        if len(run) < 2:
            continue
        pts = curve[np.array(run, dtype=int)]
        chunk = _hull_segment_to_svg_path(pts)
        if chunk:
            path_chunks.append(chunk)

    return path_chunks if path_chunks else [_closed_hull_curve_to_svg_path(curve)]


def pocket_outline_ligand_hug(
    ligand_atom_xy: list[dict],
    lig_center: np.ndarray,
    residue_xy: list[tuple[float, float]],
    margin_px: float = 12.0,
    atom_radius_px: float = 9.0,
    samples_per_atom: int = 26,
    curve_eval_points: int = 320,
    closed_spline_s: float = 0.12,
    sector_half_deg: float = 28.0,
    max_hull_res_dist_px: float = 280.0,
    continuous_band: bool = True,
) -> str:
    """
    Pocket outline: smooth closed band following the ligand silhouette (union-of-disks hull + spline).

    When continuous_band is True (default), the full contour is always drawn—irregular shape
    around rings/tails like Maestro, and the pocket is never dropped when residue positions
    do not match the old “buried sector” mask.

    When continuous_band is False, only arc segments with a nearby residue (solvent gaps) are
    drawn, matching the previous behaviour.
    """
    parts = pocket_outline_paths_ligand_hug(
        ligand_atom_xy,
        lig_center,
        residue_xy,
        margin_px=margin_px,
        atom_radius_px=atom_radius_px,
        samples_per_atom=samples_per_atom,
        curve_eval_points=curve_eval_points,
        closed_spline_s=closed_spline_s,
        sector_half_deg=sector_half_deg,
        max_hull_res_dist_px=max_hull_res_dist_px,
        continuous_band=continuous_band,
    )
    return " ".join(parts) if parts else ""


# =========================
# Atom-anchored interaction detection
# =========================

def is_backbone_atom(atom_name: str) -> bool:
    return atom_name.strip().upper() in {"N", "CA", "C", "O", "OXT"}


_POLAR_ELEMS = {"N", "O", "S"}
_HBOND_ELEMS = {"N", "O", "S", "F"}
_HYDROPHOBIC_ELEMS = {"C", "S", "F", "CL", "BR", "I"}

# Rough protein-side H-bond atom-role inference by atom name
_PROT_HBOND_DONOR_PREFIXES = (
    "N", "NE", "NH", "NZ", "ND", "OG", "OH", "SG"
)
_PROT_HBOND_ACCEPTOR_PREFIXES = (
    "O", "OD", "OE", "OG", "OH", "SD", "ND1", "NE2"
)


def _norm_elem(atom) -> str:
    elem = (getattr(atom, "element", None) or "").strip().upper()
    if elem:
        return elem
    # PDBs are often messy; fall back to atom-name first letter.
    return (getattr(atom, "name", "")[:1] or "").strip().upper()


def _is_protein_hbond_donor(atom_name: str) -> bool:
    n = atom_name.strip().upper()
    return any(n.startswith(p) for p in _PROT_HBOND_DONOR_PREFIXES)


def _is_protein_hbond_acceptor(atom_name: str) -> bool:
    n = atom_name.strip().upper()
    return any(n.startswith(p) for p in _PROT_HBOND_ACCEPTOR_PREFIXES)


def _best_interaction_type_for_pair(
    *,
    min_dist: float,
    rname: str,
    prot_atom_name: str,
    prot_elem: str,
    lig_elem: str,
    lig_rd_atom,
    prolif_hits: set[str],
) -> str:
    """
    Prioritized interaction assignment for the closest ligand-protein atom pair.
    Priority: salt_bridge > hbond > hydrophobic > contact.
    """
    prolif_s = " ".join(sorted(prolif_hits)).lower()
    lig_charge = int(lig_rd_atom.GetFormalCharge()) if lig_rd_atom is not None else 0

    residue_positive = rname in POS_RES
    residue_negative = rname in NEG_RES
    residue_hydrophobic = rname in HYDROPHOBIC_RES

    # --- Salt bridge (charged residue + opposite-charge ligand atom) ---
    # Keep a permissive fallback if formal charge is unavailable.
    has_prolif_salt = any(k in prolif_s for k in ("ionic", "salt", "cation", "anion"))
    opposite_charge_pair = (
        (residue_positive and lig_charge < 0) or
        (residue_negative and lig_charge > 0)
    )
    charged_fallback = (
        (residue_positive or residue_negative) and lig_elem in _POLAR_ELEMS
    )
    if min_dist <= 4.2 and (opposite_charge_pair or (has_prolif_salt and charged_fallback)):
        return "salt_bridge"

    # --- Hydrogen bond ---
    has_prolif_hbond = any(k in prolif_s for k in ("hbond", "hba", "hbd"))
    prot_can_donate = _is_protein_hbond_donor(prot_atom_name)
    prot_can_accept = _is_protein_hbond_acceptor(prot_atom_name)
    lig_can_hbond = lig_elem in _HBOND_ELEMS
    prot_can_hbond = prot_elem in _HBOND_ELEMS and (prot_can_donate or prot_can_accept)
    if min_dist <= 3.6 and lig_can_hbond and prot_can_hbond:
        return "hbond"
    if has_prolif_hbond and min_dist <= 3.9 and lig_can_hbond and prot_elem in _HBOND_ELEMS:
        return "hbond"

    # --- Hydrophobic contact ---
    has_prolif_hydrophobic = "hydroph" in prolif_s
    if min_dist <= 4.8 and residue_hydrophobic and lig_elem in _HYDROPHOBIC_ELEMS:
        return "hydrophobic"
    if has_prolif_hydrophobic and min_dist <= 5.0 and lig_elem in _HYDROPHOBIC_ELEMS:
        return "hydrophobic"

    return "contact"


def _dedupe_interactions_by_residue_closest(interactions: list[dict]) -> list[dict]:
    """
    Keep only the closest interaction per residue (chain+resname+resid).
    """
    best: dict[tuple[str, str, int], dict] = {}
    for it in interactions:
        key = (
            str(it.get("chain", "A")),
            str(it.get("resname", "")),
            int(it.get("resid", -1)),
        )
        prev = best.get(key)
        if prev is None or float(it.get("distance", 1e9)) < float(prev.get("distance", 1e9)):
            best[key] = it
    return list(best.values())


def detect_interactions_atom_anchored(
    u: mda.Universe,
    mol: Chem.Mol,
    ligand_resname: str,
    cutoff_contact: float = 5.0,
):
    lig = u.select_atoms(f"resname {ligand_resname}")
    if lig.n_atoms == 0:
        raise ValueError(f"No ligand found with resname '{ligand_resname}' in PDB.")

    binding = u.select_atoms("protein or nucleic")

    map_pdb_to_rd = map_pdb_ligand_atoms_to_rdkit(u, mol, ligand_resname)

    interactions = []

    # ProLIF base interactions (if available), used as additional evidence.
    prolif_by_residue: dict[str, set[str]] = {}
    if plf is not None:
        try:
            lig_plf = plf.Molecule.from_mda(lig)
            prot_plf = plf.Molecule.from_mda(binding)
            fp = plf.Fingerprint()
            ifp = fp.generate(lig_plf, prot_plf)
            for key in ifp.keys():
                # key format generally contains residue + interaction class; keep robust.
                key_s = str(key)
                for res in binding.residues:
                    rid = f"{res.resname}{int(res.resid)}"
                    if rid in key_s:
                        prolif_by_residue.setdefault(rid, set()).add(key_s)
        except Exception:
            # Keep pipeline robust if ProLIF parsing fails for some structures.
            prolif_by_residue = {}

    for res in binding.residues:
        res_atoms = res.atoms
        if res_atoms.n_atoms == 0:
            continue

        dmat = mda.lib.distances.distance_array(lig.positions, res_atoms.positions)
        min_idx = np.unravel_index(np.argmin(dmat), dmat.shape)
        lig_i = int(min_idx[0])
        res_j = int(min_idx[1])
        min_dist = float(dmat[min_idx])

        if min_dist > cutoff_contact:
            continue

        lig_atom = lig.atoms[lig_i]
        prot_atom = res_atoms[res_j]
        prot_atom_name = prot_atom.name.strip()

        rd_idx = map_pdb_to_rd.get(lig_i, -1)
        lig_rd_atom = mol.GetAtomWithIdx(int(rd_idx)) if int(rd_idx) >= 0 else None

        rname = res.resname
        lig_elem = _norm_elem(lig_atom)
        prot_elem = _norm_elem(prot_atom)
        residue_key = f"{rname}{int(res.resid)}"
        prolif_hits = prolif_by_residue.get(residue_key, set())
        itype = _best_interaction_type_for_pair(
            min_dist=min_dist,
            rname=rname,
            prot_atom_name=prot_atom_name,
            prot_elem=prot_elem,
            lig_elem=lig_elem,
            lig_rd_atom=lig_rd_atom,
            prolif_hits=prolif_hits,
        )

        interactions.append({
            "residue": residue_key,
            "resname": rname,
            "resid": int(res.resid),
            "chain": getattr(res, "segid", "") or "A",
            "res_class": classify_residue(rname),
            "type": itype,
            "distance": round(min_dist, 2),
            "ligand_atom_index": int(rd_idx),
            "protein_atom_name": prot_atom_name,
            "backbone": bool(is_backbone_atom(prot_atom_name)),
            "prolif_match": bool(prolif_hits),
        })

    # Defensive dedupe: one interaction per residue, closest wins.
    interactions = _dedupe_interactions_by_residue_closest(interactions)

    # Stable ordering: strongest types first, then shortest distances.
    priority = {"salt_bridge": 0, "hbond": 1, "hydrophobic": 2, "contact": 3}
    interactions.sort(key=lambda it: (priority.get(str(it.get("type")), 99), float(it.get("distance", 1e9))))
    return interactions


# =========================
# Main builder
# =========================

def build_diagram(
    pdb_bytes: bytes,
    ligand_resname: str | None,
    pocket_radius: float,
    svg_w: int,
    svg_h: int,
    sdf_bytes: bytes | None = None,
    model_index: int = 1,
):
    selected_pdb_text, total_models = extract_selected_model_from_text(
        pdb_bytes.decode("utf-8"),
        model_index=model_index,
    )
    selected_pdb_bytes = selected_pdb_text.encode("utf-8")

    # 1) Get ligand molecule: from SDF if provided, else from PDB (ATOM + CONECT)
    if sdf_bytes:
        mol = Chem.MolFromMolBlock(sdf_bytes.decode("utf-8"), removeHs=False)
        if mol is None:
            raise ValueError("Failed to parse ligand SDF (MolBlock).")
        if not ligand_resname:
            ligand_resname = detect_ligand_resname_from_pdb(selected_pdb_bytes)
    else:
        resname = ligand_resname
        if not resname:
            resname = detect_ligand_resname_from_pdb(selected_pdb_bytes)
        mol = ligand_mol_from_pdb(selected_pdb_bytes, resname)
        ligand_resname = resname

    # 2D depiction and atom anchors: heavy atoms only (no explicit H in SVG)
    mol_draw = Chem.RemoveHs(Chem.Mol(mol))
    if mol_draw.GetNumAtoms() == 0:
        mol_draw = Chem.Mol(mol)

    try:
        AllChem.Compute2DCoords(mol_draw)
    except (ValueError, RuntimeError):
        # Fallback for symmetric/fused ring systems: use 3D -> 2D projection (x, y)
        if mol_draw.GetNumConformers():
            conf = mol_draw.GetConformer()
            for i in range(mol_draw.GetNumAtoms()):
                p = conf.GetAtomPosition(i)
                conf.SetAtomPosition(i, (p.x, p.y, 0.0))

    lig_xy = ligand_2d_coords(mol_draw)

    # Make ligand smaller by increasing padding (key change)
    # Smaller pad → larger ligand in the diagram (Maestro-style prominence).
    lig_xy_px = normalize_to_svg_coords(lig_xy, svg_w, svg_h, pad=88)

    # Provide per-atom pixel coords for atom-anchored lines (indices = heavy-atom order)
    ligand_atom_xy = [{"x": float(x), "y": float(y)} for x, y in lig_xy_px]

    lig_center = lig_xy_px.mean(axis=0)
    svg = rdkit_svg(mol_draw, w=svg_w, h=svg_h)

    heavy_rd_index = rdkit_full_to_heavy_index_map(mol)

    # 2) Load complex
    u = universe_from_pdb_bytes(selected_pdb_bytes)
    ligand_resname = _resolve_single_ligand_resname(u, ligand_resname)

    # 3) All residues within pocket_radius of the ligand (layout + backbone use this full set)
    res_nodes = place_residues_ring(u, ligand_resname, pocket_radius)

    # 4) Interactions ranked/filtered for drawing lines (does not remove pocket residues)
    interactions = detect_interactions_atom_anchored(u, mol, ligand_resname, cutoff_contact=5.0)
    total_interactions = len(interactions)
    interactions = rank_and_filter_interactions(interactions)
    logger.debug(
        "interactions: total=%d filtered=%d | pocket residues (radius=%s Å): %d",
        total_interactions,
        len(interactions),
        pocket_radius,
        len(res_nodes),
    )

    # Remap ligand_atom_index to heavy-only indices for the frontend
    for it in interactions:
        rd = it.get("ligand_atom_index")
        if rd is not None and int(rd) >= 0:
            it["ligand_atom_index"] = heavy_rd_index.get(int(rd), -1)

    # 5) Place residues on a loop parallel to the pocket hull (not a uniform circle),
    #    so the chain follows the binding-site contour and sits just outside the pocket band.
    max_dist = max((n["dist"] for n in res_nodes), default=1.0)
    max_dist = max(max_dist, 1e-6)

    R_base = 238.0
    R_amp = 118.0
    min_node_dist = 76.0

    n_res = len(res_nodes)
    residues = []
    placed_positions: list[tuple[float, float]] = []
    rng = np.random.default_rng(7)

    pocket_xy = None
    if n_res >= 1:
        pocket_xy = compute_residue_positions_on_pocket_offset(
            ligand_atom_xy,
            lig_center,
            n_res,
            margin_px=14.0,
            atom_radius_px=9.0,
            residue_outward_px=52.0,
            curve_eval_points=520,
        )

    if pocket_xy is not None and pocket_xy.shape[0] == n_res:
        pocket_xy = _resolve_xy_overlaps(
            pocket_xy, np.asarray(lig_center, dtype=float), min_node_dist
        )
        for i, node in enumerate(res_nodes):
            x = float(pocket_xy[i, 0])
            y = float(pocket_xy[i, 1])
            placed_positions.append((x, y))
            residues.append({**node, "x": x, "y": y})
    else:
        for n in res_nodes:
            dynamic_R = R_base + R_amp * (n["dist"] / max_dist)
            angle = n["angle"] + float(rng.uniform(-0.03, 0.03))

            x = float(lig_center[0] + dynamic_R * np.cos(angle))
            y = float(lig_center[1] + dynamic_R * np.sin(angle))

            for _ in range(22):
                overlap = False
                for (px, py) in placed_positions:
                    if float(np.hypot(x - px, y - py)) < min_node_dist:
                        overlap = True
                        break
                if not overlap:
                    break
                dynamic_R += 34.0
                x = float(lig_center[0] + dynamic_R * np.cos(angle))
                y = float(lig_center[1] + dynamic_R * np.sin(angle))

            placed_positions.append((x, y))
            residues.append({**n, "x": x, "y": y, "_angle": angle, "_R": dynamic_R})

        lig_cx, lig_cy = float(lig_center[0]), float(lig_center[1])
        for _ in range(16):
            moved = False
            for i in range(len(residues)):
                r_i = residues[i]
                xi, yi = r_i["x"], r_i["y"]
                angle_i = r_i["_angle"]
                R_i = r_i["_R"]
                for j in range(len(residues)):
                    if i == j:
                        continue
                    r_j = residues[j]
                    d = float(np.hypot(xi - r_j["x"], yi - r_j["y"]))
                    if d < min_node_dist and d > 1e-6:
                        R_new = R_i + 36.0
                        residues[i]["x"] = float(lig_cx + R_new * np.cos(angle_i))
                        residues[i]["y"] = float(lig_cy + R_new * np.sin(angle_i))
                        residues[i]["_R"] = R_new
                        placed_positions[i] = (residues[i]["x"], residues[i]["y"])
                        moved = True
                        break
            if not moved:
                break

        for r in residues:
            r.pop("_angle", None)
            r.pop("_R", None)

    for r in residues:
        r.setdefault("strain_score", 0.0)

    # 6) Backbone: contiguous runs + quadratic bridges across missing sequence numbers (gap dots)
    backbone_path, backbone_gap_dots = generate_backbone_path_and_gap_dots(
        residues, np.asarray(lig_center, dtype=float).reshape(2)
    )

    # 7) Pocket: dotted boundary with a visible entrance gap (Nature / Maestro LID)
    residue_xy_list = [(float(r["x"]), float(r["y"])) for r in residues]
    pocket_outline_paths = pocket_outline_paths_maestro_entrance_gap(
        ligand_atom_xy,
        lig_center,
        residue_xy_list,
        margin_px=14.0,
        atom_radius_px=9.0,
        samples_per_atom=26,
        curve_eval_points=320,
        gap_half_deg=26.0,
    )
    outline = " ".join(pocket_outline_paths) if pocket_outline_paths else ""

    ligand_atom_solvent_exposed = compute_ligand_atom_solvent_exposed(
        u, ligand_resname, mol, mol_draw, min_protein_dist_angstrom=3.45
    )

    return {
        "svg": svg,
        "ligand_center": [float(lig_center[0]), float(lig_center[1])],
        "ligand_atom_xy": ligand_atom_xy,
        "ligand_atom_solvent_exposed": ligand_atom_solvent_exposed,
        "residues": residues,
        "backbone_path": backbone_path,
        "backbone_gap_dots": backbone_gap_dots,
        "pocket_outline_path": outline,
        "pocket_outline_paths": pocket_outline_paths,
        "interactions": interactions,
        "meta": {
            "pocket_radius": pocket_radius,
            "ligand_resname": ligand_resname,
            "svg_w": svg_w,
            "svg_h": svg_h,
            "selected_model": model_index,
            "total_models": total_models,
            "ligand_atoms": int(u.select_atoms(f"resname {ligand_resname}").n_atoms),
            "protein_atoms": int(u.select_atoms("protein").n_atoms),
            "interaction_total": total_interactions,
            "interaction_filtered": len(interactions),
        },
        "selected_model": model_index,
        "total_models": total_models,
        "ligand_atoms": int(u.select_atoms(f"resname {ligand_resname}").n_atoms),
        "protein_atoms": int(u.select_atoms("protein").n_atoms),
    }