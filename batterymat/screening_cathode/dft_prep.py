"""VASP input generation for sequential supercell delithiation voltage curves."""
import json
import os
import re
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Vacancy import kept for backward compatibility but no longer used internally.
from jarvis.analysis.defects.vacancy import Vacancy  # noqa: F401


# ---------------------------------------------------------------------------
# Voltage helper
# ---------------------------------------------------------------------------

def compute_dft_voltage(
    e_lithiated: float,
    e_delithiated: float,
    n_li: int,
    e_li_metal: float,
) -> float:
    """Compute average intercalation voltage from PBE total energies.

    Formula: V = (E_delithiated - E_lithiated - n_li * E_li_metal) / n_li
    """
    if n_li == 0:
        raise ValueError("n_li must be > 0")
    return (e_delithiated - e_lithiated - n_li * e_li_metal) / n_li


# ---------------------------------------------------------------------------
# PAW potentials from JARVIS default_potcars.json
# ---------------------------------------------------------------------------
def _load_jarvis_paw() -> Dict[str, str]:
    """Load PAW potential labels from JARVIS default_potcars.json."""
    import jarvis
    potcar_json = Path(jarvis.__file__).parent / "io" / "vasp" / "default_potcars.json"
    return json.loads(potcar_json.read_text())


_JARVIS_PAW = _load_jarvis_paw()

# ---------------------------------------------------------------------------
# Hubbard U values (Dudarev, LDAUTYPE=2)
# From Materials Project: https://docs.materialsproject.org/methodology/
# materials-methodology/calculation-details/gga+u-calculations/hubbard-u-values
# JARVIS only supports uniform U, so we maintain element-specific values here.
# ---------------------------------------------------------------------------
_HUBBARD_U = {
    "Mn": 3.9,
    "Fe": 5.3,
    "Co": 3.32,
    "Ni": 6.2,
    "V": 3.25,
    "Cr": 3.7,
    "Mo": 4.38,
    "W": 6.2,
}

_OPTB88VDW_TAGS = """\
GGA = OR
LUSE_VDW = .TRUE.
AGGAC = 0.0
"""

_TMBJ_INCAR = """\
IBRION = -1
NSW = 0
ENCUT = 520
EDIFF = 1E-6
PREC = Accurate
ISMEAR = 0
SIGMA = 0.05
METAGGA = MBJ
LASPH = .TRUE.
ICHARG = 2
ALGO = All
TIME = 0.4
LORBIT = 11
NEDOS = 2000
NELM = 1000
AMIX = 0.02
BMIX = 0.001
AMIN = 0.01
LWAVE = .TRUE.
LCHARG = .TRUE.
ISPIN = 2
LMAXMIX = 4
ISTART = 0
ISYM = 0
NCORE = 8
"""

_E_LI_METAL = {
    # Reference Li metal energies, computed in-house with the SAME PAW (Li_sv),
    # ENCUT=520, ISMEAR=1, SIGMA=0.1, k=17x17x17 as the cathode runs.
    # See dft_inputs/JVASP-913-Li/Li_sv_PBE/ and dft_inputs/JVASP-913-Li/Li_sv_optB88vdW/ (GGA=OR).
    # The old JARVIS value (-0.925) used a different PAW/basis and was wrong by ~1 eV.
    "pbe": -1.9031,
    "optb88vdw": -0.9646,
}

_LAYERED_SPACEGROUPS = {
    166,  # R-3m (LiCoO2, NMC, NCA)
    194,  # P63/mmc (O3-type, graphite)
    12,   # C2/m (Li2MnO3, Li-rich layered)
    15,   # C2/c (some Li-rich layered)
}

_PBE_INCAR = """\
IBRION = 1
NSW = 100
ISIF = 3
ENCUT = 520
EDIFF = 1E-4
PREC = Accurate
ISMEAR = 0
SIGMA = 0.05
LWAVE = .TRUE.
LCHARG = .TRUE.
NELM = 100
NELMIN = 4
ISPIN = 2
LASPH = .TRUE.
LMAXMIX = 4
ISTART = 1
ICHARG = 1
ISYM = 0
LORBIT = 11
LREAL = Auto
NCORE = 8
KPAR = 2
NSIM = 8
AMIX = 0.2
BMIX = 0.0001
AMIX_MAG = 0.4
BMIX_MAG = 0.0001
"""

_DEFAULT_MAGMOM = {
    "Mn": 4.0,
    "Fe": 5.0,
    "Ni": 2.0,
    "Co": 0.6,
    "V": 3.0,
    "Cr": 3.0,
    "Mo": 3.0,
    "W": 2.0,
}


# ---------------------------------------------------------------------------
# Layered-structure detection
# ---------------------------------------------------------------------------

def _get_spacegroup_number(jid: str, dft3d_df=None) -> Optional[int]:
    """Look up spacegroup number from JARVIS dft_3d DataFrame."""
    if dft3d_df is None:
        return None
    import pandas as pd
    for _, row in dft3d_df.iterrows():
        if row.get("jid") == jid and "spg_number" in row:
            try:
                return int(row["spg_number"])
            except (ValueError, TypeError):
                return None
    return None


def _is_layered(jid: str, dft3d_df=None) -> bool:
    """Check if a material is layered based on spacegroup whitelist."""
    spg = _get_spacegroup_number(jid, dft3d_df)
    if spg is None:
        return False
    return spg in _LAYERED_SPACEGROUPS


def _resolve_functional(jid: str, dft3d_df=None, functional: str = "auto") -> str:
    """Resolve functional choice.

    Args:
        jid:        JARVIS JID.
        dft3d_df:   Pre-loaded JARVIS-DFT DataFrame (needs spg_number column).
        functional: "auto", "pbe", or "optb88vdw".

    Returns:
        Resolved functional string: "pbe" or "optb88vdw".
    """
    if functional == "auto":
        return "optb88vdw" if _is_layered(jid, dft3d_df) else "pbe"
    if functional in ("pbe", "optb88vdw"):
        return functional
    raise ValueError(f"Unknown functional: {functional!r}. Use 'auto', 'pbe', or 'optb88vdw'.")


# ---------------------------------------------------------------------------
# INCAR / KPOINTS / POTCAR helpers
# ---------------------------------------------------------------------------

def _hubbard_u_lines(elements: List[str]) -> str:
    """Generate LDAU INCAR lines for elements that have Hubbard U values."""
    needs_u = any(el in _HUBBARD_U for el in elements)
    if not needs_u:
        return ""
    ldaul = [2 if el in _HUBBARD_U else -1 for el in elements]
    ldauu = [_HUBBARD_U.get(el, 0.0) for el in elements]
    ldauj = [0.0] * len(elements)
    lines = [
        "LDAU = .TRUE.",
        "LDAUTYPE = 2",
        f"LDAUL = {' '.join(str(v) for v in ldaul)}",
        f"LDAUU = {' '.join(f'{v:.2f}' for v in ldauu)}",
        f"LDAUJ = {' '.join(f'{v:.2f}' for v in ldauj)}",
    ]
    return "\n".join(lines) + "\n"


def write_relax_incar(directory, elements: List[str] = None,
                      isif: int = 3, nsw: int = 100,
                      magmom_str: str = None,
                      functional: str = "pbe") -> None:
    """Write relaxation INCAR to directory.

    Args:
        functional: "pbe" or "optb88vdw". When "optb88vdw", appends
                    GGA=OR, LUSE_VDW, AGGAC tags for optB88-vdW.
    """
    content = _PBE_INCAR
    content = content.replace("ISIF = 3", f"ISIF = {isif}")
    content = content.replace("NSW = 100", f"NSW = {nsw}")
    if functional == "optb88vdw":
        content += _OPTB88VDW_TAGS
    if elements:
        content += _hubbard_u_lines(elements)
    if magmom_str:
        content += f"MAGMOM = {magmom_str}\n"
    Path(directory, "INCAR").write_text(content)


def write_pbe_relax_incar(directory, elements: List[str] = None,
                          isif: int = 3, nsw: int = 100,
                          magmom_str: str = None) -> None:
    """Write PBE relaxation INCAR to directory. Alias for write_relax_incar()."""
    write_relax_incar(directory, elements=elements, isif=isif, nsw=nsw,
                      magmom_str=magmom_str, functional="pbe")


def write_tmbj_incar(directory, elements: List[str] = None,
                     magmom_str: str = None) -> None:
    """Write TB-mBJ static INCAR to directory.

    Uses ICHARG=2 (atomic superposition) so MAGMOM is respected.
    No DFT+U — MBJ is a separate functional.
    """
    content = _TMBJ_INCAR
    if magmom_str:
        content += f"MAGMOM = {magmom_str}\n"
    Path(directory, "INCAR").write_text(content)


def write_kpoints(directory, mesh: str = "3 3 3") -> None:
    """Write a Gamma-centered KPOINTS to directory."""
    content = f"Automatic mesh\n0\nGamma\n{mesh}\n0 0 0\n"
    Path(directory, "KPOINTS").write_text(content)


def write_potcar_spec(directory, elements: List[str]) -> None:
    """Write POTCAR_spec listing recommended PAW labels for elements."""
    labels = [_JARVIS_PAW.get(el, el) for el in elements]
    Path(directory, "POTCAR_spec").write_text("\n".join(labels) + "\n")


# ---------------------------------------------------------------------------
# ALIGNN energy helper
# ---------------------------------------------------------------------------

def _alignn_energy(atoms) -> float:
    """Compute ALIGNN ML force field energy for a jarvis Atoms object."""
    try:
        from alignn.ff.ff import AlignnAtomwiseCalculator, wt01_path
    except ImportError as e:
        raise ImportError(
            f"ALIGNN not available: {e}. Verify your ALIGNN installation."
        ) from e
    calc = AlignnAtomwiseCalculator(path=wt01_path(), stress_wt=0.3)
    ase_atoms = atoms.ase_converter()
    ase_atoms.calc = calc
    return ase_atoms.get_potential_energy()


# ---------------------------------------------------------------------------
# Vacancy selection
# ---------------------------------------------------------------------------

def get_next_vacancy(atoms) -> Tuple:
    """Return the lowest-energy single-Li-vacancy structure with metadata.

    Evaluates ALIGNN energy for every Li vacancy (no Wyckoff deduplication)
    to ensure the best site is selected even after symmetry is broken by
    prior vacancies.

    Args:
        atoms: jarvis.core.atoms.Atoms object.

    Returns:
        (defect_structure, removed_atom_index, info_dict)
        info_dict has keys: 'n_candidates', 'candidate_energies', 'selected_energy'
    """
    li_indices = [i for i, el in enumerate(atoms.elements) if el == "Li"]
    n_li = len(li_indices)

    if n_li == 0:
        raise ValueError("Structure has no Li atoms — cannot create vacancy.")

    if n_li == 1:
        result = atoms.remove_site_by_index(li_indices[0])
        return result, li_indices[0], {
            "n_candidates": 1,
            "candidate_energies": [],
            "selected_energy": None,
        }

    # Evaluate every Li vacancy with ALIGNN — no Wyckoff deduplication
    candidates = []
    energies = []
    for idx in li_indices:
        defect = atoms.remove_site_by_index(idx)
        e = _alignn_energy(defect)
        candidates.append((defect, idx))
        energies.append(e)

    best_i = energies.index(min(energies))
    best_defect, best_idx = candidates[best_i]

    return best_defect, best_idx, {
        "n_candidates": len(li_indices),
        "candidate_energies": energies,
        "selected_energy": energies[best_i],
    }


def get_delithiated_structure(atoms):
    """Return the lowest-energy single-Li-vacancy structure.

    Backward-compatible wrapper around get_next_vacancy().
    """
    structure, _, _ = get_next_vacancy(atoms)
    return structure


# ---------------------------------------------------------------------------
# POSCAR writer
# ---------------------------------------------------------------------------

def _write_poscar(directory: Path, atoms) -> None:
    """Write VASP POSCAR from a jarvis Atoms object."""
    from jarvis.io.vasp.inputs import Poscar
    poscar = Poscar(atoms)
    poscar.write_file(str(Path(directory, "POSCAR")))


# ---------------------------------------------------------------------------
# Supercell helpers
# ---------------------------------------------------------------------------

_DEFAULT_OUTPUT_DIR = Path(__file__).parent / "dft_inputs"


def _min_supercell_dim(atoms, min_length: float = 7.0) -> List[int]:
    """Compute smallest supercell dimensions so all lattice vectors >= min_length."""
    import math
    lattice = np.array(atoms.lattice_mat)
    norms = np.linalg.norm(lattice, axis=1)
    return [max(1, math.ceil(min_length / n)) for n in norms]


def _afm_magmom(atoms) -> List[float]:
    """Generate an AFM MAGMOM list for any structure."""
    magmom = []
    sign = 1
    for el in atoms.elements:
        mag = _DEFAULT_MAGMOM.get(el, 0.0)
        if mag != 0.0:
            magmom.append(sign * mag)
            sign *= -1
        else:
            magmom.append(0.0)
    return magmom


def _supercell_magmom(supercell, primitive, primitive_magmom, dim):
    """Map each supercell atom to its primitive-cell equivalent and assign MAGMOM."""
    prim_frac = np.array(primitive.frac_coords)
    sup_frac = np.array(supercell.frac_coords)
    dim_arr = np.array(dim, dtype=float)

    magmom = []
    for i, sf in enumerate(sup_frac):
        pf = (sf * dim_arr) % 1.0
        diffs = pf - prim_frac
        diffs = diffs - np.round(diffs)
        dists = np.linalg.norm(diffs, axis=1)
        nearest = np.argmin(dists)
        if supercell.elements[i] != primitive.elements[nearest]:
            raise ValueError(
                f"Supercell atom {i} ({supercell.elements[i]}) mapped to "
                f"primitive atom {nearest} ({primitive.elements[nearest]}) — "
                f"element mismatch. Check supercell construction."
            )
        magmom.append(primitive_magmom[nearest])
    return magmom


def _magmom_string(magmom_list):
    """Convert a MAGMOM list to compressed VASP format."""
    if not magmom_list:
        return ""
    groups = []
    current = magmom_list[0]
    count = 1
    for val in magmom_list[1:]:
        if val == current:
            count += 1
        else:
            groups.append((count, current))
            current = val
            count = 1
    groups.append((count, current))

    parts = []
    for count, val in groups:
        if val == int(val):
            val_str = str(int(val))
        else:
            val_str = str(val)
        if count == 1:
            parts.append(val_str)
        else:
            parts.append(f"{count}*{val_str}")
    return " ".join(parts)


def _parse_magmom_string(magmom_str):
    """Parse compressed VASP MAGMOM string back to a list of floats.

    Inverse of ``_magmom_string()``.  Handles ``"16*0 4*5 4*-5 80*0"``
    as well as bare values like ``"0 5 -5 0"``.
    """
    result = []
    for token in magmom_str.split():
        if "*" in token:
            count_s, val_s = token.split("*", 1)
            result.extend([float(val_s)] * int(count_s))
        else:
            result.append(float(token))
    return result


def _read_magmom_from_incar(incar_path):
    """Read and parse the MAGMOM line from an INCAR file.

    Returns a list of floats, or None if the file doesn't exist or has
    no MAGMOM tag.
    """
    p = Path(incar_path)
    if not p.exists():
        return None
    for line in p.read_text().splitlines():
        stripped = line.strip()
        if stripped.upper().startswith("MAGMOM"):
            # MAGMOM = 16*0 4*5 ...
            _, _, rhs = stripped.partition("=")
            rhs = rhs.strip()
            if rhs:
                return _parse_magmom_string(rhs)
    return None


# ---------------------------------------------------------------------------
# energies.json helpers
# ---------------------------------------------------------------------------

def _energies_path(supercell_dir: str) -> Path:
    return Path(supercell_dir) / "energies.json"


def _load_energies(supercell_dir: str) -> dict:
    p = _energies_path(supercell_dir)
    if not p.exists():
        raise FileNotFoundError(f"No energies.json in {supercell_dir}")
    return json.loads(p.read_text())


def _save_energies(supercell_dir: str, data: dict) -> None:
    _energies_path(supercell_dir).write_text(
        json.dumps(data, indent=2) + "\n"
    )


# ---------------------------------------------------------------------------
# Step directory parsing
# ---------------------------------------------------------------------------

_STEP_RE = re.compile(r"^step_(\d+)_Li(\d+)$")


def _find_latest_step(supercell_dir: str) -> Optional[Path]:
    """Find the latest step_XX_LiYY directory in supercell_dir."""
    sup = Path(supercell_dir)
    steps = []
    for d in sup.iterdir():
        if d.is_dir():
            m = _STEP_RE.match(d.name)
            if m:
                steps.append((int(m.group(1)), int(m.group(2)), d))
    if not steps:
        return None
    steps.sort(key=lambda x: x[0])
    return steps[-1][2]


def _find_supercell_dir(jid: str, output_dir: Optional[str] = None) -> Optional[str]:
    """Find the supercell_NxNxN directory for a JID.

    Matches directories named exactly ``jid`` or ``jid-<suffix>``
    (e.g. ``JVASP-144791-NMC``), so users can append abbreviations
    for easier navigation without breaking the CLI.
    """
    root = Path(output_dir) if output_dir else _DEFAULT_OUTPUT_DIR
    if not root.exists():
        return None
    # Match directories named exactly `jid` or `jid-*` (abbreviated suffix)
    candidates = [d for d in root.iterdir()
                  if d.is_dir() and (d.name == jid or d.name.startswith(jid + "-"))]
    if not candidates:
        return None
    if len(candidates) > 1:
        raise ValueError(
            f"Multiple directories match {jid}: {[d.name for d in candidates]}. "
            f"Use --output-dir to specify."
        )
    jid_dir = candidates[0]
    for d in jid_dir.iterdir():
        if d.is_dir() and d.name.startswith("supercell_"):
            return str(d)
    return None


# ---------------------------------------------------------------------------
# Load primitive structure
# ---------------------------------------------------------------------------

def _load_primitive(jid: str, output_dir: Optional[str] = None, dft3d_df=None):
    """Load primitive structure from existing POSCAR or JARVIS DB."""
    from jarvis.core.atoms import Atoms as JAtoms

    root = Path(output_dir) if output_dir else _DEFAULT_OUTPUT_DIR

    # Try existing POSCAR in old dir 1 or in any supercell step_00
    for poscar_path in [
        root / jid / "1_pbe_relax_lithiated" / "POSCAR",
    ]:
        if poscar_path.exists():
            from jarvis.io.vasp.inputs import Poscar
            return Poscar.from_file(str(poscar_path)).atoms

    # Try supercell step_00 POSCAR (already a supercell, not primitive)
    # Fall through to DB lookup

    if dft3d_df is not None:
        import pandas as pd
        jid_to_row = {row["jid"]: row for _, row in dft3d_df.iterrows()}
        if jid in jid_to_row:
            return JAtoms.from_dict(jid_to_row[jid]["atoms"])
        return None

    try:
        from jarvis.db.figshare import data as jarvis_data
        import pandas as pd
        dft3d_df = pd.DataFrame(jarvis_data("dft_3d"))
        jid_to_row = {row["jid"]: row for _, row in dft3d_df.iterrows()}
        if jid in jid_to_row:
            return JAtoms.from_dict(jid_to_row[jid]["atoms"])
    except Exception as e:
        print(f"WARNING: Cannot load {jid}: {e}")
    return None


# ---------------------------------------------------------------------------
# Sequential delithiation workflow
# ---------------------------------------------------------------------------

def generate_sequential_init(
    jid: str,
    output_dir: Optional[str] = None,
    dft3d_df=None,
    min_length: float = 7.0,
    max_atoms: int = 300,
    functional: str = "auto",
) -> str:
    """Create supercell directory with step_00 (fully lithiated).

    Args:
        jid:        JARVIS JID.
        output_dir: Root output directory. Defaults to dft_inputs/.
        dft3d_df:   Pre-loaded JARVIS-DFT DataFrame.
        min_length: Minimum lattice vector length in Angstroms.
        max_atoms:  Maximum allowed atoms in supercell. Dimensions are
                    reduced if exceeded, with a warning.
        functional: "auto", "pbe", or "optb88vdw". Auto-detects layered
                    materials and uses optB88-vdW for them.

    Returns:
        Path to the supercell directory.
    """
    root = Path(output_dir) if output_dir else _DEFAULT_OUTPUT_DIR
    root.mkdir(parents=True, exist_ok=True)

    # Load JARVIS DB once if not provided, so both functional resolution
    # and primitive loading can use it.
    if dft3d_df is None:
        try:
            from jarvis.db.figshare import data as jarvis_data
            import pandas as pd
            dft3d_df = pd.DataFrame(jarvis_data("dft_3d"))
        except Exception:
            pass  # Fall through — _load_primitive has its own local-file fallback

    resolved_functional = _resolve_functional(jid, dft3d_df, functional)

    prim_atoms = _load_primitive(jid, output_dir=output_dir, dft3d_df=dft3d_df)
    if prim_atoms is None:
        raise ValueError(f"{jid} not found in JARVIS-DFT dataset or local files.")

    if "Li" not in prim_atoms.elements:
        raise ValueError(f"{jid} has no Li atoms.")

    dim = _min_supercell_dim(prim_atoms, min_length=min_length)

    # Cap supercell size to max_atoms
    n_prim = len(prim_atoms.elements)
    lattice = np.array(prim_atoms.lattice_mat)
    norms_prim = np.linalg.norm(lattice, axis=1)
    reduced = False
    while n_prim * dim[0] * dim[1] * dim[2] > max_atoms and any(d > 1 for d in dim):
        # Reduce the axis with the largest multiplier; break ties by shortest vector
        max_d = max(d for d in dim if d > 1)
        candidates = [i for i, d in enumerate(dim) if d == max_d]
        # Among ties, pick the axis with shortest primitive vector (least impact)
        reduce_i = min(candidates, key=lambda i: norms_prim[i])
        dim[reduce_i] -= 1
        reduced = True

    if reduced:
        actual_atoms = n_prim * dim[0] * dim[1] * dim[2]
        warnings.warn(
            f"{jid}: supercell reduced to {'x'.join(str(d) for d in dim)} "
            f"({actual_atoms} atoms) to stay within max_atoms={max_atoms}. "
            f"Dilute vacancy limit may be compromised."
        )

    dim_str = "x".join(str(d) for d in dim)

    # Build supercell
    if dim == [1, 1, 1]:
        sup_atoms = prim_atoms
    else:
        sup_atoms = prim_atoms.make_supercell(dim)

    n_atoms = len(sup_atoms.elements)
    n_li = sum(1 for el in sup_atoms.elements if el == "Li")
    unique_elements = list(dict.fromkeys(sup_atoms.elements))

    # Generate AFM MAGMOM
    prim_magmom = _afm_magmom(prim_atoms)
    if dim == [1, 1, 1]:
        magmom = prim_magmom
    else:
        magmom = _supercell_magmom(sup_atoms, prim_atoms, prim_magmom, dim)

    # KPOINTS scaled inversely
    sup_kmesh = [max(1, round(3 / d)) for d in dim]
    kmesh_str = " ".join(str(k) for k in sup_kmesh)

    # NSW scales with atom count
    nsw = min(300, max(200, n_atoms * 2))

    # Create directory
    jid_dir = root / jid
    sup_dir = jid_dir / f"supercell_{dim_str}"
    step_dir = sup_dir / f"step_00_Li{n_li}"
    step_dir.mkdir(parents=True, exist_ok=True)

    # Write VASP inputs: ISIF=3 for step 0 (full relaxation)
    _write_poscar(step_dir, sup_atoms)
    write_relax_incar(
        step_dir, elements=unique_elements,
        isif=3, nsw=nsw,
        magmom_str=_magmom_string(magmom),
        functional=resolved_functional,
    )
    write_kpoints(step_dir, mesh=kmesh_str)
    write_potcar_spec(step_dir, unique_elements)

    # Initialize energies.json
    layered = _is_layered(jid, dft3d_df)
    energies_data = {
        "jid": jid,
        "dim": dim,
        "n_li_total": n_li,
        "functional": resolved_functional,
        "layered": layered,
        "e_li_metal": _E_LI_METAL[resolved_functional],
        "steps": [
            {"step": 0, "n_li": n_li, "energy": None, "removed_li_index": None},
        ],
    }
    _save_energies(str(sup_dir), energies_data)

    # Print info
    lattice = np.array(sup_atoms.lattice_mat)
    norms = np.linalg.norm(lattice, axis=1)
    print(f"{jid}: supercell {dim_str} ({n_atoms} atoms, {n_li} Li)")
    print(f"  Functional: {resolved_functional}")
    print(f"  Lattice: {norms[0]:.2f}, {norms[1]:.2f}, {norms[2]:.2f} Å")
    print(f"  KPOINTS: {kmesh_str}")
    print(f"  -> {step_dir}")

    return str(sup_dir)


def generate_next_step(supercell_dir: str) -> Optional[str]:
    """Generate VASP inputs for the next delithiation step.

    Reads CONTCAR from the latest step, removes one Li via ALIGNN-ranked
    vacancy selection, and creates the next step directory.

    Args:
        supercell_dir: Path to supercell_NxNxN directory.

    Returns:
        Path to the new step directory, or None if no Li remain.
    """
    from jarvis.io.vasp.inputs import Poscar

    sup = Path(supercell_dir)
    latest = _find_latest_step(supercell_dir)
    if latest is None:
        raise FileNotFoundError(f"No step directories found in {supercell_dir}")

    m = _STEP_RE.match(latest.name)
    current_step = int(m.group(1))
    current_n_li = int(m.group(2))

    if current_n_li == 0:
        print("No Li remaining — sequential delithiation complete.")
        return None

    # Read relaxed structure from CONTCAR (check results/ subdirectory first)
    contcar = latest / "results" / "CONTCAR"
    if not contcar.exists():
        contcar = latest / "CONTCAR"
    if not contcar.exists():
        raise FileNotFoundError(
            f"No CONTCAR in {latest} or {latest}/results/. "
            f"Run DFT for step {current_step} first."
        )
    atoms = Poscar.from_file(str(contcar)).atoms

    # Verify Li count
    actual_li = sum(1 for el in atoms.elements if el == "Li")
    if actual_li != current_n_li:
        warnings.warn(
            f"Expected {current_n_li} Li in CONTCAR but found {actual_li}."
        )

    if actual_li == 0:
        print("No Li remaining in CONTCAR — sequential delithiation complete.")
        return None

    # Remove one Li via ALIGNN-ranked vacancy selection
    defect_atoms, removed_idx, info = get_next_vacancy(atoms)
    new_n_li = actual_li - 1
    new_step = current_step + 1

    # Create new step directory
    step_dir = sup / f"step_{new_step:02d}_Li{new_n_li}"
    step_dir.mkdir(parents=True, exist_ok=True)

    n_atoms = len(defect_atoms.elements)
    unique_elements = list(dict.fromkeys(defect_atoms.elements))

    nsw = min(300, max(200, n_atoms * 2))

    # Preserve MAGMOM from previous step's INCAR, dropping the removed Li entry
    prev_magmom = _read_magmom_from_incar(latest / "INCAR")
    if prev_magmom is not None and len(prev_magmom) == n_atoms + 1:
        del prev_magmom[removed_idx]
        magmom = prev_magmom
    else:
        if prev_magmom is not None:
            warnings.warn(
                f"Previous INCAR MAGMOM has {len(prev_magmom)} entries, "
                f"expected {n_atoms + 1}. Falling back to _afm_magmom()."
            )
        else:
            warnings.warn(
                "No MAGMOM found in previous INCAR. "
                "Falling back to _afm_magmom()."
            )
        magmom = _afm_magmom(defect_atoms)

    # Read functional and layered flag from energies.json
    data = _load_energies(supercell_dir)
    step_functional = data.get("functional", "pbe")
    # Backward compat: if "layered" key missing, infer from functional
    is_layered = data.get("layered", step_functional == "optb88vdw")
    # ISIF=3 for layered (cell shape changes on delithiation), ISIF=2 otherwise
    isif = 3 if is_layered else 2

    _write_poscar(step_dir, defect_atoms)
    write_relax_incar(
        step_dir, elements=unique_elements,
        isif=isif, nsw=nsw,
        magmom_str=_magmom_string(magmom),
        functional=step_functional,
    )

    # Read KPOINTS mesh from previous step
    prev_kpoints = latest / "KPOINTS"
    if prev_kpoints.exists():
        kp_lines = prev_kpoints.read_text().strip().split("\n")
        mesh = kp_lines[3] if len(kp_lines) > 3 else "3 3 3"
    else:
        mesh = "3 3 3"
    write_kpoints(step_dir, mesh=mesh)
    write_potcar_spec(step_dir, unique_elements)

    # Update energies.json (data already loaded above for functional)
    data["steps"].append({
        "step": new_step,
        "n_li": new_n_li,
        "energy": None,
        "removed_li_index": removed_idx,
    })
    _save_energies(supercell_dir, data)

    # Print info
    print(f"Step {new_step}: {actual_li} -> {new_n_li} Li ({n_atoms} atoms)")
    print(f"  Removed Li atom index: {removed_idx}")
    if info["candidate_energies"]:
        e_str = ", ".join(f"{e:.4f}" for e in info["candidate_energies"])
        print(f"  Candidate energies: [{e_str}]")
    print(f"  -> {step_dir}")

    return str(step_dir)


# ---------------------------------------------------------------------------
# TB-mBJ static generation
# ---------------------------------------------------------------------------

def generate_static(supercell_dir: str, step: int) -> str:
    """Generate TB-mBJ static INCAR for a completed relaxation step.

    Reads the CONTCAR from the specified step and writes TB-mBJ inputs
    into a ``tmbj_step_XX_LiYY/`` directory alongside the relaxation step.

    Args:
        supercell_dir: Path to supercell_NxNxN directory.
        step:          Step number whose relaxed structure to use.

    Returns:
        Path to the TB-mBJ directory.
    """
    from jarvis.io.vasp.inputs import Poscar

    sup = Path(supercell_dir)
    data = _load_energies(supercell_dir)

    # Find the step entry
    step_entry = None
    for s in data["steps"]:
        if s["step"] == step:
            step_entry = s
            break
    if step_entry is None:
        raise ValueError(
            f"Step {step} not found in energies.json. "
            f"Available steps: {[s['step'] for s in data['steps']]}"
        )

    n_li = step_entry["n_li"]
    relax_dir = sup / f"step_{step:02d}_Li{n_li}"
    if not relax_dir.exists():
        raise FileNotFoundError(f"Relaxation directory not found: {relax_dir}")

    # Read relaxed structure from CONTCAR
    contcar = relax_dir / "results" / "CONTCAR"
    if not contcar.exists():
        contcar = relax_dir / "CONTCAR"
    if not contcar.exists():
        raise FileNotFoundError(
            f"No CONTCAR in {relax_dir} or {relax_dir}/results/. "
            f"Run DFT relaxation for step {step} first."
        )
    atoms = Poscar.from_file(str(contcar)).atoms

    unique_elements = list(dict.fromkeys(atoms.elements))

    # Preserve MAGMOM from relaxation step's INCAR (same atom count)
    prev_magmom = _read_magmom_from_incar(relax_dir / "INCAR")
    if prev_magmom is not None and len(prev_magmom) == len(atoms.elements):
        magmom = prev_magmom
    else:
        magmom = _afm_magmom(atoms)

    # Create TB-mBJ directory
    tmbj_dir = sup / f"tmbj_step_{step:02d}_Li{n_li}"
    tmbj_dir.mkdir(parents=True, exist_ok=True)

    _write_poscar(tmbj_dir, atoms)
    write_tmbj_incar(tmbj_dir, elements=unique_elements,
                     magmom_str=_magmom_string(magmom))

    # Read KPOINTS mesh from relaxation step
    prev_kpoints = relax_dir / "KPOINTS"
    if prev_kpoints.exists():
        kp_lines = prev_kpoints.read_text().strip().split("\n")
        mesh = kp_lines[3] if len(kp_lines) > 3 else "3 3 3"
    else:
        mesh = "3 3 3"
    write_kpoints(tmbj_dir, mesh=mesh)
    write_potcar_spec(tmbj_dir, unique_elements)

    print(f"TB-mBJ static for step {step} (Li{n_li})")
    print(f"  -> {tmbj_dir}")

    return str(tmbj_dir)


# ---------------------------------------------------------------------------
# Energy recording
# ---------------------------------------------------------------------------

def _read_outcar_energy(step_dir: Path) -> Optional[float]:
    """Parse final total energy from OUTCAR in results/ subdirectory."""
    outcar = step_dir / "results" / "OUTCAR"
    if not outcar.exists():
        return None
    energy = None
    for line in outcar.open():
        if "free  energy   TOTEN" in line:
            try:
                energy = float(line.split()[-2])
            except (IndexError, ValueError):
                pass
    return energy


def record_energy(
    supercell_dir: str,
    step: int,
    energy: Optional[float] = None,
) -> dict:
    """Record a DFT total energy for a completed step.

    Args:
        supercell_dir: Path to supercell_NxNxN directory.
        step:          Step number (0, 1, 2, ...).
        energy:        DFT total energy in eV. If None, reads from
                       results/OUTCAR in the step directory.

    Returns:
        Updated energies dict.
    """
    data = _load_energies(supercell_dir)

    if energy is None:
        # Auto-read from results/OUTCAR
        n_li = None
        for entry in data["steps"]:
            if entry["step"] == step:
                n_li = entry["n_li"]
                break
        if n_li is None:
            raise ValueError(
                f"Step {step} not found in energies.json. "
                f"Available steps: {[s['step'] for s in data['steps']]}"
            )
        step_dir = None
        for d in Path(supercell_dir).iterdir():
            if d.is_dir() and d.name == f"step_{step:02d}_Li{n_li}":
                step_dir = d
                break
        if step_dir is None:
            raise FileNotFoundError(f"No directory found for step {step} in {supercell_dir}")
        energy = _read_outcar_energy(step_dir)
        if energy is None:
            raise FileNotFoundError(
                f"No results/OUTCAR found in {step_dir}. "
                f"Supply energy explicitly or place DFT outputs in results/."
            )
    found = False
    for entry in data["steps"]:
        if entry["step"] == step:
            entry["energy"] = energy
            found = True
            break
    if not found:
        raise ValueError(
            f"Step {step} not found in energies.json. "
            f"Available steps: {[s['step'] for s in data['steps']]}"
        )
    _save_energies(supercell_dir, data)
    print(f"Recorded energy for step {step}: {energy:.6f} eV")
    return data


# ---------------------------------------------------------------------------
# Voltage curve computation
# ---------------------------------------------------------------------------

def compute_voltage_curve(
    supercell_dir: str,
    plot: bool = True,
    e_li_metal: Optional[float] = None,
) -> List[Dict]:
    """Compute voltage curve from recorded energies.

    Calculates both step voltages (raw, between consecutive steps) and
    convex hull equilibrium voltages.

    Args:
        supercell_dir: Path to supercell_NxNxN directory.
        plot:          If True, save voltage_curve.png.
        e_li_metal:    DFT energy of Li metal per atom (eV/atom).
                       If None, reads from energies.json.

    Returns:
        List of dicts with keys: x (Li fraction), n_li, voltage.
    """
    data = _load_energies(supercell_dir)
    n_li_total = data["n_li_total"]

    # Resolve e_li_metal: CLI arg > energies.json > default PBE
    if e_li_metal is None:
        e_li_metal = data.get("e_li_metal")
    if e_li_metal is None:
        raise ValueError(
            f"No e_li_metal available. For optB88-vdW, supply --e-li-metal "
            f"once the Li metal reference energy is computed."
        )

    # Collect steps with recorded energies, sort by n_li descending
    recorded = [s for s in data["steps"] if s["energy"] is not None]
    if len(recorded) < 2:
        raise ValueError(
            f"Need at least 2 recorded energies, have {len(recorded)}."
        )
    recorded.sort(key=lambda s: s["n_li"], reverse=True)

    # Step voltages: V = (E_low - E_high) / delta_n + e_li_metal
    # Derivation: reaction Li_lo + dn*Li(metal) -> Li_hi, dE = E_hi - E_lo - dn*mu_Li,
    # V = -dE/dn = (E_lo - E_hi)/dn + mu_Li. With mu_Li negative, this lowers V.
    results = []
    for i in range(len(recorded) - 1):
        hi = recorded[i]
        lo = recorded[i + 1]
        dn = hi["n_li"] - lo["n_li"]
        if dn <= 0:
            continue
        v = (lo["energy"] - hi["energy"]) / dn + e_li_metal
        x = lo["n_li"] / n_li_total
        results.append({
            "x": x,
            "n_li": lo["n_li"],
            "voltage": v,
            "n_li_from": hi["n_li"],
        })

    # Convex hull on formation energies
    # dE(x) = E(x) - x*E(fully_lith) - (1-x)*E(fully_delith)
    # where x = n_li / n_li_total
    e_full = None
    e_empty = None
    for s in recorded:
        if s["n_li"] == n_li_total:
            e_full = s["energy"]
        if s["n_li"] == 0:
            e_empty = s["energy"]

    hull_voltages = []
    if e_full is not None and e_empty is not None:
        # Compute formation energies
        points = []  # (x, E, formation_energy)
        for s in recorded:
            x = s["n_li"] / n_li_total
            fe = s["energy"] - x * e_full - (1 - x) * e_empty
            points.append((x, s["n_li"], s["energy"], fe))
        points.sort(key=lambda p: p[0], reverse=True)

        # Lower convex hull: start from x=1, greedily pick the point that
        # gives the lowest (most negative) formation energy per unit x
        hull = [points[0]]  # x=1 (fully lithiated)
        for p in points[1:]:
            # Remove points that are above the line from hull[-1] to p
            while len(hull) > 1:
                # Check if hull[-1] is above line from hull[-2] to p
                x0, _, _, fe0 = hull[-2]
                x1, _, _, fe1 = hull[-1]
                x2, _, _, fe2 = p
                if x0 == x2:
                    break
                # Linear interpolation of fe at x1 between x0 and x2
                t = (x0 - x1) / (x0 - x2)
                fe_interp = fe0 + t * (fe2 - fe0)
                if fe1 <= fe_interp + 1e-10:
                    break
                hull.pop()
            hull.append(p)

        # Hull voltages connect consecutive hull vertices
        for i in range(len(hull) - 1):
            x_hi, nli_hi, e_hi, _ = hull[i]
            x_lo, nli_lo, e_lo, _ = hull[i + 1]
            dn = nli_hi - nli_lo
            if dn <= 0:
                continue
            v = (e_lo - e_hi) / dn + e_li_metal
            hull_voltages.append({
                "x_from": x_hi,
                "x_to": x_lo,
                "voltage": v,
            })

    if plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            # Extract abbreviation from directory name (e.g. JVASP-42723-LFP -> LFP)
            _jid_dir_name = Path(supercell_dir).parent.name
            _abbrev_parts = _jid_dir_name.split("-")[2:]  # after JVASP-XXXXX
            _abbrev = "-".join(_abbrev_parts) if _abbrev_parts else ""
            _title_suffix = f" ({_abbrev})" if _abbrev else ""

            # --- discharge_curve.png: staircase (step) plot ---
            fig, ax = plt.subplots(figsize=(8, 5))
            if results:
                xs = [1.0]
                vs = [results[0]["voltage"]]
                for r in results:
                    xs.extend([r["n_li_from"] / n_li_total, r["x"]])
                    vs.extend([r["voltage"], r["voltage"]])
                ax.plot(xs, vs, "b-", linewidth=1.5, label="Step voltage")
            if hull_voltages:
                hx = []
                hv = []
                for hv_entry in hull_voltages:
                    hx.extend([hv_entry["x_from"], hv_entry["x_to"]])
                    hv.extend([hv_entry["voltage"], hv_entry["voltage"]])
                ax.plot(hx, hv, "r--", linewidth=2, label="Equilibrium (hull)")
            ax.set_xlabel("x in Li$_x$MO")
            ax.set_ylabel("Voltage (V)")
            ax.set_title(f"Discharge curve — {data['jid']}{_title_suffix}")
            ax.legend()
            ax.set_xlim(-0.05, 1.05)
            fig.tight_layout()
            analysis_dir = Path(supercell_dir).parents[2] / "analysis"
            if analysis_dir.exists():
                fig.savefig(str(analysis_dir / f"discharge_curve_{data['jid']}.png"), dpi=150)
                print(f"Saved discharge_curve_{data['jid']}.png in {analysis_dir}")
            else:
                fig.savefig(str(Path(supercell_dir) / "discharge_curve.png"), dpi=150)
                print(f"Saved discharge_curve.png in {supercell_dir}")
            plt.close(fig)

            # --- voltage_curve.png: line plot (raw voltages at step midpoints) ---
            fig, ax = plt.subplots(figsize=(8, 5))
            if results:
                xs = [(r["n_li_from"] + r["n_li"]) / (2 * n_li_total) for r in results]
                vs = [r["voltage"] for r in results]
                ax.plot(xs, vs, "b-o", linewidth=1.5, markersize=5, label="Step voltage")
            if hull_voltages:
                label = "Equilibrium (hull)"
                for h in hull_voltages:
                    ax.plot([h["x_to"], h["x_from"]], [h["voltage"], h["voltage"]],
                            "r--", linewidth=2, label=label)
                    label = None
            ax.set_xlabel("x in Li$_x$MO")
            ax.set_ylabel("Voltage (V)")
            ax.set_title(f"Voltage curve — {data['jid']}{_title_suffix}")
            ax.legend()
            ax.set_xlim(-0.05, 1.05)
            fig.tight_layout()
            if analysis_dir.exists():
                fig.savefig(str(analysis_dir / f"voltage_curve_{data['jid']}.png"), dpi=150)
                print(f"Saved voltage_curve_{data['jid']}.png in {analysis_dir}")
            else:
                fig.savefig(str(Path(supercell_dir) / "voltage_curve.png"), dpi=150)
                print(f"Saved voltage_curve.png in {supercell_dir}")
            plt.close(fig)
        except ImportError:
            print("matplotlib not available — skipping plot.")

    return results


# ---------------------------------------------------------------------------
# Deprecated: old workflows
# ---------------------------------------------------------------------------

def generate_supercell_inputs(*args, **kwargs):
    """Deprecated: use generate_sequential_init() instead."""
    raise DeprecationWarning(
        "generate_supercell_inputs() is removed. "
        "Use: python dft_prep.py init JVASP-XXXXX"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cli():
    import argparse

    parser = argparse.ArgumentParser(
        description="Sequential supercell delithiation for voltage curves."
    )
    sub = parser.add_subparsers(dest="command")

    # init
    p_init = sub.add_parser("init", help="Initialize supercell + step_00")
    p_init.add_argument("jids", nargs="+", help="JARVIS JIDs")
    p_init.add_argument("--output-dir", default=None)
    p_init.add_argument("--min-length", type=float, default=7.0)
    p_init.add_argument("--max-atoms", type=int, default=300,
                        help="Max atoms in supercell (default: 300)")
    p_init.add_argument("--functional", choices=["auto", "pbe", "optb88vdw"],
                        default="auto",
                        help="DFT functional (default: auto-detect)")

    # next
    p_next = sub.add_parser("next", help="Generate next step from CONTCAR")
    p_next.add_argument("jid", help="JARVIS JID")
    p_next.add_argument("--output-dir", default=None)

    # record
    p_rec = sub.add_parser("record", help="Record DFT energy for a step")
    p_rec.add_argument("jid", help="JARVIS JID")
    p_rec.add_argument("step", type=int, help="Step number")
    p_rec.add_argument("energy", type=float, nargs="?", default=None,
                        help="DFT total energy (eV). If omitted, reads from results/OUTCAR")
    p_rec.add_argument("--output-dir", default=None)

    # static
    p_static = sub.add_parser("static", help="Generate TB-mBJ static inputs for a step")
    p_static.add_argument("jid", help="JARVIS JID")
    p_static.add_argument("steps", type=int, nargs="+",
                          help="Step number(s) to generate TB-mBJ inputs for")
    p_static.add_argument("--output-dir", default=None)

    # voltage
    p_volt = sub.add_parser("voltage", help="Compute + plot voltage curve")
    p_volt.add_argument("jid", help="JARVIS JID")
    p_volt.add_argument("--output-dir", default=None)
    p_volt.add_argument("--no-plot", action="store_true")
    p_volt.add_argument("--e-li-metal", type=float, default=None,
                        help="Li metal energy (eV/atom). Default: read from energies.json")

    args = parser.parse_args()

    if args.command == "init":
        for jid in args.jids:
            generate_sequential_init(
                jid, output_dir=args.output_dir,
                min_length=args.min_length,
                max_atoms=args.max_atoms,
                functional=args.functional,
            )

    elif args.command == "next":
        sup_dir = _find_supercell_dir(args.jid, args.output_dir)
        if sup_dir is None:
            print(f"ERROR: No supercell directory found for {args.jid}. Run 'init' first.")
            sys.exit(1)
        generate_next_step(sup_dir)

    elif args.command == "record":
        sup_dir = _find_supercell_dir(args.jid, args.output_dir)
        if sup_dir is None:
            print(f"ERROR: No supercell directory found for {args.jid}.")
            sys.exit(1)
        record_energy(sup_dir, args.step, args.energy)

    elif args.command == "static":
        sup_dir = _find_supercell_dir(args.jid, args.output_dir)
        if sup_dir is None:
            print(f"ERROR: No supercell directory found for {args.jid}.")
            sys.exit(1)
        for step in args.steps:
            generate_static(sup_dir, step)

    elif args.command == "voltage":
        sup_dir = _find_supercell_dir(args.jid, args.output_dir)
        if sup_dir is None:
            print(f"ERROR: No supercell directory found for {args.jid}.")
            sys.exit(1)
        results = compute_voltage_curve(
            sup_dir, plot=not args.no_plot,
            e_li_metal=args.e_li_metal,
        )
        for r in results:
            print(f"  x={r['x']:.3f}  n_li={r['n_li']}  V={r['voltage']:.4f}")

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    _cli()
