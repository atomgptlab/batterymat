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
    z: int = 1,
) -> float:
    """Average intercalation voltage from two endpoint DFT total energies.

    Same sign convention as ``compute_voltage_curve()``:

        V = (E_delithiated - E_lithiated) / (z * n_ion) + E_ion_metal / z

    Derivation: insertion reaction  host + n M(metal) -> M_n host,
    dE = E_lith - E_delith - n*E_metal,  V = -dE / (z*n).  With E_metal < 0
    (e.g. Li PBE -1.9031 eV/atom) the metal reference *lowers* the voltage.
    ``z`` is the number of electrons per ion (1 for Li/Na/K, 2 for Mg/Ca/Zn).
    """
    if n_li == 0:
        raise ValueError("n_li must be > 0")
    return (e_delithiated - e_lithiated) / (z * n_li) + e_li_metal / z


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

# Nonlocal vdW-DF tags (VASP wiki, "Nonlocal vdW-DF functionals").
# optB88-vdW is GGA = BO with PARAM1/PARAM2 (the JARVIS-DFT functional).
# GGA = OR is optPBE-vdW: it was used, mislabelled as optB88-vdW, for the
# original JVASP-2017-LCO chain and Li_sv reference; kept as "optpbevdw".
_VDW_TAGS = {
    "optb88vdw": """\
GGA = BO
PARAM1 = 0.1833333333
PARAM2 = 0.2200000000
LUSE_VDW = .TRUE.
AGGAC = 0.0
""",
    "optpbevdw": """\
GGA = OR
LUSE_VDW = .TRUE.
AGGAC = 0.0
""",
}
_OPTB88VDW_TAGS = _VDW_TAGS["optb88vdw"]

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
    # See dft_inputs/JVASP-913-Li/Li_sv_PBE/, Li_sv_optPBEvdW/ (GGA=OR) and
    # Li_sv_optB88vdW/ (GGA=BO + PARAM1/2, confirmation rerun of Li_sv/static).
    # The old JARVIS value (-0.925) used a different PAW/basis and was wrong by ~1 eV.
    "pbe": -1.9031,
    "optpbevdw": -0.9646,
    "optb88vdw": -0.9778,   # JVASP-913-Li/Li_sv/static (GGA=BO, PARAM1/2); confirmed -0.97791 by Li_sv_optB88vdW/ (2026-09-29)
}

# ---------------------------------------------------------------------------
# Multi-ion support (--ion). Li is the default and keeps the legacy
# energies.json schema (n_li_total / e_li_metal / n_li) byte-for-byte.
# ---------------------------------------------------------------------------
_ION_CHARGE = {"Li": 1, "Na": 1, "K": 1, "Mg": 2, "Ca": 2, "Zn": 2, "Al": 3}

_E_ION_METAL = {
    "Li": _E_LI_METAL,
    # Na (bcc, Na_pv) and Mg (hcp, Mg_pv) references are computed with the
    # inputs in dft_inputs/ref-Na/ and dft_inputs/ref-Mg/ using the same
    # protocol as JVASP-913-Li (ENCUT=520, ISMEAR=1, SIGMA=0.1, dense k).
    # Fill these in once the runs finish, or pass --e-ion-metal.
    # atomgptlab 2026-09-29 (lowest TOTEN / atoms). Na PBE, Na optB88 and Mg PBE stopped
    # with "ZBRENT: fatal error in bracketing" at the EDIFFG=-1E-3 noise floor; their last
    # energies agree to <0.2 meV and |P| < 2 kB, so they are converged for voltages.
    "Na": {"pbe": -1.31054, "optb88vdw": 0.91960, "optpbevdw": None},
    "Mg": {"pbe": -1.50576, "optb88vdw": 1.10509, "optpbevdw": None},
}

_FARADAY_MAH_PER_MOL = 96485.33212 / 3.6   # mAh per mole of electrons
_AVOGADRO = 6.02214076e23

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
    if functional in ("pbe", "optb88vdw", "optpbevdw"):
        return functional
    raise ValueError(f"Unknown functional: {functional!r}. Use 'auto', 'pbe', 'optb88vdw' or 'optpbevdw'.")


# ---------------------------------------------------------------------------
# INCAR / KPOINTS / POTCAR helpers
# ---------------------------------------------------------------------------

# Relative cell-volume change between a step's input POSCAR and its CONTCAR
# above which `next` emits a warning (see generate_next_step).
VOLUME_CHANGE_WARN = 0.10


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
        functional: "pbe", "optb88vdw" (appends GGA=BO, PARAM1, PARAM2,
                    LUSE_VDW, AGGAC) or "optpbevdw" (GGA=OR, LUSE_VDW, AGGAC).
    """
    content = _PBE_INCAR
    content = content.replace("ISIF = 3", f"ISIF = {isif}")
    content = content.replace("NSW = 100", f"NSW = {nsw}")
    if functional in _VDW_TAGS:
        content += _VDW_TAGS[functional]
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

def get_next_vacancy(atoms, ion: str = "Li") -> Tuple:
    """Return the lowest-energy single-ion-vacancy structure with metadata.

    Evaluates ALIGNN energy for every vacancy of ``ion`` (no Wyckoff
    deduplication) to ensure the best site is selected even after symmetry
    is broken by prior vacancies.

    Args:
        atoms: jarvis.core.atoms.Atoms object.
        ion:   Element symbol of the working ion (default "Li").

    Returns:
        (defect_structure, removed_atom_index, info_dict)
        info_dict has keys: 'n_candidates', 'candidate_energies', 'selected_energy'
    """
    li_indices = [i for i, el in enumerate(atoms.elements) if el == ion]
    n_li = len(li_indices)

    if n_li == 0:
        raise ValueError(f"Structure has no {ion} atoms — cannot create vacancy.")

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


def theoretical_capacity(atoms, ion: str = "Li", z: Optional[int] = None) -> Tuple[float, float]:
    """Gravimetric (mAh/g) and volumetric (mAh/cm^3) capacity for removing every ``ion``.

    Q_grav = n_ion * z * F / (3.6 * M)          (M = cell mass, g/mol)
    Q_vol  = n_ion * z * F / (3.6 * V * N_A)    (V = cell volume, cm^3)
    z electrons per ion: 1 for Li/Na/K, 2 for Mg/Ca/Zn (see _ION_CHARGE).
    """
    if z is None:
        z = _ION_CHARGE[ion]
    n_ion = sum(1 for el in atoms.elements if el == ion)
    mass = atoms.composition.weight            # g/mol for the whole cell
    vol_cm3 = atoms.volume * 1e-24             # A^3 -> cm^3
    q_grav = n_ion * z * _FARADAY_MAH_PER_MOL / mass
    q_vol = n_ion * z * _FARADAY_MAH_PER_MOL / (vol_cm3 * _AVOGADRO)
    return q_grav, q_vol


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


def _ion_info(data: dict) -> Tuple[str, int, int, Optional[float]]:
    """Return (ion, z, n_ion_total, e_ion_metal) from an energies.json dict.

    Accepts both the legacy Li schema (``n_li_total`` / ``e_li_metal``) and
    the multi-ion schema (``ion`` / ``z`` / ``n_ion_total`` / ``e_ion_metal``).
    """
    ion = data.get("ion", "Li")
    z = int(data.get("z", _ION_CHARGE.get(ion, 1)))
    n_total = data.get("n_ion_total", data.get("n_li_total"))
    e_metal = data.get("e_ion_metal", data.get("e_li_metal"))
    return ion, z, n_total, e_metal


def _step_n(entry: dict) -> int:
    """Ion count of a step entry (``n_ion`` for non-Li, legacy ``n_li`` for Li)."""
    return entry["n_ion"] if "n_ion" in entry else entry["n_li"]


def _step_entry(data: dict, step: int, n: int, removed_idx) -> dict:
    """Build a step entry using the schema already used by ``data``."""
    if "ion" in data:
        return {"step": step, "n_ion": n, "energy": None, "removed_ion_index": removed_idx}
    return {"step": step, "n_li": n, "energy": None, "removed_li_index": removed_idx}


# ---------------------------------------------------------------------------
# Step directory parsing
# ---------------------------------------------------------------------------

# step_XX_<Ion>YY, e.g. step_00_Li16, step_03_Na5, step_00_Mg8
_STEP_RE = re.compile(r"^step_(\d+)_([A-Z][a-z]?)(\d+)$")


def _find_latest_step(supercell_dir: str) -> Optional[Path]:
    """Find the latest step_XX_<Ion>YY directory in supercell_dir."""
    sup = Path(supercell_dir)
    steps = []
    for d in sup.iterdir():
        if d.is_dir():
            m = _STEP_RE.match(d.name)
            if m:
                steps.append((int(m.group(1)), int(m.group(3)), d))
    if not steps:
        return None
    steps.sort(key=lambda x: x[0])
    return steps[-1][2]


def _find_supercell_dir(jid: str, output_dir: Optional[str] = None) -> Optional[str]:
    """Find the supercell_NxNxN directory for a JID.

    Matches directories named exactly ``jid`` or ``jid-<suffix>``
    (e.g. ``JVASP-144791-NMC``), so users can append abbreviations
    for easier navigation without breaking the CLI.

    An exact directory-name match always wins, so when several campaigns
    share a JID (``JVASP-2017-LCO`` and ``JVASP-2017-LCO-PBE``) pass the
    full directory name instead of the bare JID.
    """
    root = Path(output_dir) if output_dir else _DEFAULT_OUTPUT_DIR
    if not root.exists():
        return None
    # Exact match first (full directory name given, e.g. JVASP-2017-LCO-PBE)
    exact = root / jid
    if exact.is_dir():
        candidates = [exact]
    else:
        # Otherwise match `jid-*` (abbreviated suffix)
        candidates = [d for d in root.iterdir()
                      if d.is_dir() and d.name.startswith(jid + "-")]
    if not candidates:
        return None
    if len(candidates) > 1:
        raise ValueError(
            f"Multiple directories match {jid}: {sorted(d.name for d in candidates)}. "
            f"Pass the full directory name (e.g. {sorted(d.name for d in candidates)[0]}) "
            f"instead of the bare JID."
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


def _load_dft3d_df(db_cache: Optional[str] = None):
    """Return the JARVIS dft_3d table as a pandas DataFrame.

    ``db_cache`` (or the ``BATTERYMAT_DFT3D_CACHE`` environment variable)
    points to a local snapshot of dft_3d: a pickle of the list of dicts
    (``.pkl``/``.pickle``) or a JSON file. When neither is set, the table is
    downloaded from figshare via ``jarvis.db.figshare.data("dft_3d")``.
    """
    import pandas as pd
    db_cache = db_cache or os.environ.get("BATTERYMAT_DFT3D_CACHE")
    if db_cache:
        p = Path(db_cache)
        if not p.exists():
            raise FileNotFoundError(f"JARVIS dft_3d cache not found: {p}")
        if p.suffix in (".pkl", ".pickle"):
            import pickle
            with p.open("rb") as fh:
                rows = pickle.load(fh)
        else:
            rows = json.loads(p.read_text())
        return pd.DataFrame(rows)
    from jarvis.db.figshare import data as jarvis_data
    return pd.DataFrame(jarvis_data("dft_3d"))


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
    ion: str = "Li",
    suffix: Optional[str] = None,
    db_cache: Optional[str] = None,
    e_ion_metal: Optional[float] = None,
) -> str:
    """Create supercell directory with step_00 (fully intercalated).

    Args:
        jid:        JARVIS JID.
        output_dir: Root output directory. Defaults to dft_inputs/.
        dft3d_df:   Pre-loaded JARVIS-DFT DataFrame.
        min_length: Minimum lattice vector length in Angstroms.
        max_atoms:  Maximum allowed atoms in supercell. Dimensions are
                    reduced if exceeded, with a warning.
        functional: "auto", "pbe", or "optb88vdw". Auto-detects layered
                    materials and uses optB88-vdW for them.
        ion:        Working ion to remove ("Li" default; "Na", "Mg", ...).
        suffix:     Optional directory suffix: writes ``<jid>-<suffix>/``
                    (e.g. ``JVASP-2017-LCO-PBE``) instead of ``<jid>/``.
        db_cache:   Local dft_3d snapshot (pickle/JSON) used instead of the
                    figshare download. Also read from BATTERYMAT_DFT3D_CACHE.
        e_ion_metal: Metal reference energy (eV/atom) to store in
                    energies.json; overrides the built-in table.

    Returns:
        Path to the supercell directory.
    """
    root = Path(output_dir) if output_dir else _DEFAULT_OUTPUT_DIR
    root.mkdir(parents=True, exist_ok=True)

    if ion not in _ION_CHARGE:
        raise ValueError(f"Unsupported ion {ion!r}. Known: {sorted(_ION_CHARGE)}")

    # Load JARVIS DB once if not provided, so both functional resolution
    # and primitive loading can use it.
    if dft3d_df is None:
        if db_cache or os.environ.get("BATTERYMAT_DFT3D_CACHE"):
            dft3d_df = _load_dft3d_df(db_cache)   # local cache: let errors surface
        else:
            try:
                dft3d_df = _load_dft3d_df()
            except Exception:
                pass  # Fall through — _load_primitive has its own local-file fallback

    resolved_functional = _resolve_functional(jid, dft3d_df, functional)

    prim_atoms = _load_primitive(jid, output_dir=output_dir, dft3d_df=dft3d_df)
    if prim_atoms is None:
        raise ValueError(f"{jid} not found in JARVIS-DFT dataset or local files.")

    if ion not in prim_atoms.elements:
        raise ValueError(f"{jid} has no {ion} atoms.")

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
    n_li = sum(1 for el in sup_atoms.elements if el == ion)
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
    jid_dir = root / (f"{jid}-{suffix}" if suffix else jid)
    sup_dir = jid_dir / f"supercell_{dim_str}"
    step_dir = sup_dir / f"step_00_{ion}{n_li}"
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
    if e_ion_metal is None:
        e_ion_metal = _E_ION_METAL.get(ion, {}).get(resolved_functional)
    if ion == "Li":
        # Legacy schema, kept byte-identical for Li
        energies_data = {
            "jid": jid,
            "dim": dim,
            "n_li_total": n_li,
            "functional": resolved_functional,
            "layered": layered,
            "e_li_metal": e_ion_metal,
            "steps": [
                {"step": 0, "n_li": n_li, "energy": None, "removed_li_index": None},
            ],
        }
    else:
        energies_data = {
            "jid": jid,
            "ion": ion,
            "z": _ION_CHARGE[ion],
            "dim": dim,
            "n_ion_total": n_li,
            "functional": resolved_functional,
            "layered": layered,
            "e_ion_metal": e_ion_metal,   # None until the ref-<ion> run is done
            "steps": [
                {"step": 0, "n_ion": n_li, "energy": None, "removed_ion_index": None},
            ],
        }
    if suffix:
        # Tag used to disambiguate plot filenames (voltage_curve_<jid>-<tag>.png)
        energies_data["tag"] = suffix
    _save_energies(str(sup_dir), energies_data)

    # Print info
    lattice = np.array(sup_atoms.lattice_mat)
    norms = np.linalg.norm(lattice, axis=1)
    q_grav, q_vol = theoretical_capacity(sup_atoms, ion=ion)
    print(f"{jid}: supercell {dim_str} ({n_atoms} atoms, {n_li} {ion})")
    print(f"  Functional: {resolved_functional}")
    print(f"  Lattice: {norms[0]:.2f}, {norms[1]:.2f}, {norms[2]:.2f} Å")
    print(f"  KPOINTS: {kmesh_str}")
    print(f"  Capacity (all {ion} removed, z={_ION_CHARGE[ion]}): "
          f"{q_grav:.1f} mAh/g, {q_vol:.1f} mAh/cm^3")
    if e_ion_metal is None:
        print(f"  WARNING: no {ion} metal reference for {resolved_functional}; "
              f"run dft_inputs/ref-{ion}/ and pass --e-ion-metal to 'voltage'.")
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

    # Read ion / functional / layered flag from energies.json
    data = _load_energies(supercell_dir)
    ion, _, _, _ = _ion_info(data)

    m = _STEP_RE.match(latest.name)
    current_step = int(m.group(1))
    current_n_li = int(m.group(3))

    if current_n_li == 0:
        print(f"No {ion} remaining — sequential delithiation complete.")
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

    # Volume-change check: a relaxation that changes the cell volume by more
    # than VOLUME_CHANGE_WARN (fraction) is flagged for inspection. This does
    # not stop the chain; it records the condition that required a manual
    # ISIF=2 rerun on the delithiated LiCoO2 endpoint.
    poscar_in = latest / "POSCAR"
    if poscar_in.exists():
        v_in = Poscar.from_file(str(poscar_in)).atoms.volume
        dv = atoms.volume / v_in - 1.0
        if abs(dv) > VOLUME_CHANGE_WARN:
            warnings.warn(
                f"Step {current_step}: cell volume changed by {100*dv:+.1f}% during "
                f"relaxation (>{100*VOLUME_CHANGE_WARN:.0f}%). Inspect the CONTCAR; "
                f"consider rerunning this step with ISIF=2 before continuing."
            )

    # Verify ion count
    actual_li = sum(1 for el in atoms.elements if el == ion)
    if actual_li != current_n_li:
        warnings.warn(
            f"Expected {current_n_li} {ion} in CONTCAR but found {actual_li}."
        )

    if actual_li == 0:
        print(f"No {ion} remaining in CONTCAR — sequential delithiation complete.")
        return None

    # Remove one ion via ALIGNN-ranked vacancy selection
    defect_atoms, removed_idx, info = get_next_vacancy(atoms, ion=ion)
    new_n_li = actual_li - 1
    new_step = current_step + 1

    # Create new step directory
    step_dir = sup / f"step_{new_step:02d}_{ion}{new_n_li}"
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

    step_functional = data.get("functional", "pbe")
    # Backward compat: if "layered" key missing, infer from functional
    is_layered = data.get("layered", step_functional in _VDW_TAGS)
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
    data["steps"].append(_step_entry(data, new_step, new_n_li, removed_idx))
    _save_energies(supercell_dir, data)

    # Print info
    print(f"Step {new_step}: {actual_li} -> {new_n_li} {ion} ({n_atoms} atoms)")
    print(f"  Removed {ion} atom index: {removed_idx}")
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
    ion, _, _, _ = _ion_info(data)

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

    n_li = _step_n(step_entry)
    relax_dir = sup / f"step_{step:02d}_{ion}{n_li}"
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
    tmbj_dir = sup / f"tmbj_step_{step:02d}_{ion}{n_li}"
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

    print(f"TB-mBJ static for step {step} ({ion}{n_li})")
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
    ion, _, _, _ = _ion_info(data)

    if energy is None:
        # Auto-read from results/OUTCAR
        n_li = None
        for entry in data["steps"]:
            if entry["step"] == step:
                n_li = _step_n(entry)
                break
        if n_li is None:
            raise ValueError(
                f"Step {step} not found in energies.json. "
                f"Available steps: {[s['step'] for s in data['steps']]}"
            )
        step_dir = None
        for d in Path(supercell_dir).iterdir():
            if d.is_dir() and d.name == f"step_{step:02d}_{ion}{n_li}":
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
    e_ion_metal: Optional[float] = None,
) -> List[Dict]:
    """Compute voltage curve from recorded energies.

    Calculates both step voltages (raw, between consecutive steps) and
    convex hull equilibrium voltages.

    Args:
        supercell_dir: Path to supercell_NxNxN directory.
        plot:          If True, save voltage_curve.png.
        e_li_metal:    DFT energy of the working-ion metal per atom
                       (eV/atom). If None, reads from energies.json.
        e_ion_metal:   Alias of ``e_li_metal`` for non-Li ions.

    Returns:
        List of dicts with keys: x (ion fraction), n_li (ion count after the
        step), n_li_from (ion count before), voltage. The ``n_li`` key names
        are kept for backward compatibility; they hold the working-ion count.
    """
    data = _load_energies(supercell_dir)
    ion, z, n_li_total, e_metal_json = _ion_info(data)

    # Resolve metal reference: CLI arg > energies.json
    if e_ion_metal is not None:
        e_li_metal = e_ion_metal
    if e_li_metal is None:
        e_li_metal = e_metal_json
    if e_li_metal is None:
        e_li_metal = _E_ION_METAL.get(ion, {}).get(data.get("functional", "pbe"))
    if e_li_metal is None:
        raise ValueError(
            f"No {ion} metal reference energy available. Supply --e-ion-metal "
            f"(or --e-li-metal) once the {ion} metal reference run "
            f"(dft_inputs/ref-{ion}/) is computed."
        )

    # Collect steps with recorded energies, sort by ion count descending
    recorded = [s for s in data["steps"] if s["energy"] is not None]
    if len(recorded) < 2:
        raise ValueError(
            f"Need at least 2 recorded energies, have {len(recorded)}."
        )
    recorded.sort(key=_step_n, reverse=True)

    # Step voltages: V = (E_low - E_high) / (z*delta_n) + e_metal / z
    # Derivation: reaction M_lo + dn*M(metal) -> M_hi, dE = E_hi - E_lo - dn*mu_M,
    # V = -dE/(z*dn) = (E_lo - E_hi)/(z*dn) + mu_M/z. With mu_M negative, this
    # lowers V. z = 1 for Li/Na/K (reduces to the original formula), 2 for Mg.
    results = []
    for i in range(len(recorded) - 1):
        hi = recorded[i]
        lo = recorded[i + 1]
        dn = _step_n(hi) - _step_n(lo)
        if dn <= 0:
            continue
        v = (lo["energy"] - hi["energy"]) / (dn * z) + e_li_metal / z
        x = _step_n(lo) / n_li_total
        results.append({
            "x": x,
            "n_li": _step_n(lo),
            "voltage": v,
            "n_li_from": _step_n(hi),
        })

    # Convex hull on formation energies
    # dE(x) = E(x) - x*E(fully_lith) - (1-x)*E(fully_delith)
    # where x = n_ion / n_ion_total
    e_full = None
    e_empty = None
    for s in recorded:
        if _step_n(s) == n_li_total:
            e_full = s["energy"]
        if _step_n(s) == 0:
            e_empty = s["energy"]

    hull_voltages = []
    if e_full is not None and e_empty is not None:
        # Compute formation energies
        points = []  # (x, E, formation_energy)
        for s in recorded:
            x = _step_n(s) / n_li_total
            fe = s["energy"] - x * e_full - (1 - x) * e_empty
            points.append((x, _step_n(s), s["energy"], fe))
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
            v = (e_lo - e_hi) / (dn * z) + e_li_metal / z
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

            # Host formula (excludes the working ion) for x-axis labels:
            # "x in Li_x<host>". Keyed on the first suffix token so that
            # JVASP-2017-LCO-PBE still resolves to the LCO host.
            _HOST_FORMULA = {
                "LFP": "FePO$_4$",
                "LMP": "MnPO$_4$",
                "LMO": "Mn$_2$O$_4$",
                "NMC": "Mn$_3$Co$_2$Ni$_3$O$_{16}$",
                "LCO": "CoO$_2$",
                "NCO": "CoO$_2$",
                "NMO": "MnO$_2$",
                "MMO": "Mn$_2$O$_4$",
            }
            _host = _HOST_FORMULA.get(_abbrev_parts[0] if _abbrev_parts else "", "MO")
            _xlabel = f"x in {ion}$_x${_host}"
            # Plot file id: <jid> for legacy dirs, <jid>-<tag> for tagged campaigns
            _plot_id = f"{data['jid']}-{data['tag']}" if data.get("tag") else data["jid"]

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
            ax.set_xlabel(_xlabel)
            ax.set_ylabel("Voltage (V)")
            ax.set_title(f"Discharge curve — {data['jid']}{_title_suffix}")
            ax.legend()
            ax.set_xlim(-0.05, 1.05)
            fig.tight_layout()
            analysis_dir = Path(supercell_dir).parents[2] / "analysis"
            if analysis_dir.exists():
                fig.savefig(str(analysis_dir / f"discharge_curve_{_plot_id}.png"), dpi=150)
                print(f"Saved discharge_curve_{_plot_id}.png in {analysis_dir}")
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
            ax.set_xlabel(_xlabel)
            ax.set_ylabel("Voltage (V)")
            ax.set_title(f"Voltage curve — {data['jid']}{_title_suffix}")
            ax.legend()
            ax.set_xlim(-0.05, 1.05)
            fig.tight_layout()
            if analysis_dir.exists():
                fig.savefig(str(analysis_dir / f"voltage_curve_{_plot_id}.png"), dpi=150)
                print(f"Saved voltage_curve_{_plot_id}.png in {analysis_dir}")
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
    p_init.add_argument("--functional", choices=["auto", "pbe", "optb88vdw", "optpbevdw"],
                        default="auto",
                        help="DFT functional (default: auto-detect)")
    p_init.add_argument("--ion", default="Li", choices=sorted(_ION_CHARGE),
                        help="Working ion to remove (default: Li). Later commands "
                             "read the ion from energies.json.")
    p_init.add_argument("--suffix", default=None,
                        help="Directory suffix: write <JID>-<suffix>/ (e.g. LCO-PBE)")
    p_init.add_argument("--db-cache", default=None,
                        help="Local JARVIS dft_3d snapshot (.pkl or .json) instead of "
                             "the figshare download; also env BATTERYMAT_DFT3D_CACHE")
    p_init.add_argument("--e-ion-metal", type=float, default=None,
                        help="Metal reference energy (eV/atom) to store in energies.json")

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
    p_volt.add_argument("--e-ion-metal", type=float, default=None,
                        help="Working-ion metal energy (eV/atom); same as --e-li-metal, "
                             "for Na/Mg chains whose reference is computed later")

    args = parser.parse_args()

    if args.command == "init":
        for jid in args.jids:
            generate_sequential_init(
                jid, output_dir=args.output_dir,
                min_length=args.min_length,
                max_atoms=args.max_atoms,
                functional=args.functional,
                ion=args.ion,
                suffix=args.suffix,
                db_cache=args.db_cache,
                e_ion_metal=args.e_ion_metal,
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
            e_ion_metal=args.e_ion_metal,
        )
        ion_label = _ion_info(_load_energies(sup_dir))[0].lower()
        for r in results:
            print(f"  x={r['x']:.3f}  n_{ion_label}={r['n_li']}  V={r['voltage']:.4f}")

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    _cli()
