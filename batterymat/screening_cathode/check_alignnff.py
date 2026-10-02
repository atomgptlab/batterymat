"""Pre-flight check of ALIGNN-FF before launching dft_prep.py chains.

`dft_prep.py next` ranks vacancies by ALIGNN-FF (alignnff_wt01) single-point
ENERGIES only. This checks that the installed stack reproduces reference
energies for DFT-relaxed LiCoO2 (JVASP-2017) in the 2x2x2 supercell: the full
cell and all 8 single-Li-vacancy cells.

Reference values (2026-09-28) agree to <0.05 meV across alignn 2026.5.20
(atomgptlab), alignn 2025.4.1 (Rockfish `bmat`, Mac) and the Mac enumeration
results. The ~59 meV spread among the 8 translation-equivalent vacancies is a
property of the model's neighbour graph, not an install fault.

Forces are printed for information only. alignn 2025.4.1 returns unphysical
forces (fmax ~24 eV/A here) with correct energies, so it is fine for vacancy
ranking but NOT for ALIGNN-FF relaxations (enumeration level b).

Usage:  python check_alignnff.py --cache /path/to/jarvis_dft3d_cache.pkl
Exit code 0 = energies reproduced, 1 = do not launch.
"""
import argparse
import pickle

import numpy as np
from jarvis.core.atoms import Atoms

from dft_prep import _alignn_energy

REF_FULL = -147.85197
REF_VAC = [-143.78147, -143.75793, -143.76217, -143.79759,
           -143.73881, -143.77621, -143.78065, -143.74913]
TOL = 1e-3  # eV per cell

ap = argparse.ArgumentParser()
ap.add_argument("--cache", required=True)
args = ap.parse_args()

d = pickle.load(open(args.cache, "rb"))
prim = Atoms.from_dict([e for e in d if e["jid"] == "JVASP-2017"][0]["atoms"])
sup = prim.make_supercell([2, 2, 2])
li = [i for i, el in enumerate(sup.elements) if el == "Li"]

e_full = float(_alignn_energy(sup))
e_vac = [float(_alignn_energy(sup.remove_site_by_index(i))) for i in li]
dev = max([abs(e_full - REF_FULL)] + [abs(a - b) for a, b in zip(e_vac, REF_VAC)])

from alignn.ff.ff import AlignnAtomwiseCalculator, wt01_path
a = prim.ase_converter()
a.calc = AlignnAtomwiseCalculator(path=wt01_path(), stress_wt=0.3)
fmax = float(np.abs(a.get_forces()).max())

print(f"max energy deviation from reference: {1000 * dev:.3f} meV/cell   (fail > {1000 * TOL:.0f})")
print(f"lowest-energy vacancy: Li {int(np.argmin(e_vac))}   (reference: Li {int(np.argmin(REF_VAC))})")
print(f"fmax on DFT-relaxed LiCoO2: {fmax:.3f} eV/A   (info only; > 5 means do not use for relaxations)")
ok = dev < TOL and int(np.argmin(e_vac)) == int(np.argmin(REF_VAC))
print("OK" if ok else "FAIL: ALIGNN-FF energies differ from reference; do not launch chains")
raise SystemExit(0 if ok else 1)
