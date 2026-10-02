"""Task 7 (R1.M3): unrelaxed-input robustness of the tier-1 model.

Random 200-structure subset (seed 0) of the 761 held-out test structures.
Conditions:
  relaxed          JARVIS (OptB88vdW-relaxed) structure as used for training
  noise0.05_dN     every atom displaced by isotropic Gaussian noise, sigma 0.05 A
                   per Cartesian component (3 independent draws, N = 0,1,2)
  noise0.10_dN     same with sigma 0.10 A
  alignnff_relaxed structure re-relaxed with ALIGNN-FF (alignnff_wt01, FIRE,
                   fmax 0.05 eV/A, <=200 steps, cell + positions)

Run with ``--skip-ff`` to omit condition (c). NOTE (2026-09-24): in this
environment (alignn 2025.4.1, torch 2.2.1, CPU) the alignnff_wt01 checkpoint
returns forces of 12-100 eV/A on JARVIS-relaxed structures (LiCoO2: 24 eV/A,
LiFePO4: 21 eV/A, Si: 0.5 eV/A) and FIRE does not converge (LiCoO2: fmax 10.7
eV/A after 200 steps, 93 s); the current default ALIGNN-FF checkpoint
(v12.2.2024_dft_3d_307k) cannot be downloaded (figshare WAF returns a 0-byte
zip). Condition (c) was therefore skipped; see README.md.
Outputs: robustness.csv (summary), robustness_per_structure.csv (all preds)
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from predict import ARCHIVE, MODEL_DIR, Tier1Predictor  # noqa: E402

N_SUB, SEED = 200, 0
SIGMAS = (0.05, 0.10)
N_DRAWS = 3
FF_PATH = "/Users/jaelee/software/alignn/alignn/ff/alignnff_wt01"
FF_TRIAL = 20            # time the first 20 relaxations
FF_MAX_MIN = 150         # skip the rest if the projected total exceeds this


def perturb(atoms_dict, sigma, rng):
    from jarvis.core.atoms import Atoms
    a = Atoms.from_dict(atoms_dict)
    cart = a.cart_coords + rng.normal(0.0, sigma, size=(a.num_atoms, 3))
    return Atoms(lattice_mat=a.lattice_mat, elements=a.elements, coords=cart, cartesian=True)


def ff_relax(atoms_dict):
    from jarvis.core.atoms import Atoms
    from alignn.ff.ff import ForceField
    a = Atoms.from_dict(atoms_dict)
    ff = ForceField(jarvis_atoms=a, model_path=FF_PATH, logfile=None)
    relaxed, energy, forces = ff.optimize_atoms(optimizer="FIRE", logfile=None, trajectory=None,
                                                steps=200, fmax=0.05, optimize_lattice=True, interval=None)
    return relaxed, float(np.abs(forces).max())


def main(skip_ff=False):
    test = json.load(open(os.path.join(MODEL_DIR, "Test_results.json")))
    idprop = {d["jid"]: d for d in json.load(open(os.path.join(ARCHIVE, "id_prop.json")))}
    rng = np.random.default_rng(SEED)
    sub = [test[i] for i in sorted(rng.choice(len(test), N_SUB, replace=False))]
    p = Tier1Predictor()
    rows = []
    t0 = time.time()
    for k, r in enumerate(sub):
        j = r["id"]; ad = idprop[j]["atoms"]
        row = {"jid": j, "n_atoms": len(ad["elements"]), "label": r["target_out"][0],
               "pred_relaxed": p.predict_one(ad)}
        for s in SIGMAS:
            for d in range(N_DRAWS):
                row[f"pred_noise{s:.2f}_d{d}"] = p.predict_one(perturb(ad, s, rng))
        rows.append(row)
        if (k + 1) % 50 == 0:
            print(f"perturbation phase {k+1}/{N_SUB} {time.time()-t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(HERE, "robustness_per_structure.csv"), index=False)

    # ALIGNN-FF relaxation phase
    ff_pred, ff_fmax, ff_time, ff_dvol = [], [], [], []
    tff = time.time(); skipped_after = None
    for k, r in enumerate([] if skip_ff else sub):
        j = r["id"]; ad = idprop[j]["atoms"]
        t1 = time.time()
        try:
            from jarvis.core.atoms import Atoms
            rel, fmax = ff_relax(ad)
            v0 = Atoms.from_dict(ad).volume
            ff_pred.append(p.predict_one(rel)); ff_fmax.append(fmax); ff_dvol.append(rel.volume / v0 - 1)
        except Exception as exc:
            print(f"[warn] FF relax failed for {j}: {exc}")
            ff_pred.append(np.nan); ff_fmax.append(np.nan); ff_dvol.append(np.nan)
        ff_time.append(time.time() - t1)
        print(f"FF {k+1}/{N_SUB} {j} {ff_time[-1]:.1f}s fmax={ff_fmax[-1]:.3f} dV={ff_dvol[-1]*100:+.2f}% "
              f"V {df.loc[k,'pred_relaxed']:.3f}->{ff_pred[-1]:.3f}", flush=True)
        if k + 1 == FF_TRIAL:
            proj = np.mean(ff_time) * N_SUB / 60
            print(f"projected FF total {proj:.1f} min for {N_SUB} structures", flush=True)
            if proj > FF_MAX_MIN:
                skipped_after = FF_TRIAL
                print("too slow; stopping ALIGNN-FF phase after 20 structures", flush=True)
                break
    n_ff = len(ff_pred)
    df["pred_alignnff_relaxed"] = ff_pred + [np.nan] * (N_SUB - n_ff)
    df["ff_fmax_eV_A"] = ff_fmax + [np.nan] * (N_SUB - n_ff)
    df["ff_dvol_frac"] = ff_dvol + [np.nan] * (N_SUB - n_ff)
    df["ff_seconds"] = ff_time + [np.nan] * (N_SUB - n_ff)
    df.to_csv(os.path.join(HERE, "robustness_per_structure.csv"), index=False)

    # summary
    out = []
    def add(name, col, note=""):
        m = df[col].notna()
        out.append({"condition": name, "n": int(m.sum()),
                    "mae_vs_label_V": float(np.mean(np.abs(df.loc[m, col] - df.loc[m, "label"]))),
                    "mean_abs_shift_vs_relaxed_V": float(np.mean(np.abs(df.loc[m, col] - df.loc[m, "pred_relaxed"]))),
                    "max_abs_shift_vs_relaxed_V": float(np.max(np.abs(df.loc[m, col] - df.loc[m, "pred_relaxed"]))),
                    "mae_vs_label_on_same_subset_relaxed_V": float(np.mean(np.abs(df.loc[m, "pred_relaxed"] - df.loc[m, "label"]))),
                    "note": note})
    add("relaxed (JARVIS)", "pred_relaxed")
    for s in SIGMAS:
        cols = [f"pred_noise{s:.2f}_d{d}" for d in range(N_DRAWS)]
        for c in cols:
            add(f"noise sigma={s:.2f} A, draw {c[-1]}", c)
        # pooled over draws
        stacked = pd.DataFrame({"pred": np.concatenate([df[c].values for c in cols]),
                                "label": np.tile(df["label"].values, N_DRAWS),
                                "rel": np.tile(df["pred_relaxed"].values, N_DRAWS)})
        out.append({"condition": f"noise sigma={s:.2f} A, pooled {N_DRAWS} draws", "n": len(stacked),
                    "mae_vs_label_V": float(np.mean(np.abs(stacked["pred"] - stacked["label"]))),
                    "mean_abs_shift_vs_relaxed_V": float(np.mean(np.abs(stacked["pred"] - stacked["rel"]))),
                    "max_abs_shift_vs_relaxed_V": float(np.max(np.abs(stacked["pred"] - stacked["rel"]))),
                    "mae_vs_label_on_same_subset_relaxed_V": float(np.mean(np.abs(stacked["rel"] - stacked["label"]))),
                    "note": ""})
    if skip_ff:
        out.append({"condition": "ALIGNN-FF relaxed", "n": 0, "mae_vs_label_V": np.nan, "mean_abs_shift_vs_relaxed_V": np.nan,
                    "max_abs_shift_vs_relaxed_V": np.nan, "mae_vs_label_on_same_subset_relaxed_V": np.nan,
                    "note": "SKIPPED: alignnff_wt01 forces unphysical in this env (12-100 eV/A on relaxed cells, FIRE non-convergent); v12.2.2024 checkpoint not downloadable (figshare WAF)"})
    else:
        note = f"FIRE fmax=0.05 eV/A, <=200 steps, cell+positions; mean {np.nanmean(ff_time):.1f} s/structure"
        if skipped_after:
            note += f"; STOPPED after {skipped_after} structures (too slow on CPU)"
        add("ALIGNN-FF relaxed", "pred_alignnff_relaxed", note)
    summ = pd.DataFrame(out)
    summ.to_csv(os.path.join(HERE, "robustness.csv"), index=False)
    print(summ.to_string(index=False))
    print(f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main(skip_ff="--skip-ff" in sys.argv)
