"""Regenerate the JARVIS-DFT lithium screening distributions (SI Fig. S2).

Left: ALIGNN-FF average voltage over the 7,193-entry lithium pool.
Right: theoretical gravimetric capacity per formula unit (all Li extracted).

Usage: python make_screening_distribution.py [--cache jarvis_dft3d_cache.pkl]
Writes screening_distribution.png next to this script.
"""
import argparse
import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import screen_cathode as sc  # noqa: E402

plt.rcParams.update({"font.size": 14, "axes.labelsize": 16, "axes.titlesize": 16,
                     "xtick.labelsize": 13, "ytick.labelsize": 13})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", default=None, help="pickled JARVIS dft_3d list")
    ap.add_argument("-o", "--output", default=os.path.join(os.path.dirname(__file__), "screening_distribution.png"))
    args = ap.parse_args()
    li = pd.read_csv(sc._DEFAULT_LI_MIN)
    if args.cache:
        dft3d = pd.DataFrame(pickle.load(open(args.cache, "rb")))
    else:
        from jarvis.db.figshare import data as jarvis_data
        dft3d = pd.DataFrame(jarvis_data("dft_3d"))
    m = sc.load_and_merge(li, dft3d)
    print(f"{len(m)} Li entries; avg V median {m.avg_voltage.median():.2f} V; "
          f"Q_grav median {m.q_grav.median():.0f} mAh/g")
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    ax[0].hist(m.avg_voltage.clip(-2, 8), bins=80, color="#4C72B0")
    ax[0].axvspan(3.0, 4.5, color="0.85", zorder=0)
    ax[0].set_xlabel("ALIGNN-FF average voltage (V vs. Li/Li$^+$)")
    ax[0].set_ylabel("Structures")
    ax[0].set_title(f"(a) Average voltage, n = {len(m):,}")
    ax[1].hist(m.q_grav.clip(0, 1500), bins=80, color="#DD8452")
    ax[1].axvline(100, color="k", ls="--", lw=1.5)
    ax[1].set_xlabel("Theoretical capacity per formula unit (mAh/g)")
    ax[1].set_ylabel("Structures")
    ax[1].set_title("(b) Gravimetric capacity")
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output, dpi=300)
    print("wrote", args.output)


if __name__ == "__main__":
    main()
