"""Cathode candidate screening from Li_min.csv + JARVIS-DFT metadata."""
import sys
import warnings
from typing import Optional
import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# Filtering thresholds (spec Section 3)
# ---------------------------------------------------------------------------
FILTERS = {
    "avg_voltage_min": 3.0,
    "avg_voltage_max": 4.5,
    "max_grav_cap_min": 20.0,
    "ehull_max": 0.05,
    "max_voltage_max": 5.5,
}


def extract_jid(name: str) -> str:
    """Extract JARVIS JID from a Li_min.csv name field.

    Example: 'Li_JVASP-47127_Li3SiNi3O8.json' -> 'JVASP-47127'
    """
    try:
        return name.split("_")[1]
    except IndexError as exc:
        raise ValueError(
            f"Cannot extract JID from name {name!r}. "
            "Expected format: '<Ion>_<JID>_<formula>.json' "
            "(e.g. 'Li_JVASP-47127_Li3SiNi3O8.json')."
        ) from exc


def load_and_merge(
    li_min_df: pd.DataFrame,
    dft3d_df: Optional[pd.DataFrame],
) -> pd.DataFrame:
    """Merge Li_min.csv DataFrame with JARVIS-DFT metadata on JID.

    Args:
        li_min_df: DataFrame loaded from Li_min.csv.
        dft3d_df: DataFrame from jarvis.db.figshare.data('dft_3d').
                  Pass None to trigger a RuntimeError (offline guard).

    Returns:
        Merged DataFrame with ehull, optb88vdw_bandgap, formula, atoms added.
    """
    if dft3d_df is None:
        raise RuntimeError(
            "JARVIS-DFT dataset is None. Fetch it with "
            "jarvis.db.figshare.data('dft_3d') or provide a cached copy."
        )

    li_min_df = li_min_df.copy()

    if "name" not in li_min_df.columns:
        raise ValueError(
            "li_min_df is missing the required 'name' column. "
            "Expected a DataFrame loaded from Li_min.csv with a 'name' column "
            "containing entries like 'Li_JVASP-47127_Li3SiNi3O8.json'."
        )

    li_min_df["jid"] = li_min_df["name"].apply(extract_jid)

    merged = pd.merge(li_min_df, dft3d_df[["jid", "ehull", "optb88vdw_bandgap", "formula", "atoms"]], on="jid", how="inner")

    missing_jids = set(li_min_df["jid"]) - set(merged["jid"])
    if missing_jids:
        print(
            f"WARNING: {len(missing_jids)} JIDs from Li_min.csv not found in "
            f"JARVIS-DFT dataset and were dropped: {sorted(missing_jids)}",
            file=sys.stdout,
        )

    return merged.reset_index(drop=True)


def filter_cathode_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """Apply cathode-specific voltage, capacity, and stability filters.

    Filters (spec Section 3):
      - avg_voltage: 3.0 – 4.5 V (inclusive)
      - max_grav_cap: > 20 mAh/g (strict)
      - ehull: <= 0.05 eV
      - max_voltage: <= 5.5 V
    """
    mask = (
        (df["avg_voltage"] >= FILTERS["avg_voltage_min"])
        & (df["avg_voltage"] <= FILTERS["avg_voltage_max"])
        & (df["max_grav_cap"] > FILTERS["max_grav_cap_min"])
        & (df["ehull"] <= FILTERS["ehull_max"])
        & (df["max_voltage"] <= FILTERS["max_voltage_max"])
    )
    return df[mask].reset_index(drop=True)


def _safe_norm(series: pd.Series) -> pd.Series:
    """Min-max normalize; returns 0.5 for all if constant."""
    lo, hi = series.min(), series.max()
    if hi == lo:
        return pd.Series(0.5, index=series.index)
    return (series - lo) / (hi - lo)


def rank_candidates(df: pd.DataFrame) -> pd.DataFrame:
    """Add composite score and sort candidates descending.

    Score = (1/3)*norm(avg_voltage) + (1/3)*norm(max_grav_cap) - (1/3)*norm(ehull)

    If fewer than 2 candidates, returns unranked with a warning.
    """
    if len(df) < 2:
        print(
            "WARNING: Fewer than 2 candidates after filtering — "
            "composite score not computed.",
            file=sys.stdout,
        )
        return df.copy()

    ranked = df.copy()
    ranked["score"] = (
        (1 / 3) * _safe_norm(ranked["avg_voltage"])
        + (1 / 3) * _safe_norm(ranked["max_grav_cap"])
        - (1 / 3) * _safe_norm(ranked["ehull"])
    )
    return ranked.sort_values("score", ascending=False).reset_index(drop=True)


# Default paths relative to repo root
_REPO_ROOT = __file__  # batterymat/screening_cathode/screen_cathode.py
import os as _os
_DEFAULT_LI_MIN = _os.path.join(
    _os.path.dirname(_os.path.dirname(_os.path.abspath(_REPO_ROOT))),
    "average_voltage", "Li_min.csv",
)
_DEFAULT_OUTPUT = _os.path.join(
    _os.path.dirname(_os.path.abspath(_REPO_ROOT)),
    "cathode_candidates_ranked.csv",
)


def run_screening(
    li_min_df: Optional[pd.DataFrame] = None,
    dft3d_df: Optional[pd.DataFrame] = None,
    output_path: Optional[str] = None,
) -> pd.DataFrame:
    """Full screening pipeline: load → merge → filter → rank → export.

    Args:
        li_min_df: Pre-loaded Li_min.csv DataFrame (or None to read from disk).
        dft3d_df: Pre-loaded JARVIS-DFT DataFrame (or None to fetch from figshare).
        output_path: Where to write the ranked CSV. Defaults to
                     batterymat/screening_cathode/cathode_candidates_ranked.csv.

    Returns:
        Ranked candidates DataFrame.
    """
    if li_min_df is None:
        li_min_df = pd.read_csv(_DEFAULT_LI_MIN)

    if dft3d_df is None:
        try:
            from jarvis.db.figshare import data as jarvis_data
            dft3d_df = pd.DataFrame(jarvis_data("dft_3d"))
        except Exception as e:
            raise RuntimeError(
                f"Failed to fetch JARVIS-DFT dataset: {e}. "
                "Check your internet connection or provide dft3d_df directly."
            ) from e

    merged = load_and_merge(li_min_df, dft3d_df)
    filtered = filter_cathode_candidates(merged)
    ranked = rank_candidates(filtered)

    # Ensure 'score' column is always present in the output CSV
    if "score" not in ranked.columns:
        ranked = ranked.copy()
        ranked["score"] = float("nan")

    save_path = output_path or _DEFAULT_OUTPUT
    ranked.to_csv(save_path, index=False)
    print(f"Saved {len(ranked)} cathode candidates to {save_path}")

    return ranked


if __name__ == "__main__":
    candidates = run_screening()
    print(candidates[["jid", "formula", "avg_voltage", "max_grav_cap", "ehull", "score"]].to_string(index=False))
