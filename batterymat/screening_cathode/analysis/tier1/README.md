# Tier-1 ALIGNN average-voltage model: revision analyses

Everything in this directory uses the **archived tier-1 model** (ALIGNN
`alignn_atomwise`, 4 ALIGNN + 4 GCN layers, 256 hidden, cutoff 8 Å, 12
neighbours, CGCNN atom features, seed 123, 80/10/10 split with
`keep_data_order`, 250 epochs) trained on the ALIGNN-FF average Li
intercalation voltage of 7,610 JARVIS-DFT structures:

```
/Users/jaelee/Desktop/archive/hpc_runs/bmat_result/bmat_inverse/Li_250/
  voltage/best_model.pt, config.json, ids_train_val_test.json,
  Train/Val/Test_results.json, history_train.json, history_val.json, mad
  id_prop.json          structures + targets (7,610)
  summary.csv           ALIGNN-FF labels the model was trained on (7,623 rows)
  jarvis_dft3d_cache.pkl  JARVIS dft_3d snapshot used for training structures
```

Reviewer items covered: R2.M2 (tier-1 vs tier-2 over the pool), R2.M7
(tier-1 as the front of the funnel), R2.M8 (screening confusion matrix,
residual bounds), R1.M3 (unrelaxed-input robustness), R1.m4 (five benchmark
cathodes), R1.m5 (loss curves), plus the regenerated parity plot.

All numbers below are copied from the JSON/CSV files produced by the scripts;
nothing is estimated. Python: `/Users/jaelee/miniconda3/envs/alignn/bin/python`
(torch 2.2.1, dgl 1.1.1, alignn 2025.4.1, CPU). Run every script from this
directory.

## Files

| File | What it is |
|---|---|
| `predict.py` | Inference harness (`Tier1Predictor`), loads `best_model.pt` with the archived config and the identical graph construction (k-nearest, cutoff 8.0, max_neighbors 12, use_canonize, line graph, cutoff_extra 3.0). `--sanity` reproduces `Test_results.json`. |
| `test_reproduction.csv` | Archived vs recomputed prediction for all 761 test structures. |
| `test_metrics.py` / `test_metrics.json` | Held-out test metrics, 3.0–4.5 V confusion matrix, residual percentiles, bootstrap CI. |
| `residuals.png` | Residual histogram (2.5/97.5 percentiles dashed) + test parity with the screening window shaded. |
| `alignn_parity.py` / `alignn_parity.png` / `.pdf` | Regenerated manuscript parity plot (761 test structures, dashed y=x, dotted test-set mean baseline, MAE/R²/MAD in panel, 300 dpi, 6.2×6.2 in). |
| `loss_curves.py` / `loss_curves.png` / `loss_curves.json` | Train vs validation MSE per epoch from `history_train.json` / `history_val.json`. |
| `run_pool.py` / `tier1_pool_predictions.csv` | Tier-1 prediction for every `Li_` JID in `average_voltage/Li_min.csv` (7,193), with `v_tier2` (Li_min.csv `avg_voltage`), `v_label` (summary.csv `avg_voltage_V`), split partition and tier-2 `max_voltage`. |
| `pool_stats.py` / `pool_stats.json` / `pool_parity.png` | Agreement statistics between tier-1, tier-2 (Li_min.csv) and the training labels (summary.csv). |
| `front_ranking.py` / `tier1_front_ranking.csv` / `front_ranking_comparison.json` | Tier-1-as-front ranking with the live `screen_cathode.py` filters and score; overlap with `cathode_candidates_ranked.csv`. |
| `front_ranking_plot.py` / `front_ranking.png` | Rank-vs-rank and voltage-vs-voltage plots for the two fronts. |
| `benchmarks.py` / `benchmarks_tier1.csv` | Tier-1 prediction, partition membership, labels, DFT and experimental averages for the five DFT benchmark cathodes. |
| `robustness.py` / `robustness.csv` / `robustness_per_structure.csv` | Unrelaxed-input test on a random 200-structure subset of the held-out test set. |

## 1. Harness sanity check (mandatory)

`python predict.py --sanity` rebuilds every test structure from `id_prop.json`
and compares with the archived `Test_results.json` (which was written by the
end-of-training model; the best validation epoch is epoch 250 = the last one,
so `best_model.pt` is that model):

* 761/761 test predictions reproduced; **max |Δ| = 1.7e-6 V**, mean |Δ| =
  2.2e-7 V, 0 structures above the 1e-3 V tolerance (`test_reproduction.csv`).
* Throughput on CPU: ~0.10–0.14 s per structure (7,193-structure pool took
  940 s).

Held-out test metrics recomputed from `Test_results.json` (n = 761):

| Metric | Value |
|---|---|
| MAE | **0.1672 V** |
| RMSE | 0.3198 V |
| R² | **0.9431** |
| Pearson r / Spearman ρ | 0.9712 / 0.9672 |
| MAD of the test-set labels (mean-predictor baseline on the test set) | **1.0325 V** |
| MAD of all 7,610 labels (value in the archived `mad` file) | **1.2444 V** |
| MAE / MAD(test) | 0.162 |

The two MAD values differ because the archived `mad` file is computed on the
whole dataset, while the mean-predictor baseline on the test split alone is
1.03 V. Both are reported; the parity plot prints the test-set value.

## 2. Tier-1 over the full Li pool (R2.M2)

Pool: 7,193 unique `Li_` JIDs from `Li_min.csv`; all 7,193 are in the JARVIS
cache and all were predicted (`tier1_pool_predictions.csv`).

**Important:** every pool JID is in the model's dataset (train 5,756 / val 699
/ test 738; none unseen), because the tier-1 training set was built from the
same JARVIS intercalation set. "Not in the training partition" below means
val + test (1,437 JIDs).

Also important: `Li_min.csv avg_voltage` (tier-2 as used in the paper) and
`summary.csv avg_voltage_V` (the labels the tier-1 model was trained on) are
**two different ALIGNN-FF runs with a systematic offset**: label − Li_min =
−0.628 V mean (median −0.593 V, std 0.805 V); only 7.4 % of JIDs agree within
0.1 V. The tier-1 model therefore reproduces its own labels closely but sits
~0.6 V below the Li_min.csv values.

| Pair (n) | Spearman ρ | Pearson r | MAE (V) | mean diff (V) | 3.0–4.5 V window agreement |
|---|---|---|---|---|---|
| tier-1 vs tier-2 (Li_min), all 7,193 | **0.868** | 0.877 | **0.776** | −0.633 | 73.0 % |
| tier-1 vs tier-2, val+test only (1,437) | **0.842** | 0.851 | 0.745 | −0.625 | 74.0 % |
| tier-1 vs tier-2, test only (738) | 0.823 | 0.809 | 0.799 | −0.701 | 73.4 % |
| tier-1 vs tier-2, restricted to tier-2 window (3,034) | 0.538 | 0.486 | 0.717 | −0.657 | 67.0 % |
| tier-1 vs training label (summary.csv), all | 0.994 | 0.994 | 0.064 | −0.005 | 97.6 % |
| tier-1 vs training label, val+test (1,437) | 0.964 | 0.964 | 0.186 | −0.003 | 93.7 % |
| tier-1 vs training label, test (738) | 0.966 | 0.971 | 0.165 | −0.017 | 94.3 % |
| training label vs tier-2 (Li_min), all | 0.868 | 0.878 | 0.770 | −0.628 | 73.3 % |
| training label vs tier-2, val+test | 0.833 | 0.847 | 0.749 | −0.621 | 74.1 % |

(mean diff = first − second; full detail in `pool_stats.json`; figure
`pool_parity.png`.) The tier-1/tier-2 disagreement is almost entirely the
label-set offset: label-vs-tier-2 (0.868 / 0.770 V) is indistinguishable from
tier-1-vs-tier-2 (0.868 / 0.776 V).

## 3. Tier-1 as the front of the funnel (R2.M7)

`front_ranking.py` imports the live `screen_cathode.py` (filters: 3.0 ≤ V ≤
4.5, q_grav > 100 mAh/g, ehull ≤ 0.05 eV, max_voltage ≤ 5.5 V, and the
redox-metal requirement added on 2026-09-24; score =
(norm V + norm q_grav − norm ehull)/3) and replaces `avg_voltage` by the tier-1
prediction. q_grav/ehull come from `cathode_candidates_ranked.csv` for its
682 JIDs and from the JARVIS cache otherwise (recomputed values agree with the
CSV to 6e-14 mAh/g and 0 eV). The tier-2 `max_voltage` is kept. Rebuilding
the tier-2 ranking from the same pool frame reproduces
`cathode_candidates_ranked.csv` exactly (682 JIDs, same order).

| Quantity | Value |
|---|---|
| Survivors, tier-1 front | **571** (train 413 / val 73 / test 85) |
| Survivors, tier-2 front (`cathode_candidates_ranked.csv`) | 682 |
| Common survivors | 507 (Jaccard 0.680; 64 tier-1-only, 175 tier-2-only) |
| Spearman ρ of ranks over the 507 common survivors | **0.852** |
| Top-12 overlap | **6 / 12** (JVASP-117295, 122347, 118606, 31181, 49957, 11735) |
| Top-50 overlap | **24 / 50** |

Tier-1 top-12: JVASP-117295 Li3FeO3 (tier-2 rank 1), 143995 Li2VF6 (not a
tier-2 survivor, V_tier2 = 4.71), 122347 KLi4FeO5 (3), 118606 Li2NiO2 (9),
31181 Li3RuO4 (5), 19051 Li6WN4 (V_tier2 = 2.45), 49957 Li2RhO3 (10), 11735
Li2MnO3 (2), 113201 Li2Co3NiO8 (V_tier2 = 4.56), 34836 LiCo2O4 (4.72), 119348
LiCoNiO4 (18), 45853 Li2VF6 (4.72). Because tier-1 runs ~0.6 V low, it admits
64 materials that tier-2 puts above 4.5 V and drops 175 that tier-2 puts in
3.0–3.3 V (see `front_ranking.png`, right panel).

Where named materials land (rank = None means filtered out; reason given):

| JID | Formula | Partition | V_tier1 | V_tier2 | rank tier-1 | rank tier-2 | Filter status |
|---|---|---|---|---|---|---|---|
| JVASP-2017 (LCO) | LiCoO2 | train | 3.418 | 3.842 | **231** | 263 | pass / pass |
| JVASP-117419 | LiCr2P2O8 | train | 3.199 | 4.420 | – | – | q_grav = 89.1 (both) |
| JVASP-141543 | Rb2LiFeF6 | val | 4.190 | 4.079 | – | – | q_grav = 77.1 (both) |
| JVASP-117295 | Li3FeO3 | train | 4.436 | 4.394 | **1** | 1 | pass / pass |
| JVASP-116849 | LiV2F7 | train | 3.582 | 4.187 | **108** | 86 | pass / pass |
| JVASP-154749 | Rb2LiFeF6 | test | 4.203 | 4.253 | – | – | q_grav = 77.1 (both) |
| JVASP-96563 | Li2MgMn3O8 | train | 3.663 | 4.230 | **50** | 46 | pass / pass |
| JVASP-95533 | LiCoH24C8N8O12 | train | 1.939 | 4.105 | – | – | tier-1: V = 1.94 and q_grav = 54.7; tier-2: q_grav |
| JVASP-97428 | LiH4SNO4 | train | 3.785 | 4.011 | – | – | no redox metal (both) |
| JVASP-80802 (control) | LiY2Ga | train | 0.611 | 1.043 | – | – | V, ehull = 0.313, no redox metal |
| JVASP-77457 (control) | LiCa2Ag | train | 0.964 | 1.063 | – | – | V, ehull = 0.121, no redox metal |
| JVASP-85076 (control) | LiCaTl2 | train | 0.294 | 1.105 | – | – | V, q_grav = 58.8, no redox metal |
| JVASP-81829 (control) | LiMnIr2 | train | 0.768 | 1.120 | – | – | V, q_grav = 60.1, ehull = 0.712 |

Note that 5 of the 8 "top" prospective JIDs (117419, 141543, 154749, 95533,
97428) are no longer survivors of the *tier-2* ranking either under the
current filters (q_grav per formula unit > 100 mAh/g and the redox-metal
requirement); they were selected from an earlier cell-normalised ranking.

## 4. Screening confusion matrix and residual bounds (R2.M8)

Held-out test set (n = 761), in-window = 3.0 ≤ V ≤ 4.5 by label vs by tier-1
prediction (`test_metrics.json`):

| | pred in | pred out |
|---|---|---|
| **label in** (385) | TP = **366** | FN = **19** |
| **label out** (376) | FP = **24** | TN = **352** |

Precision = **0.938**, recall = **0.951**, F1 = 0.945, accuracy = 0.943.
FN/FP JIDs are listed in the JSON.

Residuals (pred − label): mean −0.012 V, std 0.320 V, median −0.013 V.
|residual| percentiles: p50 0.058, p90 0.453, **p95 0.677**, **p99 1.416**,
max 2.023 V. Signed 2.5/97.5 % bounds [−0.707, +0.664] V; 0.5/99.5 % bounds
[−1.497, +1.225] V. 79.5 % of test structures are within 0.25 V and 90.9 %
within 0.5 V. Bootstrap (10,000 resamples, seed 0) 95 % CI on the test MAE:
**[0.148, 0.187] V**. Figure: `residuals.png`.

## 5. Loss curves (R1.m5)

`history_train.json` / `history_val.json` hold one entry per epoch
`[total_loss, graphwise_loss, 0, 0, 0, 0]`, where the loss is the MSE summed
over the batches of that epoch (190 train batches, 23 val batches of 32;
`drop_last`). `loss_curves.png` shows both the raw sums (left) and the
per-batch means (right). Numbers (`loss_curves.json`):

| Epoch | train (summed) | val (summed) |
|---|---|---|
| 1 | 133.96 | 13.08 |
| 50 | 20.00 | 6.97 |
| 100 | 11.87 | 6.31 |
| 150 | 9.59 | 5.87 |
| 200 | 8.20 | 6.00 |
| 250 | 7.18 | 5.75 |

Best validation loss at epoch 250 (the final epoch). Final per-batch RMSE:
train 0.194 V, val 0.500 V (val/train MSE ratio 6.6). From the archived
per-structure outputs, train MAE = 0.038 V (n = 6,080) vs val MAE = 0.250 V
(n = 736) vs test MAE = 0.167 V: the model fits the training set much more
tightly than held-out data, but the validation curve is still flat/decreasing
at epoch 250 (no upturn), so 250 epochs is not past the overfitting point by
the validation criterion.

## 6. Five benchmark cathodes (R1.m4)

`benchmarks_tier1.csv` (V in volts; DFT and experiment from the main
analysis):

| JID | Material | Partition | tier-1 | training label | tier-2 (Li_min) | DFT avg | Exp |
|---|---|---|---|---|---|---|---|
| JVASP-42723 | LiFePO4 (LFP) | **train** | 3.246 | 3.267 | 3.490 | 3.60 | 3.45 |
| JVASP-116897 | LiMnPO4 (LMP) | **train** | 3.095 | 3.099 | 3.166 | 3.91 | ~4.1 |
| JVASP-141792 | LiMn2O4 (LMO) | **train** | 3.597 | 3.598 | 4.065 | 4.08 | 4.1 |
| JVASP-144791 | Li4Mn3Co2Ni3O16 (NMC) | **val** | 3.700 | 3.689 | 4.178 | 4.40 | ~3.7 |
| JVASP-2017 | LiCoO2 (LCO) | **train** | 3.418 | 3.449 | 3.842 | 4.18 | 3.9–4.2 |

Four of the five benchmarks were in the training partition and NMC was in the
validation partition; none was held out. The archived per-structure
prediction is recoverable only for val/test (NMC: 3.7004 V archived vs 3.7004
V recomputed); `Train_results.json` comes from a shuffled loader.

## 7. Unrelaxed-input robustness (R1.M3)

Random 200-structure subset (numpy seed 0) of the 761 held-out test
structures. Conditions: (a) the JARVIS (OptB88vdW-relaxed) structure; (b)
every atom displaced by isotropic Gaussian noise with σ = 0.05 Å and 0.10 Å
per Cartesian component, 3 independent draws each; (c) ALIGNN-FF relaxation.

| Condition | n | MAE vs label (V) | mean |shift| vs relaxed-input prediction (V) | max |shift| (V) |
|---|---|---|---|---|
| relaxed (JARVIS) | 200 | 0.161 | 0.000 | 0.000 |
| noise sigma=0.05 A, draw 0 | 200 | 0.170 | 0.048 | 0.279 |
| noise sigma=0.05 A, draw 1 | 200 | 0.173 | 0.050 | 0.379 |
| noise sigma=0.05 A, draw 2 | 200 | 0.173 | 0.051 | 0.343 |
| noise sigma=0.05 A, pooled 3 draws | 600 | 0.172 | 0.050 | 0.379 |
| noise sigma=0.10 A, draw 0 | 200 | 0.196 | 0.086 | 0.526 |
| noise sigma=0.10 A, draw 1 | 200 | 0.196 | 0.085 | 0.396 |
| noise sigma=0.10 A, draw 2 | 200 | 0.198 | 0.098 | 0.442 |
| noise sigma=0.10 A, pooled 3 draws | 600 | 0.197 | 0.089 | 0.526 |
| ALIGNN-FF relaxed | 0 | – | – | – |
Subset relaxed-input MAE = 0.161 V (n = 200). Pooled rows are the 3 draws stacked (n = 600). Per-structure predictions are in `robustness_per_structure.csv`.

**Condition (c) was not completed.** With `alignn.ff` in this environment
(alignn 2025.4.1, torch 2.2.1, CPU) the `alignnff_wt01` checkpoint returns
forces of 12–100 eV/Å on DFT-relaxed cells (LiCoO2 24.3 eV/Å, LiFePO4 20.8
eV/Å, Li7Mn5O12 42.1 eV/Å; Si 0.5 eV/Å), and FIRE does not converge (LiCoO2:
fmax 10.7 eV/Å after 200 positions-only steps, 93 s; with `ExpCellFilter` the
first 10 steps of the first structure took 100 s at fmax 20–85 eV/Å). Setting
`add_reverse_forces=False` only halves the forces. The current default
ALIGNN-FF checkpoint (`v12.2.2024_dft_3d_307k`) could not be obtained: the
figshare download returns a 0-byte zip (WAF block). A relaxation with these
forces would not be an ALIGNN-FF relaxation, so (c) is reported as skipped
rather than with meaningless numbers. If the `v12.2.2024_dft_3d_307k` weights
are placed in `/Users/jaelee/software/alignn/alignn/ff/v12.2.2024_dft_3d_307k/`
(`best_model.pt` + `config.json`), set `FF_PATH` in `robustness.py` and run it
without `--skip-ff`.

## Regeneration

```bash
cd /Users/jaelee/Desktop/work/batterymat_jae/batterymat_jae/screening_cathode/analysis/tier1
PY=/Users/jaelee/miniconda3/envs/alignn/bin/python
$PY predict.py --sanity                # 761-structure reproduction check (~90 s)
$PY test_metrics.py                    # test_metrics.json, residuals.png, test_reproduction.csv
$PY alignn_parity.py                   # alignn_parity.png/.pdf
$PY loss_curves.py                     # loss_curves.png/.json
$PY run_pool.py                        # tier1_pool_predictions.csv (~15 min CPU)
$PY pool_stats.py                      # pool_stats.json, pool_parity.png
$PY front_ranking.py                   # tier1_front_ranking.csv, front_ranking_comparison.json
$PY front_ranking_plot.py              # front_ranking.png
$PY benchmarks.py                      # benchmarks_tier1.csv
$PY robustness.py --skip-ff            # robustness.csv, robustness_per_structure.csv (~5-20 min)
```

`front_ranking.py` must be rerun whenever `cathode_candidates_ranked.csv` or
the filters in `screen_cathode.py` change. All scripts print their headline
numbers; long jobs are best run with `nohup ... &`.
