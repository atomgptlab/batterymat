# Ranking / enrichment analysis (existing predictions only)

Cathode-window pool: **2944** unique materials (max V <= 5.5, min V > 1.0, dedup by jid)  [paper reports 2890; delta likely an additional completion constraint or pool-snapshot difference - RECONCILE BEFORE QUOTING]

## Where the five DFT-validated cathodes fall in the tier-2 ranking
| Cathode | Tier-2 rank (of pool) | Percentile | Tier-2 avg V | Exp avg V |
|---|---|---|---|---|
| LFP (JVASP-42723) | 1613 | 45.2 | 3.49 | 3.45 |
| LMP (JVASP-116897) | 2021 | 31.4 | 3.17 | 4.1 |
| LMO (JVASP-141792) | 728 | 75.3 | 4.07 | 4.15 |
| NMC (JVASP-144791) | 546 | 81.5 | 4.18 | 3.8 |
| LCO (JVASP-2017) | 1110 | 62.3 | 3.84 | 3.9 |

Spearman(tier-2 avg V, experimental avg V) over the 5 in-pool anchors: **0.00** (N is small - state this honestly; it motivates the prospective run).

_Tier-1 single-shot predictions not merged (no --tier1-csv). Pull the surrogate's pool predictions from the training cluster to add Spearman(tier-1, tier-2) over the full pool._

