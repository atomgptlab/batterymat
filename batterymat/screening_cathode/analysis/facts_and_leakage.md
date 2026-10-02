# Dataset facts + five-cathode leakage check

- Label pool (screening_anode/Li.csv): **7193** unique JVASP ids (paper states 7,610 structures in the training pool - reconcile the delta: attempted vs completed label generation).
- Filtered pool (average_voltage/Li_min.csv): 7193 ids.

## Five benchmark cathodes vs the label pool (exact membership)
| JVASP id | Material | In training/label pool? |
|---|---|---|
| JVASP-42723 | LFP (LiFePO4) | YES - in label pool |
| JVASP-116897 | LMP (LiMnPO4) | YES - in label pool |
| JVASP-141792 | LMO (LiMn2O4) | YES - in label pool |
| JVASP-144791 | NMC (Li[Ni,Mn,Co]O2) | YES - in label pool |
| JVASP-2017 | LCO (LiCoO2) | YES - in label pool |

## Structural near-duplicates of the five cathodes inside the label pool
| Benchmark cathode | Same reduced formula in pool (n) | Near-duplicates (StructureMatcher) |
|---|---|---|
| LFP (LiFePO4) (JVASP-42723) | 19 | none |
| LMP (LiMnPO4) (JVASP-116897) | 38 | JVASP-43659 |
| LMO (LiMn2O4) (JVASP-141792) | 6 | JVASP-34813 |
| NMC (Li[Ni,Mn,Co]O2) (JVASP-144791) | 1 | JVASP-40649 |
| LCO (LiCoO2) (JVASP-2017) | 2 | none |

Interpretation for the rebuttal: exact membership answers whether the five cathodes were in-sample for the *label pool*; whether they were in the surrogate's train vs held-out TEST split additionally requires the split record from the training run (on the cluster - see README).
