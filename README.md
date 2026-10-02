# BatteryMat: Computational Design Pipeline for Battery Materials

BatteryMat is a computational pipeline for high-throughput screening and analysis of battery electrode materials. It combines ML force fields (ALIGNN) with DFT calculations to identify promising cathode and anode materials, predict their voltage profiles, and compute ion migration barriers. Originated at NIST by Kamal Choudhary and Michael Woodcox; developed by Jaehyung Lee at Johns Hopkins University.

## Setup

```bash
pip install -e .
```

Key dependencies: `jarvis-tools`, `alignn`, `numpy`, `scipy`. Python >= 3.8 required.

## Pipeline Stages

### Stage 1: Anode Screening (`screening_anode/`)

**Purpose:** Screen ~76,000 materials from the JARVIS-DFT database for intercalation-type anodes across 6 ion types (Li, Na, K, Ca, Mg, Al, Zn).

**How it works (`screen.py`):**
1. Loads all structures from `jarvis.db.figshare.data("dft_3d")` (~76k entries)
2. Filters for structures containing the target ion element
3. For each candidate, iteratively removes ion atoms one at a time
4. Computes energy at each step using the ALIGNN ML force field (`alignnff_wt01` weights from figshare)
5. Calculates voltage: `V = (E_removed - E_full - E_unary_ion) / charge`
6. Exports results as JSON, converted to CSV via `json_to_csv.py`

**Output:** `Li.csv`, `Na.csv`, `K.csv`, `Ca.csv`, `Mg.csv`, `Al.csv`, `Zn.csv` -- one per ion type.

### Stage 2: Voltage Filtering (`average_voltage/`)

**Purpose:** Apply physics-based thresholds to narrow candidates from Stage 1.

**Filters (`Battery_Filtering.ipynb`):**
- **max_voltage <= 5.5 V** -- electrochemical stability window
- **ehull <= 0.05 eV** -- thermodynamic stability (energy above convex hull, from JARVIS-DFT metadata)
- **min_voltage > 1 V** -- practical operating range

**Output:** `Li_min.csv`, `Na_min.csv`, etc. -- filtered candidate lists used as input to Stage 3.

### Stage 3: Cathode Screening & DFT Preparation (`screening_cathode/`)

Two sub-stages:

#### Stage 3a: Cathode Ranking (`screen_cathode.py`)

**Purpose:** Rank Li-ion cathode candidates from `Li_min.csv` by a composite score.

**Pool:** only the `Li_` rows of `Li_min.csv` (7,193 JIDs); rows for other working ions on the same structure are dropped.

**Filters:**
- avg_voltage: 3.0--4.5 V
- q_grav > 100 mAh/g (theoretical gravimetric capacity per formula unit, `n_Li*F/(3.6*M)`; the original `max_grav_cap` column, which scales with cell size, is kept as `max_grav_cap_cell`)
- ehull <= 0.05 eV/atom
- max_voltage <= 5.5 V
- formula contains a redox-active transition metal (Ti V Cr Mn Fe Co Ni Cu Nb Mo Ru Rh W)

**Composite score:** `(1/3)*norm(avg_voltage) + (1/3)*norm(q_grav) - (1/3)*norm(ehull)`; 682 candidates survive. The superseded per-cell ranking (71 entries) is archived as `analysis/cathode_candidates_ranked_v1_cellnorm.csv`.

**Output:** `cathode_candidates_ranked.csv`

#### Stage 3b: Sequential Supercell Delithiation (`dft_prep.py`)

**Purpose:** Generate VASP DFT inputs for precise voltage curves via iterative Li removal from supercells.

**Workflow loop:**
```
init -> (run DFT on HPC -> record energy -> next step) x N -> voltage curve
```

1. **`init JVASP-XXXXX`** -- Looks up material in JARVIS DB, builds auto-sized supercell (all lattice vectors >= 7 A), auto-detects functional (PBE vs optB88-vdW based on spacegroup), writes VASP inputs for step_00 (fully lithiated)
2. **Run DFT externally** on HPC -- user adds POTCAR (species from POTCAR_spec) and submits. DFT outputs (CONTCAR, OUTCAR, etc.) should be placed in a `results/` subdirectory within each step (e.g. `step_00_Li16/results/OUTCAR`). These outputs are gitignored — the repo only carries inputs (POSCAR, INCAR, KPOINTS, POTCAR_spec).
3. **`record JVASP-XXXXX <step> <energy>`** -- Logs total energy from OUTCAR into `energies.json`
4. **`next JVASP-XXXXX`** -- Reads CONTCAR from previous step, uses ALIGNN to rank all symmetry-inequivalent Li vacancies by energy, removes the lowest-energy Li, writes VASP inputs for the next step
5. **`voltage JVASP-XXXXX`** -- Computes V(x) curve with both raw step voltages and convex hull equilibrium voltages, saves plot

**Optional:** `static JVASP-XXXXX <steps>` -- generates TB-mBJ static calculation inputs for band gap determination at specified delithiation steps.

```bash
python dft_prep.py init JVASP-XXXXX                         # Create supercell + step_00 (auto-detect functional)
python dft_prep.py init JVASP-XXXXX --functional pbe        # Force PBE
python dft_prep.py init JVASP-XXXXX --functional optb88vdw  # Force optB88-vdW
python dft_prep.py init JVASP-XXXXX --max-atoms 200         # Cap supercell at 200 atoms
python dft_prep.py next JVASP-XXXXX                         # Generate next step from CONTCAR
python dft_prep.py record JVASP-XXXXX 0 -453.0              # Log DFT energy for step 0
python dft_prep.py record JVASP-XXXXX 0                     # Auto-read energy from results/OUTCAR
python dft_prep.py static JVASP-XXXXX 0 1                   # Generate TB-mBJ inputs for steps 0 and 1
python dft_prep.py voltage JVASP-XXXXX                      # Compute + plot voltage curve
python dft_prep.py voltage JVASP-XXXXX --e-li-metal -0.95   # Override Li metal energy
```

**Voltage formula:**
```
V = (E_delithiated - E_lithiated + n_li * E_li_metal) / n_li
```
where `E_li_metal` is computed in-house with the same Li_sv PAW, ENCUT, and k-point density as the cathode runs (see `dft_inputs/JVASP-913-Li/Li_sv_PBE/`, `Li_sv_optB88vdW/` and `Li_sv_optPBEvdW/`):
- **PBE:** `-1.9031 eV/atom`
- **optB88-vdW** (`GGA = BO`, `PARAM1 = 0.1833333333`, `PARAM2 = 0.22`): `-0.9778 eV/atom`
- **optPBE-vdW** (`GGA = OR`): `-0.9646 eV/atom`. Before 2026-09-28 the generator wrote this tag set and labelled it optB88-vdW; the original LCO chain (`JVASP-2017-LCO/`) and four prospective runs used it and have been rerun with true optB88-vdW (`*-B88/` directories).

The tabulated JARVIS value (`-0.925 eV/atom`) is an optB88-vdW energy; combined with PBE+U cathode energies it shifts every PBE voltage by about +1 V, which is why the reference is recomputed per functional (it differs from the in-house optB88-vdW value by only 0.05 eV/atom).

**Key algorithms:**
- **Supercell auto-sizing:** Per-axis multiplier = `ceil(7 / lattice_vector_length)`, capped at 300 atoms. Uses diagonal-only supercell matrices (JARVIS `make_supercell` limitation).
- **Vacancy selection:** ALIGNN ranks every Li site individually (no symmetry deduplication, since symmetry is broken after the first vacancy). Removes the lowest-energy Li at each step.
- **Convex hull voltage:** Formation energies `dE(x) = E(x) - x*E(1) - (1-x)*E(0)`, lower convex hull identifies thermodynamically stable compositions, equilibrium voltages computed from hull vertex pairs. The raw step voltages (computed per individual Li removal) are noisy because many intermediate Li configurations are metastable — they sit above the convex hull and would phase-separate into a mixture of neighboring stable compositions at thermodynamic equilibrium. The hull identifies which compositions are truly stable; between two hull vertices the system exists as a two-phase mixture at constant voltage, producing the flat plateaus seen in experimental galvanostatic discharge curves. This is why the equilibrium hull voltage, not the raw step voltage, is the correct quantity to compare against experiment. For example, LMO produces two hull plateaus at 4.00 V and 4.17 V matching the known two-step discharge of spinel LiMn₂O₄, and LCO (optB88-vdW+U) produces plateaus at 4.18/4.23/4.62/4.78 V, reproducing the measured rise between its two two-phase regions (3.92 and 4.50 V) at a uniform offset of about +0.27 V.

**Functional auto-detection:**
- Layered cathodes (spacegroups 166/R-3m, 194/P63/mmc, 12/C2/m, 15/C2/c) -> **optB88-vdW** (PBE overestimates interlayer spacing in van der Waals bonded layers)
- All others (olivines, spinels) -> **PBE**
- Override with `--functional pbe|optb88vdw|optpbevdw`

**DFT settings:**
- **PBE relaxation:** ENCUT=520, EDIFF=1E-4 (EDIFFG defaulted to EDIFF*10), PREC=Accurate, ISMEAR=0 SIGMA=0.05, ISPIN=2 with material-specific AFM MAGMOM, ISYM=0, LREAL=Auto, NCORE=8, KPAR=2, AMIX=0.2 BMIX=0.0001 AMIX_MAG=0.4 (gentle mixing for DFT+U on magnetic TM oxides). ISIF=3 for step 0 (full cell+ion relaxation); ISIF=2 for steps 1+ on non-layered materials (ions only — avoids Pulay stress artifacts); ISIF=3 kept for layered materials (c-axis must relax on delithiation).
- **DFT+U (Dudarev, LDAUTYPE=2):** Applied to PBE relaxations. Element-specific U_eff values from [Materials Project](https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/gga+u-calculations/hubbard-u-values): Mn=3.9, Fe=5.3, Co=3.32, Ni=6.2, V=3.25, Cr=3.7, Mo=4.38, W=6.2 eV.
- **optB88-vdW tags:** GGA=BO, PARAM1=0.1833333333, PARAM2=0.22, LUSE_VDW=.TRUE., AGGAC=0.0 (appended to PBE base INCAR). GGA=OR with the same vdW tags is optPBE-vdW (`--functional optpbevdw`).
- **TB-mBJ static:** METAGGA=MBJ, ALGO=All, ICHARG=2, NELM=1000, no DFT+U.
- **KPOINTS:** Gamma-centered, scaled inversely with supercell: `max(1, round(3/n_i))` per axis.
- **POTCAR:** PAW labels from JARVIS `default_potcars.json` (Li_sv, Mn_pv, Fe_pv, Co, Ni_pv, etc.).

**Capacity formulas:**
- **Gravimetric:** `Q = n_Li * F / (3.6 * M)` (mAh/g), where F = 96485 C/mol and M = molecular weight of the host.
- **Volumetric:** `Q = n_Li * F / (3.6 * V * N_A)` (mAh/cm³), where V = cell volume (cm³).
- Gravimetric capacity is purely compositional (depends only on atomic masses and number of Li), so ALIGNN-FF and DFT values are identical. Volumetric capacity differs between methods because DFT relaxes the cell volume at each delithiation step while ALIGNN-FF uses the unrelaxed JARVIS structure throughout.

**Current materials under study:**

| JID | Formula | Abbreviation | Structure | Functional | DFT avg V | Experiment | Status |
|-----|---------|-------------|-----------|------------|-----------|------------|--------|
| JVASP-42723 | LiFePO4 | LFP | Olivine | PBE | 3.60 V | 3.43 V | **complete** |
| JVASP-116897 | LiMnPO4 | LMP | Olivine | PBE | 3.91 V | ~4.1 V | 16/17 (step_16, the fully delithiated MnPO₄ endpoint, did not converge) |
| JVASP-141792 | LiMn2O4 | LMO | Spinel | PBE | 4.08 V | 4.05 V | **complete** |
| JVASP-144791 | Li4Mn3Co2Ni3O16 | NMC variant (Li:TM = 1:2) | Layered (P1 in JARVIS) | PBE | 4.40 V (4.18 V with the cell relaxed at x = 0.5 and 0; optB88-vdW 4.76 / 4.45 V) | ~3.7 V (NMC-111) | **complete** |
| JVASP-2017 | LiCoO2 | LCO | Layered | optB88-vdW | 4.21 V over x = 1→0.5 (4.45 V over full delithiation; PBE+U 3.81 / 3.84 V) | ~4.05 V over x = 1→0.5 | **complete** (`JVASP-2017-LCO-B88/`) |

**Functional choice rationale:**
- **LFP, LMP (olivine), LMO (spinel):** PBE+U. These are 3D-bonded frameworks with no van der Waals gaps. Standard PBE with Hubbard U correction on the transition metal d-orbitals is sufficient and well-benchmarked for these structure types.
- **NMC variant (layered):** PBE+U. The JVASP-144791 entry is layered in structure but is stored in JARVIS-DFT with P1 symmetry, so the spacegroup rule routes it to PBE+U with fixed-cell steps after step 0. An optB88-vdW+U rerun (`JVASP-144791-NMC-vdW/`) and cell-relaxed endpoints for both functionals (`extra_2026-10-01/`) show that the vdW functional raises the voltage by 0.25–0.38 V, that the fixed cell accounts for about 0.2 V of the residual, and that the jump at x = 0.5 is present in every variant.
- **LCO (layered):** optB88-vdW+U. With PBE+U on the same supercell the voltage profile is flat (3.82 and 3.81 V over the two measured two-phase regions) and misses the 0.58 V rise that optB88-vdW+U reproduces, so the dispersion correction is needed for the Li-poor, CoO₂-like states.

### Periodic Trend Visualization (`periodic_trend/`)

**Purpose:** Interactive periodic table that color-codes elements by aggregated properties from cathode screening CSVs. Shows per-element statistics across the candidate pool -- for example, mean voltage of all materials containing a given element, or how many candidates contain each element.

**Usage (run from `periodic_trend/` directory):**
```bash
cd batterymat/periodic_trend/
python ptable.py ../screening_cathode/cathode_candidates_ranked.csv -p avg_voltage           # Mean voltage per element
python ptable.py ../screening_cathode/cathode_candidates_ranked.csv --agg count              # Element frequency (682 candidates, 40 non-Li elements)
python ptable.py ../screening_cathode/cathode_candidates_ranked.csv -p q_grav --agg max # Best capacity per element
python ptable.py ../screening_cathode/cathode_candidates_ranked.csv -p ehull --log -o ehull.html  # Log scale, custom output
```

**Arguments:**
- `csv_path` -- positional: path to CSV (e.g. `../screening_cathode/cathode_candidates_ranked.csv`, `../../average_voltage/Li_min.csv`)
- `-p` / `--property` -- column to aggregate: `avg_voltage`, `max_voltage`, `q_grav`, `max_grav_cap_cell`, `max_vol_cap`, `ehull`, `optb88vdw_bandgap`, `score` (required unless `--agg count`)
- `--agg` -- aggregation function: `mean` (default), `median`, `max`, `min`, `count`
- `--log` -- log color scale
- `-o` / `--output` -- output HTML file (default: `ptable.html`)
- `--include-li` -- include Li in the visualization (excluded by default since it appears in all 682 candidates and dominates the color scale)

**Output:** Interactive HTML periodic table (Bokeh + Plasma colormap). Elements with data are color-coded by the aggregated property value; elements without data are grayed out. Hover tooltips show element name and value.

**Dependencies:** `bokeh`, `bokeh_sampledata`, `matplotlib`.

### Stage 4: NEB Barriers (`neb_calc/`)

**Purpose:** Nudged Elastic Band (NEB) calculations for ion migration barrier energies -- determines how easily ions can move through the crystal lattice.

**Materials:** NiMnO and V2O5, with POSCAR image files for multiple ion types (Li, Na, K, Ca, Mg, Al, Zn). Each directory contains intermediate images along migration paths between equivalent sites. Uses Henkelman's NEB method scripts.

### Stage 5: Benchmarking (`benchmarks/`)

**Purpose:** Validate computational predictions against experimentally known cathode materials.

**Materials:** NMC, V2O5, TiSe2, MnO2, Chevrel phase (Mo6S8). DFT results stored as POSCARs under `PBE/` and `PBE-D3/` subdirectories. These provide reference voltage curves and structural parameters for comparison.

## Directory Structure

```
batterymat/
├── screening_anode/          # Stage 1: ML screening across 6 ion types
│   ├── screen.py             # Main screening script
│   ├── json_to_csv.py        # Convert JSON results to CSV
│   └── *.csv                 # Results per ion type
├── average_voltage/          # Stage 2: Threshold filtering
│   ├── Battery_Filtering.ipynb
│   └── *_min.csv             # Filtered candidates
├── screening_cathode/        # Stage 3: Cathode DFT pipeline
│   ├── screen_cathode.py     # 3a: Ranking by composite score
│   ├── dft_prep.py           # 3b: Sequential delithiation
│   ├── analysis/             # All plots and analysis
│   │   ├── voltage_curve_JVASP-XXXXX.png   # Per-material line plots
│   │   ├── discharge_curve_JVASP-XXXXX.png # Per-material staircase plots
│   │   ├── voltage_summary.png             # ALIGNN-FF vs DFT vs Experiment (avg/max V)
│   │   ├── capacity_summary.png            # ALIGNN-FF vs DFT vs Experiment (vol. capacity)
│   │   ├── make_summary_plots.py           # Regenerates the two summary plots
│   │   └── screening_cathode_analysis.md   # DFT hull voltage + vol. capacity vs experiment
│   └── dft_inputs/           # VASP input directories
│       └── JVASP-XXXXX-ABBREV/
│           └── supercell_NxNxN/
│               ├── energies.json
│               ├── step_00_LiXX/   (POSCAR, INCAR, KPOINTS, POTCAR_spec)
│               ├── step_01_LiXX/   # DFT outputs (CONTCAR, OUTCAR, etc.) live
│               ├── tmbj_step_00_LiXX/  #   locally in each step's results/ dir
│               └── ...             #   but are gitignored (not on GitHub)
├── periodic_trend/           # Interactive periodic table visualization
│   └── ptable.py             # CLI tool + plotting functions
├── neb_calc/                 # Stage 4: Ion migration barriers (NEB)
└── benchmarks/               # Stage 5: Validation vs known materials
tests/
└── test_dft_prep.py          # Tests for dft_prep.py
```

## Key Reference Values & Sources

| Parameter | Value | Source |
|-----------|-------|--------|
| E_li_metal (PBE) | -1.9031 eV/atom | In-house: `dft_inputs/JVASP-913-Li/Li_sv_PBE/`, Li_sv PAW, ENCUT=520, k=17×17×17 |
| E_li_metal (optB88-vdW) | -0.9778 eV/atom | In-house: `dft_inputs/JVASP-913-Li/Li_sv_optB88vdW/`, Li_sv PAW, GGA=BO + PARAM1/PARAM2, ENCUT=520, k=17×17×17 |
| E_li_metal (optPBE-vdW) | -0.9646 eV/atom | In-house: `dft_inputs/JVASP-913-Li/Li_sv_optPBEvdW/`, GGA=OR (original LCO chain only) |
| DFT+U values | Element-specific (see table above) | [Materials Project Hubbard U values](https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/gga+u-calculations/hubbard-u-values) |
| Crystal structures | JARVIS-DFT database (~76k materials) | [JARVIS-DFT](https://jarvis.nist.gov/), accessed via `jarvis.db.figshare.data("dft_3d")` |
| PAW potentials | JARVIS `default_potcars.json` | `jarvis.io.vasp.inputs` module |
| ALIGNN weights | `alignnff_wt01` | [ALIGNN (usnistgov)](https://github.com/usnistgov/alignn) -- weights auto-downloaded from figshare via `wt01_path()` |

## External Data

- Google Drive shared data: https://drive.google.com/drive/u/2/folders/1vf42QPRKc_kFXo5kmnWvDKzRFErUDdU3
- Google Colab demo: https://colab.research.google.com/gist/knc6/e7e3aecf3b748c0df1d58dc7911da950/
- Theoretical capacities reference: https://github.com/ndrewwang/BotB

## Literature References

1. Aydinol, M.K. et al. (1997). Ab initio calculation of the intercalation voltage of lithium-transition-metal oxide electrodes for rechargeable batteries. *J. Power Sources*, 68, 664. ([link](https://ceder.berkeley.edu/publications/jps-68-664-1997.pdf))
2. Choudhary, K. & DeCost, B. (2021). Atomistic Line Graph Neural Network for improved materials property predictions. *npj Computational Materials*. -- ML force field used in screening
3. Maxisch, T. et al. (2006). Modeling the Voltage Profile for LiFePO4. *J. Electrochem. Soc.* ([link](https://iopscience.iop.org/article/10.1149/1.2006607))
4. Wang, A. et al. (2021). First principles investigation into the interwoven nature of voltage and mechanical properties of the Li NMC-811 cathode. *J. Power Sources*. ([link](https://www.sciencedirect.com/science/article/pii/S0378775321011162))
5. Ab initio investigation of V2O5 for beyond lithium ion battery cathodes. *J. Power Sources* (2020). ([link](https://doi.org/10.1016/j.jpowsour.2020.228096))
6. Hybrid density functional theory modeling of Ca, Zn, and Al ion batteries using the Chevrel phase Mo6S8 cathode. *Phys. Chem. Chem. Phys.* (2017). ([link](https://doi.org/10.1039/c7cp03378h))
7. Density Functional Theory Modeling of MnO2 Polymorphs as Cathodes for Multivalent Ion Batteries. *J. Phys. Chem. C* (2018). ([link](https://doi.org/10.1021/acs.jpcc.8b00918))
8. Studies of Spinel-to-Layered Structural Transformations in LiMn2O4. *J. Phys. Chem. C* (2017). ([link](https://pubs.acs.org/doi/10.1021/acs.jpcc.7b00929))
9. TiSe2 cathode for beyond Li-ion batteries. *J. Power Sources* (2019). ([link](https://www.sciencedirect.com/science/article/pii/S0378775319307980))
10. Tailoring the Morphology of LiCoO2: A First Principles Study. *Chem. Mater.* (2009). ([link](https://pubs.acs.org/doi/epdf/10.1021/cm9008943))
11. An Overview and Future Perspectives of Aluminum Batteries. *Adv. Mater.* (2016). ([link](https://onlinelibrary.wiley.com/doi/full/10.1002/adma.201601357))

