# periodic_trend

Interactive periodic table that color-codes elements by aggregated properties from cathode screening CSVs. Forked from [arosen93/ptable_trends](https://github.com/arosen93/ptable_trends), updated for modern pandas/bokeh.

**Dependencies:** `bokeh`, `bokeh_sampledata` (periodic table element positions, split out from bokeh in newer versions), `matplotlib` (Plasma colormap + color normalization).

## CLI

Run from this directory:

```bash
cd batterymat_jae/periodic_trend/
python ptable.py ../screening_cathode/cathode_candidates_ranked.csv -p avg_voltage           # Mean voltage per element
python ptable.py ../screening_cathode/cathode_candidates_ranked.csv --agg count              # Element frequency across 682 candidates
python ptable.py ../screening_cathode/cathode_candidates_ranked.csv --agg count --include-li # Include Li in count
python ptable.py ../screening_cathode/cathode_candidates_ranked.csv -p q_grav --agg max # Best capacity per element
python ptable.py ../screening_cathode/cathode_candidates_ranked.csv -p ehull --log -o ehull.html # Log scale, custom output
```

## Arguments

- `csv_path` — positional: path to CSV file (e.g. `../screening_cathode/cathode_candidates_ranked.csv`, `../../average_voltage/Li_min.csv`)
- `-p` / `--property` — CSV column to aggregate (required unless `--agg count`). Available columns in `cathode_candidates_ranked.csv`: `avg_voltage`, `max_voltage`, `max_grav_cap`, `max_vol_cap`, `ehull`, `optb88vdw_bandgap`, `score`
- `--agg` — aggregation: `mean` (default), `median`, `max`, `min`, `count`
- `--log` — log color scale (fails if any value is negative)
- `-o` / `--output` — output HTML file (default: `ptable.html`)
- `--include-li` — include Li (excluded by default since it appears in all 682 candidates and dominates the visualization)

## Element extraction

`aggregate_by_element()` parses elements from the `atoms` column (serialized dict with `elements` key, parsed via `ast.literal_eval()`). Falls back to regex extraction from the `name` column (format: `Li_JVASP-XXXXX_LiMnPO4.json`). The `atoms` column is present in `cathode_candidates_ranked.csv`; `Li_min.csv` uses the `name` fallback.

## Functions

- `plot_ptable_trend(data_elements, data_list, ...)` — renders the Bokeh periodic table. Unchanged from upstream except pandas `.loc[]` fixes and `p.width` update for modern bokeh. Supports `log_scale`, custom `bokeh_palette`, `alpha`, `cbar_height`. `save_plot=True` saves HTML; `save_plot=False` opens in browser.
- `aggregate_by_element(csv_path, prop, agg, include_li)` — returns `(elements_list, values_list)` tuple ready for `plot_ptable_trend()`.

## Python API

```python
from batterymat_jae.periodic_trend.ptable import aggregate_by_element, plot_ptable_trend

elems, vals = aggregate_by_element("../screening_cathode/cathode_candidates_ranked.csv", prop="avg_voltage")
plot_ptable_trend(data_elements=elems, data_list=vals)
```

## Data context

`cathode_candidates_ranked.csv` contains 682 cathode candidates with 41 unique elements (40 excluding Li). The visualization shows per-element aggregated statistics, not individual material properties. For example, `-p avg_voltage` with `--agg mean` shows the mean voltage across all materials containing each element.
