# RGrowthCurves

Semi-automatic quantification of yeast growth-curve parameters — maximum specific
growth rate (Mu), doubling time, lag phase and area under the curve (AUC) — from
multi-well plate-reader tables.

**Base R only, no packages required.**

## Repository contents

| File | Description |
|---|---|
| `RGrowthCurves.R` | The script (v0.6.0) |
| `inn.csv` | Sample input file |
| `inn_report.csv` | Sample output report produced from `inn.csv` |

## Usage

From the command line:

```bash
Rscript RGrowthCurves.R input.csv [interval_min]
```

- `input.csv` — growth table (see *Input format* below).
- `interval_min` — optional; minutes between reads. Only used when the file has
  no time column.

Alternatively, open the script in RStudio, edit the `SETTINGS` block at the top
and `source()` the file.

## Input format

A "vertical" table: **one column per well**, one row per read. The header row
contains well/condition names.

- An optional first column with time is recognised when its header starts with
  `Time`, `t`, or `Время` (Cyrillic). Time values may be given in minutes,
  hours, or as `h:mm[:ss]`; the unit is guessed from the header
  (e.g. `Time, h`, `Время (час)`).
- If there is no time column, reads are assumed to be spaced `INTERVAL_MIN`
  minutes apart.
- Separator (`;`, `,` or TAB) and decimal mark (point or comma) are detected
  automatically.
- Encodings: UTF-8 (with or without BOM) or Windows-1251 (typical for Excel on
  Russian Windows).
- The strings `NA`, `OVRFLW`, `OVER`, `-`, `*` and empty cells are treated as
  missing values and skipped during fitting.

## Output

All output is written to the `RGrowthCurves_output/` folder:

| Output | Description |
|---|---|
| `<input>_report.csv` | One row per well with all parameters (UTF-8 with BOM, opens correctly in Excel) |
| `plots/NNN_<well>.png` | One diagnostic plot per well, for manually checking the fit |
| `<input>_overview_NN.png` | Overview pages, 4×4 wells per page (grid is configurable), all panels on the same y-scale |

### Report columns

| Column | Meaning |
|---|---|
| `conditions` | Well/condition name (column header) |
| `Mu` | Maximum specific growth rate, 1/h — steepest slope of log2(OD) among all "usable" sliding windows |
| `tau12, h` | Doubling time, h |
| `starting OD` | Mean OD in the `START_OD_MIN` range |
| `OD10h_minusStartingOD`, `OD10h_abs` | OD at 10 h (±30 min average), relative to the starting OD and absolute; the time point is set by `OD_AT_H` |
| `lag, min` | Start of the maximum-slope window |
| `lag_corrected, min` | Tangent lag: intersection of the exponential-phase fit line with the starting-OD level |
| `r2` | R² of the best window fit |
| `RGScript_version` | Script version used |
| `Mu_10`, `Mu_12`, `Mu_14` | Maximum Mu among windows **ending** before 10/12/14 h — i.e. Mu as if the plate had been stopped at that time (comparable between runs of different length). Empty if the run is shorter than the cutoff |
| `Mu_w250` | Mu computed with an extra 250-min window (the v0.4 default), for comparison with older runs |
| `AUC, OD*h` | Area under the OD(t) curve for the whole run (trapezoid rule) |
| `AUC_0-10h, OD*h` | AUC from 0 to `AUC_CUTOFF_H` (10 h) — comparable between runs |
| `AUC_net, OD*h`, `AUC_0-10h_net, OD*h` | AUC after subtracting the baseline (plate + medium + inoculum) |
| `baseline OD` | Lowest smoothed OD (running median of 5) from time 0 to the start of the exponential-phase fit (at least `BASELINE_MIN`) |
| `fit window, min` | Window width used for the main Mu |
| `Mu_cutoff_rule` | Rule used for the `Mu_10/12/14` columns (`end` or `start`) |
| `QC` | Quality-control flag (see below) |

### QC flags

- `ok` — no issues detected.
- `no growth (dOD<0.05)` — the well grew less than `MIN_DELTA_OD`; growth-rate
  columns are left empty.
- `no window with R2>=0.95` — growth was detected, but no window passed the
  R² threshold.
- `fit window reaches OD>1, check reader linearity` — the best window includes
  readings above OD 1, where many readers leave the linear range.
- `run shorter than X h` — the run ended before the largest `MU_CUTOFFS_H`
  value, so some `Mu_*` columns are empty.

### Diagnostic plots

Each `plots/*.png` shows, for one well:

- **black circles** — original OD values;
- **blue circles** — log2(OD);
- **green dots** — Mu ×10 in each sliding window (offset to the bottom of the
  plot; grey = windows rejected by the R² threshold);
- **orange dots** — R² of each window fit;
- **red line** — the best (maximum-slope) exponential fit, with the fit window
  shaded;
- **vertical blue line** — start of the maximum-slope window (`lag, min`);
- **red dot** — corrected lag (`lag_corrected, min`);
- **dashed green line** — log2 of the starting OD;
- the legend gives Mu (1/h), doubling time (h) and corrected lag (min).

The overview pages show raw OD curves with the fitted exponential segment in
red; wells with a QC flag are drawn in grey.



## Settings

All settings are in the `SETTINGS` block at the top of the script:

| Setting | Default | Description |
|---|---|---|
| `INPUT_FILE` | `inn.csv` | Input file (overridden by command-line argument 1) |
| `INTERVAL_MIN` | `5` | Minutes between reads (ignored if a time column is present; overridden by argument 2) |
| `WINDOW_MIN` | `150` | Sliding-window width for the main Mu, min |
| `EXTRA_WINDOWS_MIN` | `250` | Additional window widths, reported as `Mu_w250` etc. for comparison with older runs |
| `FIRST_FIT_MIN` | `0` | Windows may not start before this time, min |
| `MIN_R2` | `0.95` | Windows with lower R² are ignored when searching for the maximum slope |
| `BLANK_OD` | `0` | Known background OD subtracted before log2 and AUC (0 = off). Leave 0 if the plate has no empty-medium wells |
| `BASELINE_MIN` | `120` | Net-AUC baseline: lowest smoothed OD from 0 to the start of the exponential-phase fit (at least this many min); no blank wells needed |
| `START_OD_MIN` | `c(5, 35)` | Time range (min) averaged as the starting OD (v0.4 used rows 2–8) |
| `OD_AT_H` / `OD_AT_HALF_MIN` | `10` / `30` | Time point (h) for the `OD10h` columns and the ± minutes averaged around it |
| `MU_CUTOFFS_H` | `c(10, 12, 14)` | Time cutoffs for the `Mu_10` etc. columns |
| `MU_CUTOFF_RULE` | `"end"` | `"end"`: only windows *ending* before the cutoff count (Mu as if the plate had been stopped at that time — comparable between runs). `"start"`: windows *starting* before the cutoff (v0.4 behaviour; these windows use data up to cutoff + window width) |
| `AUC_CUTOFF_H` | `10` | AUC is also reported for 0…this time |
| `MIN_DELTA_OD` | `0.05` | Wells growing less than this are flagged `no growth` |
| `OUT_DIR` | `RGrowthCurves_output` | Output folder |
| `OUT_SEP` / `OUT_DEC` | `";"` / `"."` | Report separator and decimal mark (use `","` for Russian-locale Excel) |
| `MAKE_PLOTS` | `TRUE` | Set to `FALSE` to skip per-well and overview plots |
| `OVERVIEW_GRID` | `c(4, 4)` | Rows × columns of wells per overview page |

## Method notes

- **Mu** is the maximum slope of log2(OD) over all sliding windows of
  `WINDOW_MIN` minutes, converted to 1/h; only windows with positive slope and
  R² ≥ `MIN_R2` are considered, and wells flagged as "no growth" are excluded.
- **Lag** is reported twice: as the start of the maximum-slope window, and as
  the tangent lag — the intersection of the fitted line with the starting-OD
  level (`lag_corrected, min`).
- **AUC** uses the trapezoid rule, in OD·h. Net AUC subtracts a per-well
  baseline estimated as the lowest smoothed OD early in the run (running median
  of 5 reads removes single-read spikes from bubbles or condensation), so no
  blank wells are required. A known fixed background can instead be subtracted
  via `BLANK_OD`.
- Report column names are kept compatible with v0.4; new columns are appended
  at the end.

## Changes from v0.4

- New output layout: everything goes to `RGrowthCurves_output/`; the report is
  a CSV (`<input>_report.csv`) instead of `report.txt`.
- New metrics: AUC (whole run, fixed-cutoff, and baseline-corrected versions),
  `Mu_10/12/14` cutoff-restricted growth rates, OD at 10 h, per-well baseline.
- Per-well diagnostic PNGs plus 4×4 overview pages.
- Automatic detection of separator, decimal mark, encoding (UTF-8/Windows-1251)
  and time column (minutes, hours, `h:mm:ss`, Cyrillic headers).
- Missing-value tokens (`OVRFLW`, `OVER`, `-`, `*`) handled.
- QC column with explicit flags instead of silently empty results.
- Doubling time, corrected lag and script version are now reported.
- Base R only — no packages required.

## License

This script is released under the [MIT License](LICENSE): it is free to use,
copy, modify and distribute, also in commercial and closed-source projects, as
long as the copyright notice is retained. The software is provided "as is",
without warranty of any kind.

## Citation

The license places no restrictions on use, but if you use this script (in any
version) in work leading to a publication, **please cite the paper where an
earlier version of the script was used**:

> Kashko ND, Muravyov G, Karavaeva I, Glagoleva ES, Logacheva MD,
> Garushyants SK, Knorre DA (2025) Inheritance bias of deletion-harbouring
> mtDNA in yeast: The role of copy number and intracellular selection.
> *PLoS Genetics* 21(6): e1011737.
> https://doi.org/10.1371/journal.pgen.1011737
