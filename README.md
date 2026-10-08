# batch_detector

Chance-corrected cohort and response effects, to look at before you pool
cohorts.

When you combine public cohorts to get more samples, the cohorts usually
differ from each other (sequencer, site, protocol, patient mix) by more than
the biology you care about. `batch_detector.py` measures this. For the cohort
labels and for the response labels it runs PERMANOVA on Aitchison distances,
subtracts the amount of R² you would get from random labels, and reports the
two effects side by side with a GO / CAUTION / NO-GO pooling recommendation
and a power estimate. It is one Python file and works on any samples × features
table (microbiome counts, or other omics data that make sense under a CLR
transform or as plain Euclidean features).

Version 2.0.0. MIT license.

The `paper/` folder holds the scripts and result files behind the paper this
tool came from. You do not need it to use the tool.

## Install

```
pip install -r requirements.txt
```

Needs numpy, pandas, scipy and matplotlib. The versions pinned there are the
ones the paper's analysis ran with (Python 3.9). The tool and its tests also
pass on newer releases (checked with Python 3.13, numpy 2.1, pandas 2.2,
scipy 1.17, matplotlib 3.10), so the pins are a safe default, not a
requirement.

## Quick start

```
bash examples/run_example.sh
```

This runs the tool on a small synthetic dataset (120 samples, 3 cohorts, 50
features, a strong cohort effect and no response effect) and writes the
results to `outputs/example/`. It takes about 20 seconds. The only setting
changed from the defaults is `--n-grid 60,120,240`, which leaves out n = 480
in the power grid to save time. Everything else (999 permutations, 100 power
reps) is default.

The command it runs:

```
python batch_detector.py \
    --input  examples/example_features.tsv \
    --labels examples/example_metadata.tsv \
    --batch  cohort --clr --n-grid 60,120,240 \
    --output outputs/example
```

Expected result (the numbers are reproducible because the data and the seed
are fixed):

```
Response PERMANOVA : R2 = 0.0055   p = 0.660 (ns)
Batch PERMANOVA    : R2 = 0.4077   p = 0.001 (**)
Response delta_R2  : -0.00289  (E0=0.00840)
Batch delta_R2     : +0.39093  (E0=0.01681)
Batch:signal ratio : 73.99x
rho null percentile: 100.0
Recommendation     : NO-GO  (v2_chance_calibrated)
Min n (80% power)  : 240 same-protocol samples
```

Here `p = 0.001 (**)` is the smallest p-value 999 permutations can give, and
the label comes from the script's own cutoffs (it only prints `***` below
0.001). The synthetic data are made by `examples/make_example_data.py`.

## Input format

Both inputs are tab-separated, with a header row and sample IDs in the first
column.

**Feature table** (`--input`): one row per sample, one column per feature.

```
sample_id   taxon_01   taxon_02   ...
cohortA_01  1520       88         ...
cohortA_02  1311       143        ...
```

Counts or relative abundances, not negative. If you pass `--clr` the tool adds
a pseudocount and applies the centered log-ratio transform. Leave `--clr` off
if the table is already transformed (negative values are then allowed, but the
power analysis is skipped). Missing values are set to 0 and features with zero
variance are dropped.

**Labels table** (`--labels`): one row per sample.

```
sample_id   cohort    response
cohortA_01  cohortA   R
cohortA_02  cohortA   NR
```

`--response-col` names the outcome column (default `response`). The outcome
can use any two labels. More than two labels still gives the PERMANOVA numbers
but no power analysis.

**Cohorts** (`--batch`): either the name of a column in the labels table
(`--batch cohort`), or the path to a separate two-column TSV with the sample ID
first and the cohort second. Leave `--batch` out for a single-cohort run,
which reports the response effect and power only.

Samples are matched by ID across the files. Samples missing a label are
dropped and the tool logs how many.

## Output

These files go into the `--output` directory (the power file depends on the mode).

| file | contents |
|---|---|
| `batch_detector_report.json` | all statistics and the recommendation |
| `batch_detector_figure.pdf` | R² bars (left) and, with `--legacy-power`, the power curve (right). In the default mode the right panel only points to the power grid file. |
| `batch_detector_power_grid.tsv` | default power results, one row per (n, target ΔR²) |
| `batch_detector_power.tsv` | only with `--legacy-power`: power per n |

Fields of the JSON report:

| field | meaning |
|---|---|
| `tool`, `version`, `timestamp` | which version produced the report and when |
| `input_file`, `labels_file`, `response_col`, `batch_arg` | the arguments you passed |
| `n_samples`, `n_features`, `n_batch_groups` | size of the data after alignment and filtering |
| `clr_applied`, `pseudocount`, `zero_handling` | how the data were transformed |
| `permanova_n_perms` | number of permutations |
| `response_R2`, `response_F`, `permanova_p_response` | PERMANOVA of the response |
| `response_E0` | chance level of R², (g − 1) / (n − 1) for g groups |
| `response_delta_R2` | `response_R2 − response_E0`, the excess over chance |
| `batch_R2`, `batch_F`, `permanova_p_batch`, `batch_E0`, `batch_delta_R2` | the same for the cohort labels |
| `batch_signal_ratio` | `batch_R2 / response_R2` (the raw ratio, shown for reference) |
| `joint_margin_model` | cohort + response in one model: marginal R² and p for each term. Response is permuted within cohorts only. |
| `rho_permutation_null` | the observed ratio `rho`, its null when both label sets are shuffled independently, and `percentile_of_observed` in that null |
| `recommendation_version` | `v2_chance_calibrated` (default) or `legacy_v1` |
| `pooling_recommendation`, `recommendation_reason` | GO, CAUTION or NO-GO, and a sentence on why |
| `power_analysis_version` | `v2_hypothetical_effect_grid` (default) or `legacy_v1` |
| `min_n_for_80pct_power` | smallest n on the grid with at least 80% power (see below), or null |
| `n_for_80pct_power_by_target_delta_R2` | the same, for each target ΔR² (v2 only) |
| `dirichlet_bootstrap_n_reps`, `simulation_validation` | only filled with `--legacy-power` / `--simulate` |
| `permanova_response_detail`, `permanova_batch_detail` | sums of squares and permutation summaries |

**Why ΔR² and not R².** PERMANOVA R² is above zero even for random labels,
and by more when there are more groups or fewer samples. The excess over the
chance level is the number to compare. The raw ratio `batch_signal_ratio` also
grows with n when there is no response effect at all, which is why the tool
checks it against a shuffle null (`rho_permutation_null`).

**Reading the recommendation** (default rule):

| call | when | meaning |
|---|---|---|
| NO-GO | response ΔR² ≤ 0, or response p ≥ 0.05, or rho > 4, or rho is at or above the 95th percentile of its shuffle null | either there is no response effect to gain from pooling, or the cohort effect is much larger than the response effect |
| CAUTION | significant response effect and 0.25 < rho ≤ 4 | pooling may work, but check any batch correction yourself and run the tool again on the corrected data |
| GO | significant response effect and rho ≤ 0.25 (or a single cohort) | cohort effect is small compared with the response effect |

The cutoffs 0.25 and 4 are conservative rules of thumb. They were not
derived from the paper's simulation, which found no clean ratio-based
threshold: batch correction only helped at the highest simulated signal
level, whatever the ratio. Treat the call as a prompt to look at the numbers,
not as a verdict.

`min_n_for_80pct_power` comes from the power grid. The default (v2) grid asks
how many same-protocol samples would be needed to detect a response effect of
a given size (target ΔR²) with 80% power. It is the n found for the smallest
target ΔR² that reaches 80% on the grid. It stays null if no grid point
reaches 80%.

## Main options

```
python batch_detector.py --help
```

| option | what it does |
|---|---|
| `--clr` | apply the CLR transform before computing distances |
| `--zero-handling {pseudocount,multiplicative_replacement}` | how zeros are replaced before CLR. Default is a flat pseudocount; `--pseudocount` sets its size (or the delta for multiplicative replacement). Default 1e-6. |
| `--legacy-power` | use the v1 power analysis, which resamples the effect actually observed at several n. The default v2 analysis uses hypothetical effect sizes, which does not go circular when the observed effect is near zero. This is also the only mode that fills the power panel of the figure. |
| `--legacy-recommendation` | use the v1 rule: raw batch:signal ratio with cutoffs 3 and 8 and no chance correction |
| `--simulate` | with `--legacy-power`, repeat the power estimate at the current n as a check |
| `--n-perms N` | PERMANOVA permutations (default 999) |
| `--n-grid`, `--delta-grid` | n values and target ΔR² values for the v2 power grid |
| `--power-reps`, `--power-perms-inner`, `--calib-reps`, `--calib-iters` | size of the v2 power simulation |
| `--seed` | random seed (default 42) |

Known rough edges:

- The power analysis treats the label that sorts first as the positive class.
  It does not change the power estimate.
- The Dirichlet power simulation needs non-negative inputs.

## Runtime

Default settings, Apple M3 Pro, 18 GB RAM:

| data | samples | wall clock |
|---|---|---|
| 4 melanoma/GI microbiome cohorts (paper) | 283 | 3.6 min (217.2 s) |
| 4 colorectal cancer 16S cohorts (paper) | 570 | 13.0 min (781.1 s) |

The cost is dominated by the power grid, so for a first look you can lower
`--power-reps` and `--n-grid`.

## Tests

```
pip install -r requirements-dev.txt
pytest tests
```

The tests take a few seconds.

## Citation

If you use this tool, please cite the paper:

- Garg A. *Cohort Effects, Not Treatment Response, Dominate Public
  Immunotherapy Cohorts: Evidence from Microbiome and Tumor Genomic Data.*
  arXiv:XXXX.XXXXX.
- A shorter version is to appear in the IEEE BIBM 2026 Workshop Proceedings.

`CITATION.cff` has the software citation. A Zenodo DOI will be added after the first release.

## License

MIT, see `LICENSE`.
