# Reproducing the paper

Scripts, processed inputs and result files behind

> Garg A. *Cohort Effects, Not Treatment Response, Dominate Public
> Immunotherapy Cohorts: Evidence from Microbiome and Tumor Genomic Data.*
> arXiv:XXXX.XXXXX. (A shorter version is to appear in the IEEE BIBM 2026
> Workshop Proceedings.)

You do not need any of this to use `batch_detector.py`.

## Setup

```
pip install -r requirements.txt     # run inside paper/
```

The pins are the versions the analysis ran with (Python 3.9, numpy 2.0.2,
pandas 2.3.3, scipy 1.13.1, scikit-learn 1.6.1, xgboost 2.1.4).
scikit-bio is optional.

**Run every script from this folder (`paper/`).** The scripts use paths
relative to the project root (`results/...`, `metadata/...`, `scripts/...`),
and `paper/` keeps that layout. For example:

```
cd paper
python scripts/revision/task3_zero_handling_sensitivity.py
```

Scripts write their outputs next to the files already shipped, so a rerun
overwrites the shipped copy with a fresh one. Compare them with `git diff`.

## Layout

```
paper/
├── scripts/            analysis scripts (scripts/revision/ = the revision and review rounds)
│   └── revision/figures/   figure scripts (several versions, see the figure table)
├── metadata/           response labels, Lee 2022 site labels, one download manifest
├── data/               DATA_SOURCES.md (origin and terms of every input), get_duvallet_crc.sh
└── results/
    ├── ml/             processed abundance tables (inputs) and results of the earlier analysis phases
    └── revision/       results of the revision rounds, final figures, and the working notes:
        ├── SUMMARY.md                  every number added during the revision, with its source file
        ├── final/FINAL_CHECKS.md       last checks before the paper was locked (incl. why the rho cutoffs are heuristic)
        ├── figures/FIGURES.md          log of every figure, with sources and fixes
        ├── review/REVIEW_TASKS.md      first review round
        └── review2/REVIEW2_TASKS.md    second review round
```

`SUMMARY.md` and the `*_TASKS.md` / `*_CHECKS.md` files are working notes that
were written while the analysis was going on. They are kept because they
record where every number came from.

## Data: what is and is not in this folder

Included: the genus-level tables for the four microbiome cohorts (computed
by us from public reads), the response labels, the Lee 2022 site labels, and
all result files. See `data/DATA_SOURCES.md` for the source and terms of each.

**Not included**, and how to get them:

| missing input | needed by | how to get it |
|---|---|---|
| Raw FASTQ files (about 1 TB) | the Kraken2 step | downloaded from ENA by `scripts/process_batch.py` / `scripts/download_lee2022.sh`; accessions are in the table below |
| Kraken2 reports `results/kraken_reports/` | `build_matrix*.py`, `steps_1_4_c4_integration.py`, `phase3c_taxonomic_sensitivity.py`, `review_task4*`, `review2_task2`, `c4_sensitivity_*.py` | rerun the Kraken2 pipeline below |
| Tumor feature tables (Riaz 2017, Liu 2019, Hugo 2016) | `task1_*`, `task2_*`, `tumor_*.py`, `liu_single_cohort.py` | `python scripts/cbioportal_download.py`, then `python scripts/process_cbioportal.py` (needs network) |
| Colorectal cancer tables (Duvallet 2017) | `task1_*`, `task2_*` | `bash data/get_duvallet_crc.sh` |
| `metadata/response_sheet.xlsx` | `review2_task4_ena_lookup.py` | see `data/DATA_SOURCES.md` |

`task1_*` and `task2_*` read both the tumor and the colorectal tables, so
they only run after both have been rebuilt. The tumor and colorectal
tables were left out because of their licenses, not size
(`data/DATA_SOURCES.md`).

Accessions of the four microbiome cohorts:

| cohort | study | BioProject | platform |
|---|---|---|---|
| C1 | Frankel 2017 | PRJNA397906 | HiSeq 2000 |
| C2 | Peng 2020 (GI cancer) | PRJNA615114 | NovaSeq 6000 |
| C3 | Matson 2018 | PRJNA399742 | NextSeq 500 |
| C4 | Lee 2022 | PRJEB43119 | five sites |

## Taxonomic profiling (Kraken2)

All 283 samples went through the same two steps. Paired-end reads are
trimmed with fastp (default options), then classified with Kraken2 against the
`k2_standard_08gb` index.

- Kraken2 v2.17.1.
- Database: `k2_standard_08gb`, the March 2023 build of the pre-built 8 GB
  capped standard index (<https://benlangmead.github.io/aws-indexes/k2>).
  `hash.k2d` has MD5 `9d57fbd9f8bcee3d646991df60ed1f66`; check yours against
  it.
- fastp: default settings. The version used for the earlier runs was not
  recorded; the one installed at the end of the project was 1.1.0.

For a sample `SRRxxxxxxx` (this is what `scripts/process_batch.py` does, with
the database files in the current directory):

```
fastp -i SRRxxxxxxx_1.fastq.gz -I SRRxxxxxxx_2.fastq.gz \
      -o SRRxxxxxxx_1.trimmed.fastq.gz -O SRRxxxxxxx_2.trimmed.fastq.gz \
      -h results/fastp/SRRxxxxxxx_fastp.html -j results/fastp/SRRxxxxxxx_fastp.json \
      --thread 4

kraken2 --db . --threads 4 --paired \
        --report results/kraken_reports/SRRxxxxxxx_report.txt \
        --output results/kraken_outputs/SRRxxxxxxx_output.txt \
        SRRxxxxxxx_1.trimmed.fastq.gz SRRxxxxxxx_2.trimmed.fastq.gz
```

The scripts delete the FASTQ files after each sample to save disk space. The
reports are then turned into genus tables with `scripts/build_matrix.py`
(C1), `scripts/build_matrix_3cohort.py` (C1 to C3) and
`scripts/steps_1_4_c4_integration.py` (adds C4 and writes the n = 283 tables).
For C4 use `scripts/download_lee2022.sh`; for C2 build the manifest with
`scripts/build_manifest.py --prefix SRR11413`.

## Tables

| Table | Content | Script | Result file(s) |
|---|---|---|---|
| I | cohorts | cohort assembly (`build_manifest.py`, `process_batch.py`, `steps_1_4_c4_integration.py`) | `results/ml/n283_4cohort/response_labels_n283.tsv` |
| II | chance-corrected PERMANOVA, n = 39 to 283 | `scripts/revision/task1_chance_calibrated_permanova.py` | `results/revision/task1_chance_calibrated_permanova.tsv` (rho null: `task1_rho_permutation_null.tsv`) |
| III | joint model, response permuted within cohort | `scripts/revision/task2_joint_model_restricted_perm.py` | `results/revision/task2_joint_model_margin_permanova.tsv` |
| IV | repeated stratified 5-fold CV | `scripts/revision/task5_repeated_cv.py` | `results/revision/task5_{c1,n118,n283,c4}_permutation_summary.tsv`, `..._fold_aucs.tsv` |
| V | response PERMANOVA within each cohort | `scripts/revision/review_task2_percohort_permanova.py` | `results/revision/review/percohort_permanova_C1_C4.tsv` |
| VI | zero-handling methods | `scripts/revision/task3_zero_handling_sensitivity.py` | `results/revision/task3_summary_range.tsv` (all rows: `task3_zero_handling_sensitivity.tsv`) |
| VII | taxonomic level | `scripts/revision/task1_chance_calibrated_permanova.py` (phylum and species tables come from `scripts/phase3c_taxonomic_sensitivity.py`) | `results/revision/task1_chance_calibrated_permanova.tsv` (phylum, genus, species rows at n = 118) |
| VIII | nested LOOCV | `scripts/nested_cv.py` (n = 118), `scripts/nested_cv_n283.py` (n = 283) | `results/ml/nested_cv/nested_cv_summary.tsv` and `permutation_test_nested.tsv`; `results/ml/n283_4cohort/nested_cv_n283_summary.tsv` and `nested_cv_n283_perm_summary.tsv`. The p-values 0.248 and 0.109 for n = 118 are `(n_perm_gte_obs + 1) / (n_perms + 1)` from `permutation_test_nested.tsv` (see `results/revision/review/REVIEW_TASKS.md`, item 5) |
| IX | six batch-correction methods | `scripts/phase1_correction_shootout.py` | `results/ml/batch_correction_shootout/method_comparison.tsv` |
| X | cross-cohort transfer | `scripts/phase05_cross_cohort_holdout.py` | `results/ml/cross_cohort/holdout_summary.tsv` (`results/revision/task7_cross_cohort_full_table.tsv` is a reformatted copy made by hand, no script) |
| XI | C4 sensitivity analysis | `scripts/c4_sensitivity_analysis.py`, `scripts/c4_sensitivity_v4v5.py` | `results/ml/lee2022/sensitivity/sensitivity_summary.tsv` |
| XII | simulated change in AUC | `scripts/phase2_simulation.py` | `results/ml/simulation/grid_results.tsv` (3 cohorts, 40 per cohort, cohort R² = 0.08) |
| XIII | simulated power in C4 | `scripts/revision/task8b_c4_power_v2.py` | `results/revision/task8_c4_power_v2.tsv`; upper bound in `task8_c4_upper_bound.tsv`; 30-genus check from `scripts/revision/task8c_c4_power_top30.py` in `results/revision/final/task8c_c4_power_top30.tsv` |

Other numbers quoted in the text:

| quantity | script | result file |
|---|---|---|
| depth and classification rate per cohort | `scripts/revision/review_task4_depth_classification.py` | `results/revision/review/depth_classification_summary.tsv` |
| rarefaction to 2,098,257 reads | `scripts/revision/review_task4b_rarefaction.py` | `results/revision/review/rarefied_permanova.tsv` |
| rarefaction to 5,000,000 reads | `scripts/revision/review2_task2_rarefaction_5M.py` | `results/revision/review2/rarefied_5M_permanova.tsv`, `dropped_below_5M.tsv` |
| PERMDISP, recruiting site in C4 | `scripts/revision/review2_task1_permdisp_c4_site.py` | `results/revision/review2/permdisp_c4_site.tsv` |
| n = 281 after dropping two samples | `scripts/revision/review2_task3_drop_two_samples.py` (written afterwards; it reproduces the shipped file exactly) | `results/revision/review2/permanova_n281_dropped2.tsv` |
| C4 alone, recruiting site | `scripts/revision/task4_c4_alone.py` | `results/revision/task4_c4_alone_permanova.tsv` |
| genus meta-analysis (Hedges' g, DerSimonian-Laird) | `scripts/revision/task6_meta_analysis_hedges_g.py` | `results/revision/task6_meta_analysis_hedges_g.tsv` |
| sign flips of elastic-net coefficients (64.6% vs 58.8%) | `scripts/phase3b_analysis.py` | `results/ml/phase3b/inversion_permutation_test.tsv` |
| Dirichlet bootstrap CIs, PERMDISP for batch | `scripts/phase0_dirichlet_bootstrap.py`, `scripts/statistical_foundations.py`, `scripts/steps_1_4_c4_integration.py` | `results/ml/batch_diagnostics/` |
| tumor analyses | `scripts/tumor_ml.py`, `scripts/tumor_batch.py`, `scripts/liu_single_cohort.py` | `results/ml/tumor/` |
| `batch_detector` on the n = 283 data and on the CRC data | `batch_detector.py` (top level) | `results/revision/task9_n283/`, `results/revision/task9_duvallet/` |

## Figures

The paper numbers its figures in the order they appear. The file names
come from the order in which the figures were first made, so the two do not
line up one to one.

| Paper figure | File | Script | Reads |
|---|---|---|---|
| Fig. 1 dose-response | `results/revision/figures/fig1_dose_response.{pdf,png}` | `scripts/revision/figures/fig1_dose_response.py` | `results/revision/task1_chance_calibrated_permanova.tsv` |
| Fig. 2 ordination | `results/revision/figures/fig2_ordination_v2.{pdf,png}` | `scripts/revision/figures/fig2_ordination_v2.py` | the CLR tables for n = 283 and C4, their labels and sites, `task1_...tsv`, `task4_c4_alone_permanova.tsv` |
| Fig. 3 CV against permutation null | `results/revision/figures/fig5_cv_vs_null_v2.{pdf,png}` | `scripts/revision/figures/fig5_cv_vs_null_v2.py` | `results/revision/task5_*_{fold_aucs,permutation_aucs,permutation_summary}.tsv` |
| Fig. 4 meta-analysis volcano | `results/revision/final/fig6_volcano_v2.{pdf,png}` | `scripts/revision/figures/fig6_volcano_v2.py` | `results/revision/task6_meta_analysis_hedges_g.tsv` |
| Fig. 5 excess variance, all datasets | `results/revision/figures/fig3_chance_calibrated_summary_v2.{pdf,png}` | `scripts/revision/figures/fig3_chance_calibrated_summary_v2.py` | `task1_...tsv`, `task4_c4_alone_permanova.tsv` |
| Fig. 6 simulation grid | `results/revision/figures/fig7_simulation_grid_v3.{pdf,png}` | `scripts/revision/figures/fig7_simulation_grid_v3.py` | `results/ml/simulation/grid_results.tsv`, `calibration_deltas.tsv` |
| Fig. 7 power in C4 | `results/revision/final/fig4_power_curves_v5.{pdf,png}` | `scripts/revision/figures/fig4_power_curves_v5.py` | `results/revision/task8_c4_power_v2.tsv` |

The other versions of the figure scripts (`fig*_v1`, `_v2`, ...) are earlier
drafts and are kept for the record; `results/revision/figures/FIGURES.md` says
what changed between versions. The figure scripts only read result files, so
they run on the shipped data.

## What runs on the shipped data

Each script was run once from a fresh copy of this folder in a clean Python
3.9 environment with the pinned versions (90 second limit per script). Result:

| outcome | scripts |
|---|---|
| **Finished** (seconds to a minute) | `phase05_cross_cohort_holdout`, `phase3b_feature_stability`, `statistical_foundations`, `revision/task3`, `task4`, `task6`, `task8c`, `review_task2`, `review2_task1`, `review2_task3`, `revision/task10`, and all `revision/figures/fig*.py` |
| **Started fine, long jobs** (did not finish in 90 s) | `nested_cv`, `nested_cv_c4`, `nested_cv_n283`, `phase0_dirichlet_bootstrap`, `phase1_correction_shootout`, `phase3b_analysis`, `c4_sensitivity_analysis`, `revision/task8`, `revision/task8b` |
| **Need arguments** | `revision/task5_repeated_cv.py --dataset {c1,n118,n283,c4}`; `phase2_simulation.py --mode {timing,calibrate,full}` (the full grid takes about 15 minutes) |
| **Stop at once because an input is not shipped** | needs the Kraken2 reports: `build_matrix`, `build_matrix_3cohort`, `steps_1_4_c4_integration`, `phase3c_taxonomic_sensitivity`, `c4_sensitivity_v4v5`, `revision/review_task4_depth_classification`, `review_task4b_rarefaction`, `review2_task2_rarefaction_5M`; needs the tumor and CRC tables: `revision/task1`, `task2`, `tumor_ml`, `tumor_batch`, `liu_single_cohort`; needs the raw cBioPortal JSON: `process_cbioportal` |

**Warning: do not run the scripts that need the Kraken2 reports without the
reports.** `steps_1_4_c4_integration.py`, `review_task4*.py` and
`review2_task2_rarefaction_5M.py` do not check that their inputs exist. They
write empty tables over `results/ml/n283_4cohort/*`, `results/ml/lee2022/*` and
the `results/revision/review*/` tables. If that happens, restore the files
with `git checkout -- results metadata`.

## Reproducing the shipped numbers: what to expect

The files regenerated by those runs were compared with the shipped ones.

- Identical to the shipped files: 101 TSV/JSON result files, including the
  tables written by `task3`, `task4`, `task6`, `task8c`, `review_task2` and
  `review2_task1`, the random-forest results of
  `phase05_cross_cohort_holdout.py`, and `phase3b_feature_stability.py`.
  `review2_task3` gives the same numbers as `permanova_n281_dropped2.tsv`
  (it writes to a separate `_rerun` file).
  `batch_detector.py` run on the n = 283 data reproduces
  `results/revision/task9_n283/batch_detector_report.json` in every field.
- **Elastic-net results differ in the third decimal.** For
  `phase05_cross_cohort_holdout.py` (Table X) the elastic-net AUCs moved by up
  to 0.006 and permutation p-values by up to 0.016 (for example C2 to C1: AUC
  0.550 shipped, 0.544 on rerun, p 0.338 and 0.354). Random-forest rows were
  identical. The same happens in `statistical_foundations.py` (the n = 39
  elastic-net AUC is 0.447 shipped and 0.450 on rerun). The cause was not
  tracked down (library build or solver convergence); no conclusion changes.
- **Run order matters** for a few shared tables:
  1. `statistical_foundations.py` writes Phase 0 PERMDISP rows for n = 39, 79
     and 118 and the effect-size columns of `permanova_comparison.tsv`. It also
     writes `bootstrap_cis.tsv` with the old **naive** bootstrap, which is
     biased (see the notes in `SUMMARY.md`).
  2. `phase0_dirichlet_bootstrap.py` then overwrites `bootstrap_cis.tsv` with
     the corrected Dirichlet intervals used in the paper.
  3. `steps_1_4_c4_integration.py` appends the n = 283 rows to
     `permdisp_results.tsv` and `permanova_comparison.tsv`.
  Running only step 1 therefore leaves the shipped n = 283 rows and the
  corrected intervals missing. Running steps 1 and 2 in a fresh copy gave
  corrected intervals within 0.003 of the shipped ones; the small difference
  comes from the n = 39 elastic-net AUC being recomputed in step 1 (see the
  previous point).

## Changes made to the scripts when they were copied here

- Scripts are unchanged except as listed here. No script contained an absolute path.
- `scripts/prepare_duvallet_input.py` is `results/batch_detector_external_validation/prepare_input.py` from the
  original project. It had `/tmp/microbiomeHD_crc` hard coded; it now reads
  `microbiomeHD_crc/` and writes `results/batch_detector_external_validation/input_feature_matrix.tsv` and
  `input_labels.tsv`, the file names the analysis scripts expect.
- `scripts/revision/review2_task3_drop_two_samples.py` is new (see above).
- The HUMAnN3 scripts, the older LOOCV scripts and the first drafts of the
  figure and meta-analysis code are not included; they did not produce numbers
  in the final paper.
