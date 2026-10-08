# Review-response tasks

Rules followed: no existing result overwritten; new outputs in
`results/revision/review/`; every new number appended to `SUMMARY.md`
Section 14 with its source.

---

## 1. Cohort inclusion accounting

Full table: `results/revision/review/cohort_screening.tsv`. Source for the
11-study breakdown: `metadata/response_sheet.xlsx` (624 rows, matches the
Orletskaia & Olekhnovich compilation exactly: sum of the 11 per-study counts
= 624). Cancer type/sequencing-type/accession facts for the 7 NOT-included
studies were checked against public sources (WebSearch) where possible;
anything not confirmed with confidence is marked "Not confirmed via
available sources" rather than guessed. **The actual reason each of the 7
was excluded from this project is "not documented" for all seven** — a
project-wide grep for each study name across every `.py`/`.md`/`.tsv`/`.log`
file found zero references. This project's own history only documents an
*exhaustive search for non-melanoma cohorts* (which surfaces Gunjur 2024 as
a systematic-review source, not as a cohort it tried to download) — it does
not document evaluating or rejecting Gopalakrishnan 2018, Heshiki 2020,
Liu 2022, McCulloch 2022, Spencer 2021, or Tsakmaklis 2023 at all.

**C1-C4 accessions actually used** (verified against project metadata files
during this review, not just taken from earlier notes):

| cohort | study | accession | verified from |
|---|---|---|---|
| C1 | Frankel 2017 | PRJNA397906 | `metadata/sample_metadata.tsv`, `metadata/run_to_sample.tsv` |
| C2 | Peng 2020 | PRJNA615114 | **Not found in any project file** — see note below |
| C3 | Matson 2018 | PRJNA399742 | `metadata/PRJNA399742_wgs_runs.tsv` |
| C4 | Lee 2022 | PRJEB43119 | `scripts/download_lee2022.sh`, `scripts/steps_1_4_c4_integration.py` |

**Important finding on C2's accession**: a project-wide search found
`PRJNA615114` in **zero** files this project generated or downloaded — it
was never recorded anywhere in `metadata/` or `scripts/`. To close this gap,
The ENA API was queried directly during this review (`api.crossref.org` is for
papers; used `ebi.ac.uk/ena/portal/api/filereport` instead) for one of our
own C2 run accessions, `SRR11413606`: its `study_accession` is confirmed
**PRJNA615114**. So the number is correct, but it was not independently
verifiable from this project's own records before this check — worth adding
a one-line provenance note to the paper's methods.

Two small sample-count discrepancies surfaced and are flagged, not resolved
(would need more digging to reconcile): the Orletskaia compilation lists
Matson_2019 at 38 samples (this project's C3 has 39) and Lee_2022 at 164
(this project's C4 has 165). Off-by-one in either direction; not chased down
further here.

**Paper-ready exclusion summary (3-4 sentences):** *"Of the 11 studies in the
Orletskaia & Olekhnovich (2026) compilation (624 samples total), four had raw
sequencing reads that this project located, downloaded, and processed as
independent cohorts (Frankel 2017, Peng 2020, Matson 2018, Lee 2022; C1-C4,
n=283 combined). The remaining seven studies (Gopalakrishnan 2018, Gunjur
2024, Heshiki 2020, Liu 2022, McCulloch 2022, Spencer 2021, Tsakmaklis 2023;
341 samples) were not evaluated for inclusion in this project's own
processing pipeline, and no record of that evaluation exists in project
history. At least one (Gopalakrishnan 2018) used a substantially 16S-based
protocol incompatible with this project's shotgun-WGS Kraken2 pipeline;
sequencing type, public availability, and FMT status for the others were not
confirmed with confidence from public sources in this review and would need
independent verification before a claim of completeness."*

---

## 2. Per-cohort PERMANOVA

Script: `scripts/revision/review_task2_percohort_permanova.py`. Output:
`results/revision/review/percohort_permanova_C1_C4.tsv`. Same settings as
Task 4 (genus level, Aitchison distance, CLR pseudocount=1e-6, 999 perms).
C1 and C4 pulled from already-computed `task1_chance_calibrated_permanova.tsv`
/ `task4_c4_alone_permanova.tsv` (not recomputed, per the "don't overwrite /
don't needlessly recompute" spirit of the rules).

| cohort | n | R2 | E0 | delta_R2 | p |
|---|---|---|---|---|---|
| C1 | 39 | 0.026745 | 0.026316 | +0.000430 | 0.419 |
| C2 | 40 | 0.027937 | 0.025641 | +0.002296 | 0.275 |
| C3 | 39 | 0.023814 | 0.026316 | -0.002502 | 0.635 |
| C4 | 165 | 0.006631 | 0.006098 | +0.000533 | 0.277 |

All four cohorts, individually, show response ΔR² within ±0.0025 of zero and
p > 0.27 — the null holds within every single cohort on its own, not just
after pooling.

---

## 3. Strict response definition (CR/PR vs PD)

**Outcome: labels could not be obtained for any of the 4 cohorts in this
pass — reported honestly rather than approximated.** Details per cohort:

- **C1 (Frankel 2017)**: confirmed via PMC (PMC5602478) that the paper's own
  **Table 2** reports per-patient RECIST category (Response
  [CR/PR combined] / Stable / Progression) for all 39 patients — usable in
  principle. **Blocker**: this project has no patient-ID-to-run-accession
  crosswalk (`metadata/run_to_sample.tsv` only has SRA sample/experiment
  accessions, not the paper's own patient IDs like "P7"). Extracting Table 2
  via automated tools risks transcription errors on a 39-row clinical table,
  and there is no way to verify the extraction against a local crosswalk.
  Not completed this pass; would need a manual pull of Table 2 plus an NCBI
  BioSample-attribute lookup to link patient IDs to SRR accessions.
- **C2 (Peng 2020)**: per-patient Supplementary Table S1 is behind the AACR
  paywall (already documented in this project's prior session notes as
  inaccessible). Not obtained.
- **C3 (Matson 2018)**: the two local files that should hold this
  (`metadata/Matson_Supplement.pdf`, `metadata/Matson_Supplemental_Tables.zip`)
  are **not actually PDF/ZIP files** — `file` identifies both as HTML text,
  and `unzip` confirms the .zip is not a valid archive. These are failed
  downloads (likely a paywall/redirect page saved with the wrong extension)
  that were sitting in the project directory undetected until this check.
  Not obtained locally; not re-fetched from Science's site in this pass.
- **C4 (Lee 2022)**: confirmed via web search that the paper collected
  "standardized radiological response assessments (RECIST v1.1)" and that
  supplementary tables contain patient-level RECIST data. Not fetched or
  parsed in this pass (same crosswalk problem as C1 would apply — no
  patient-ID mapping exists locally for the 5-site PRIMM cohort either).

No PERMANOVA reruns were performed for this task since no cohort produced
usable strict labels. If this is pursued further, the concrete next
step is: pull Frankel Table 2 (open access, lowest-friction of the four),
resolve patient-ID-to-SRR mapping via NCBI BioSample attributes for those 39
samples, and rerun Task 2's C1 PERMANOVA with CR/PR vs PD (SD dropped) as a
pilot before attempting the other three.

---

## 4. Depth and classification rate

### 4a. Per-cohort median depth and classification rate

Script: `scripts/revision/review_task4_depth_classification.py`. Sources:
`results/fastp/**/*_fastp.json` (post-fastp read counts),
`results/kraken_reports/**/*_report.txt` (Kraken2 percentage columns for
"unclassified" [rank U] and "Homo" [genus, taxid 9605]). Output:
`results/revision/review/depth_classification_per_sample.tsv` (per-sample)
and `depth_classification_summary.tsv` (per-cohort).

| cohort | n | median read pairs (post-fastp) | median % classified excl. human |
|---|---|---|---|
| C1 | 39 | 41,080,015 | 52.46% |
| C2 | 40 | 25,263,955 | 63.26% |
| C3 | 39 | 35,742,645 | 75.00% |
| C4 | 160 | 20,614,361 | 55.12% |

278/283 samples had both files available; 5 C4 samples (`ERR6275667`,
`ERR6275672`, `ERR6275675`, `ERR6275676`, `ERR6279623`) had no fastp JSON on
disk (their intermediate files appear to have been cleaned up by the
project's storage-efficient download-process-delete pipeline pattern before
this check could read them) — median for C4 is over the 160 available, not
all 165.

### 4b. Rarefaction to a common depth

Script: `scripts/revision/review_task4b_rarefaction.py`. This required
re-parsing the raw Kraken2 reports for **integer** genus-level read counts
(`n_reads_clade`, column 2) — the existing `X_genus_raw.tsv` files used
everywhere else in this project store Kraken2's *percentage* column, which
cannot be rarefied (rarefaction is subsampling of integer counts). Same
genus name-cleaning and Homo-exclusion logic as `scripts/build_matrix_3cohort.py`,
for direct comparability. Per-sample total genus-read counts saved to
`results/revision/review/rarefaction_per_sample_totals.tsv`.

**Rarefaction depth used: 2,098,257 reads** — the *true minimum* total
genus-assigned read count across all 283 samples (lowest was `SRR5930514` at
2,098,257; the next-lowest, `ERR6279628` at ~2.1M-ish, was not a >3x outlier
below it, so no sample needed to be dropped as a low-depth outlier). All 283
samples had a total at or above this depth, so **0 samples were dropped** for
being below the rarefaction floor.

Rarefied count matrix: `results/revision/review/X_genus_raw_rarefied.tsv`.
Pooled genus PERMANOVA on the rarefied data (same settings — Aitchison
distance, CLR pseudocount=1e-6, 999 perms):

| dataset | factor | groups | R2 | E0 | delta_R2 | p |
|---|---|---|---|---|---|---|
| genus_rarefied_n283 | response | 2 | 0.003552 | 0.003546 | **+0.000005** | 0.336 |
| genus_rarefied_n283 | cohort | 4 | 0.081041 | 0.010638 | **+0.070403** | 0.001 |

Rarefaction changes essentially nothing: response ΔR² is even closer to
exactly zero than the unrarefied result (+0.000005 vs. the original
+0.00033), and cohort ΔR² remains large and significant (p=0.001). Depth
normalization was never the confound — the result was already depth-robust.

---

## 5. Table VII p-values (nested-LOOCV, n=118)

Source: `results/ml/nested_cv/permutation_test_nested.tsv`, columns
`n_perm_gte_obs`, `n_perms`, `empirical_pval` (N=100 permutations).

| model | n_perm_gte_obs | file's empirical_pval (n_ge/n_perms) | standard (n_ge+1)/(n_perms+1) |
|---|---|---|---|
| ElasticNet | 99 | 0.990 | 0.990 |
| RandomForest | 24 | 0.240 | 0.248 |
| Gradient Boosting (XGBoost) | 10 | 0.100 | 0.109 |

The file's own `empirical_pval` column is the simple fraction (n_ge/n_perms),
giving exactly 0.990/0.240/0.100 to 3 decimals. The standard
Phipson-and-Smyth-corrected form, (n_ge+1)/(n_perms+1) — used everywhere else
in this project's permutation tests, including `scripts/revision/lib.py`'s
`permanova()` — gives 0.990/0.248/0.109. Both are legitimate; recommend using
the corrected values (0.990/0.248/0.109) for consistency with the rest of the
paper's permutation tests, which all use that correction.

---

## 6. Runtime hardware

The machine used throughout this work — including the `batch_detector.py`
reruns that took 3.6 min (n=283, 217.2s) and 13.0 min (Duvallet, 781.1s) —
is this same persistent environment, queried directly during this review:

- **CPU**: Apple M3 Pro (arm64)
- **RAM**: 18 GB
- **OS**: macOS 26.2 (Darwin kernel 25.2.0)

---

## 7. Simulation figure (Fig 7) — v3

`results/revision/figures/fig7_simulation_grid_v3.{pdf,png}`. Same
`results/ml/simulation/grid_results.tsv` data — no numbers changed. Two
changes, plus one fix found while checking the PNG:

1. Axes relabeled "signal target $R^2$" / "cohort target $R^2$" (were
   "signal $f^2$" / "cohort $f^2$") — matches SUMMARY.md Section 13's finding
   that no R2-to-Cohen's-f2 conversion is ever applied in this codebase.
2. `percentile_norm` added back as a sixth panel, with the note
   "*percentile normalization assumes the same signal direction in every
   cohort" set as a figure-level footnote (keeps the panel title short: just
   "percentile norm*").
3. **Fix found on inspection**: adding percentile_norm's cells (which run up
   to delta_auc ≈ +0.95, a known generator artifact — see project notes)
   into the shared color-scale range washed out all five other panels to
   near-uniform white. Fixed by keeping the color scale's vmin/vmax computed
   from the original 5 methods only (unchanged from v2: ±0.0911);
   percentile_norm's own cells still render, they simply clip to the
   colorbar's end color (with `extend="both"` triangle caps added to the
   colorbar to signal this), which itself communicates "far outside the
   other methods' range" rather than hiding it.

Outline kept on the real operating point (signal target R2=0.007, cohort
target R2=0.08) in all six panels, including the new one. Real-operating-point
cell, now including percentile_norm:

| method | delta_auc |
|---|---|
| none | -0.0719 |
| mean_centering | -0.0590 |
| location_scale | -0.0364 |
| quantile_mapping | -0.0406 |
| cohort_covariate | -0.0658 |
| percentile_norm | -0.0746 |

At the real operating point (unlike elsewhere in its grid), percentile_norm
is *not* an outlier — it's in-range with the other five, all clearly
negative. Its extreme cells are confined to the higher signal/lower-signal
combinations elsewhere in the grid, consistent with its documented generator
artifact depending on signal-direction consistency, which higher signal_f2
cells simulate but the real (near-zero-signal) operating point does not.

PNG opened and checked after saving: no label collisions, no clipping, axes
and titles all legible.
