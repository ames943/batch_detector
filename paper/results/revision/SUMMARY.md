# BIBM camera-ready / arXiv revision — chance-calibrated re-analysis

All numbers below are traced to a specific output file. Nothing here overwrites
any existing result file — everything new lives under `results/revision/`.
Scripts: `scripts/revision/task*.py`. Updated tool: `scripts/batch_detector.py`
(v2.0.0, not yet published — see Section 9).

---

## 0. Quick facts

**a. Kraken2 version** (`kraken2 --version`):
```
Kraken version 2.17.1
Copyright 2013-2023, Derrick Wood (dwood@cs.jhu.edu)
```

**b. Phase-2 simulation Ledoit-Wolf covariance**: estimated from the **n=118**
matrix (top-75 most prevalent genera), not n=283. Source: `scripts/phase2_simulation.py`
lines 48-49, 59, 97 (`CLR_PATH = "results/ml/n118_3cohort/X_genus_clr.tsv"`,
`TOP_P = 75`, `lw = LedoitWolf().fit(clr_top)`). n=283 did not exist yet when
Phase 2 was run.

**c. n=283 classification (EN 0.517, RF 0.541, XGB 0.524)**: **nested** tuning
(inner 5-fold grid search per outer fold — `scripts/nested_cv_n283.py`
`inner_cv_select()`), not fixed parameters. Sensitivity/specificity from
`results/ml/n283_4cohort/nested_cv_n283_summary.tsv`:

| model | AUC | accuracy | sensitivity | specificity |
|---|---|---|---|---|
| ElasticNet | 0.5173 | 0.477 | 0.4531 | 0.4968 |
| RandomForest | 0.5409 | 0.5477 | 0.3906 | 0.6774 |
| XGBoost | 0.5235 | 0.5477 | 0.4141 | 0.6581 |

**d. Phase-4 meta-analysis per-cohort SEs**: computed as
`std(LOOCV fold coefficients) / sqrt(n)` (`scripts/phase4_meta_analysis.py`
`loocv_coefs()`, line `se_coefs = np.maximum(fold_mat.std(axis=0)/np.sqrt(n), MIN_SE)`),
where `n` = number of LOOCV folds = cohort size. **These SEs are far too small**,
for two compounding reasons: (1) LOOCV folds share n-1 of n training samples,
so they are highly correlated, not independent replicate estimates — their
empirical std already understates true between-sample coefficient variability;
(2) dividing by `sqrt(n)` again treats them *as if* independent, which they are
not, shrinking the SE a second time. This directly explains why Task 6 (proper
Hedges' g SEs) collapses mean I² from 99.5% → 14.4% and FDR-significant genera
from 4 → 0 (same data, correct uncertainty).

---

## 1. Chance-calibrated PERMANOVA

File: `task1_chance_calibrated_permanova.tsv`. E0 = (g-1)/(n-1). All PERMANOVA
runs use 999 permutations, Aitchison distance (Euclidean on CLR) for microbiome/
Duvallet, Euclidean-on-standardized-features for tumor (matches
`scripts/tumor_batch.py`'s original method — TMB/mutation flags are not
compositional counts, so CLR does not apply there).

| dataset | level | n | factor | g | R2 | E0 | ΔR2 | p |
|---|---|---|---|---|---|---|---|---|
| genus_n39 | genus | 39 | response | 2 | 0.02675 | 0.02632 | **+0.00043** | 0.419 |
| genus_n79 | genus | 79 | response | 2 | 0.01100 | 0.01282 | **-0.00182** | 0.736 |
| genus_n79 | genus | 79 | cohort | 2 | 0.08491 | 0.01282 | +0.07209 | 0.001 |
| genus_n118 | genus | 118 | response | 2 | 0.00679 | 0.00855 | **-0.00176** | 0.874 |
| genus_n118 | genus | 118 | cohort | 3 | 0.07677 | 0.01709 | +0.05968 | 0.001 |
| genus_n283 | genus | 283 | response | 2 | 0.00388 | 0.00355 | **+0.00033** | 0.284 |
| genus_n283 | genus | 283 | cohort | 4 | 0.10563 | 0.01064 | +0.09499 | 0.001 |
| phylum_n118 | phylum | 118 | response | 2 | 0.00660 | 0.00855 | -0.00195 | 0.629 |
| phylum_n118 | phylum | 118 | cohort | 3 | 0.05691 | 0.01709 | +0.03982 | 0.001 |
| species_n118 | species | 118 | response | 2 | 0.00804 | 0.00855 | -0.00051 | 0.595 |
| species_n118 | species | 118 | cohort | 3 | 0.07402 | 0.01709 | +0.05693 | 0.001 |
| tumor_n223 | gene_flags | 223 | response | 2 | 0.01200 | 0.00450 | +0.00749 | 0.003 |
| tumor_n223 | gene_flags | 223 | cohort | 3 | 0.14273 | 0.00901 | +0.13372 | 0.001 |
| duvallet_n570 | genus | 570 | response | 2 | 0.00618 | 0.00176 | +0.00443 | 0.001 |
| duvallet_n570 | genus | 570 | cohort | 4 | 0.24456 | 0.00527 | +0.23929 | 0.001 |

**Headline reframing**: at every one of our own microbiome n's (39/79/118/283),
response ΔR² is essentially zero — **and actually negative at n=79 and n=118**
(observed R² falls *below* its own chance expectation). Cohort ΔR² stays large
(6-9.5 percentage points) and significant (p=0.001) throughout. Tumor and
Duvallet both show small but real, significant positive response ΔR²
(+0.0075 and +0.0044) — useful contrast cases where "signal exists" is
genuinely true, and Section 9 shows the tool still (correctly) says NO-GO for
Duvallet because the batch effect swamps it anyway.

**Sanity check** (`task1_sanity_check_E0.tsv`) — mean permuted R² vs. E0, genus n=118,
999 label-shuffles:

| factor | g | mean permuted R2 | E0 |
|---|---|---|---|
| response | 2 | 0.008651 | 0.008547 |
| cohort | 3 | 0.017068 | 0.017094 |

Confirms the formula: mean permuted R² and E0 agree to the 3rd decimal.

**rho = R2_cohort/R2_response permutation null** (`task1_rho_permutation_null.tsv`,
999 independent shuffles of each factor):

| dataset | observed rho | null mean | null median | percentile of observed |
|---|---|---|---|---|
| genus_n79 | 7.72x | 1.06 | 1.01 | 100.0 |
| genus_n118 | 11.31x | 2.09 | 2.04 | 100.0 |
| genus_n283 | 27.25x | 3.12 | 3.07 | 100.0 |
| phylum_n118 | 8.62x | 2.44 | 2.00 | 98.9 |
| species_n118 | 9.21x | 2.08 | 2.07 | 100.0 |
| tumor_n223 | 11.90x | 2.20 | 2.11 | 100.0 |
| duvallet_n570 | 39.55x | 3.18 | 3.06 | 100.0 |

Every observed rho sits at or above the ~99th-100th percentile of what
independent random shuffles of the same two factors would produce — the
batch:signal dominance is not an artifact of small-sample noise in the ratio.

---

## 2. Joint model with restricted permutation

File: `task2_joint_model_margin_permanova.tsv`. No R/vegan available in this
environment (`Rscript` not found) — implemented from scratch in Python
(`scripts/revision/lib.py::margin_permanova_two_factor`), a Type-III/margin
two-factor PERMANOVA via Gower-centered hat-matrix partitioning
(McArdle & Anderson 2001), equivalent to
`adonis2(dist ~ cohort + response, by="margin", permutations=how(nperm=999, blocks=cohort))`.
Response permutations are restricted within cohort strata; cohort permutations
are unrestricted (cohort is the blocking factor and cannot be tested under its
own restriction — this mirrors the plain single-factor batch PERMANOVA in
Section 1).

| dataset | n | R2(response\|cohort) | p | R2(cohort\|response) | p |
|---|---|---|---|---|---|
| genus_n79 | 79 | 0.01112 | 0.590 | 0.08503 | 0.001 |
| genus_n118 | 118 | 0.00582 | 0.959 | 0.07581 | 0.001 |
| genus_n283 | 283 | 0.00346 | 0.312 | 0.10521 | 0.001 |
| tumor_n223 | 223 | 0.00713 | 0.025 | 0.13786 | 0.001 |
| duvallet_n570 | 570 | 0.00620 | 0.001 | 0.24457 | 0.001 |

Margin R² values are close to the plain single-factor R² from Section 1 (as
expected — cohort and response are only weakly associated), and the margin
response p-values track the plain response p-values closely, confirming the
restricted-permutation implementation is behaving sensibly.

---

## 3. Zero-handling sensitivity

Files: `task3_zero_handling_sensitivity.tsv`, `task3_summary_range.tsv`. Redid
genus n=39/79/118/283 three ways: pseudocount ε=0.5, ε=1.0, multiplicative
replacement (implemented directly, no skbio dependency in `batch_detector.py`;
`scripts/revision/lib.py` uses skbio when available with an identical manual
fallback).

**Range of ΔR² across all three zero-handling methods:**

| dataset | factor | min ΔR2 | max ΔR2 |
|---|---|---|---|
| genus_n39 | response | -0.0061 | -0.0019 |
| genus_n79 | response | -0.0028 | -0.0025 |
| genus_n79 | cohort | +0.0856 | +0.0887 |
| genus_n118 | response | -0.0015 | +0.0010 |
| genus_n118 | cohort | +0.0772 | +0.1176 |
| genus_n283 | response | -0.0009 | +0.0003 |
| genus_n283 | cohort | +0.1070 | +0.1844 |

**Conclusion**: cohort ΔR² stays large (7.7-18.4 percentage points) and highly
significant (p=0.001 in every cell) regardless of zero-handling choice.
Response ΔR² stays within ±0.001-0.006 of zero everywhere — never
distinguishable from chance. The paper's central finding is not a pseudocount
artifact.

---

## 4. C4 alone

File: `task4_c4_alone_permanova.tsv`. Lee 2022, n=165, 5 recruiting sites
(PRIMM-UK=55, PRIMM-NL=55, Manchester=25, Leeds=18, Barcelona=12).

| factor | g | R2 | E0 | ΔR2 | p |
|---|---|---|---|---|---|
| response | 2 | 0.006631 | 0.006098 | +0.000533 | 0.277 |
| site | 5 | 0.094009 | 0.024390 | +0.069619 | 0.001 |

Even within the single largest cohort, recruiting site alone produces a
significant, large batch-like effect (ΔR²=+0.070, p=0.001) while response
signal is indistinguishable from chance (ΔR²=+0.0005, p=0.277).

---

## 5. Repeated stratified 5-fold CV (new primary protocol)

Files: `task5_{c1,n118,n283,c4}_fold_aucs.tsv`, `task5_*_permutation_summary.tsv`,
`task5_*_permutation_aucs.tsv`, `task5_*_consensus_params.json`.
10 repeats x stratified 5-fold, inner 5-fold hyperparameter tuning
(ElasticNet/RF/XGBoost grids identical to `scripts/nested_cv_n283.py`),
same leak-free per-fold prevalence->variance->top-100-|point-biserial-r|
feature selection, per-fold mean-centering batch correction for the two
multi-cohort datasets (n118, n283). C1 and C4 are single-cohort, no
correction applied (matches existing `nested_cv_c4.py`).

**Permutation-test methodology note** (documented per the task instructions):
re-tuning the full inner grid for all 50 outer folds under every one of
N=100 permuted label draws would multiply cost ~100x for negligible benefit.
Following the precedent already set by this codebase's own
`scripts/nested_cv_n283.py` (which permutation-tests with **fixed consensus
hyperparameters**, not full re-tuning, under every permutation), the null
here uses the same approach: the most-frequent best hyperparameters from the
50 observed outer folds, with feature selection and batch correction still
freshly re-done per fold per permutation (no leakage introduced).

**Observed: mean AUC ± SD across the 50 outer folds**, and **permutation test (N=100, fixed consensus params)**:

| dataset | n | model | mean AUC | SD | perm mean | perm SD | p |
|---|---|---|---|---|---|---|---|
| C1 | 39 | ElasticNet | 0.390 | 0.180 | 0.485 | 0.131 | 0.752 |
| C1 | 39 | RandomForest | 0.431 | 0.229 | 0.495 | 0.126 | 0.653 |
| C1 | 39 | XGBoost | 0.392 | 0.215 | 0.491 | 0.129 | 0.762 |
| n=118 | 118 | ElasticNet | 0.376 | 0.090 | 0.491 | 0.068 | 0.950 |
| n=118 | 118 | RandomForest | **0.546** | 0.107 | 0.488 | 0.072 | 0.248 |
| n=118 | 118 | XGBoost | 0.541 | 0.117 | 0.488 | 0.070 | 0.248 |
| n=283 | 283 | ElasticNet | 0.508 | 0.063 | 0.502 | 0.041 | 0.455 |
| n=283 | 283 | RandomForest | 0.537 | 0.060 | 0.502 | 0.043 | 0.228 |
| n=283 | 283 | XGBoost | 0.533 | 0.053 | 0.504 | 0.040 | 0.198 |
| C4 | 165 | **ElasticNet** | **0.571** | 0.099 | 0.504 | 0.051 | **0.109** |
| C4 | 165 | RandomForest | 0.502 | 0.104 | 0.502 | 0.059 | 0.485 |
| C4 | 165 | XGBoost | 0.531 | 0.120 | 0.491 | 0.056 | 0.267 |

Total wall time: C1 15.2 min, n=118 17.1 min, n=283 19.2 min, C4 16.1 min
(observed run + 100-permutation null, each; consensus hyperparameters and
per-fold AUCs in `task5_*_consensus_params.json` / `task5_*_fold_aucs.tsv`).

**All 12 model x dataset combinations remain non-significant at α=0.05** under
the new repeated-CV protocol — consistent with, and slightly more favorable
in absolute AUC than, the original LOOCV numbers (e.g. n=118 RF 0.546 vs.
LOOCV's 0.542; C4 ElasticNet 0.571 vs. the original single LOOCV run's 0.570).
C4 ElasticNet at p=0.109 remains the closest-to-significant result across the
whole paper, same as under LOOCV — the new protocol does not change which
result is "least null," it just adds a second, independent CV design that
confirms the same ranking with tighter per-fold SDs (5 samples per outer fold
vs. 1 in LOOCV, so AUC is no longer computed on n=1 test points).

---

## 6. Meta-analysis redo (Hedges' g, proper SEs)

Files: `task6_meta_analysis_hedges_g.tsv`, `task6_per_cohort_hedges_g.tsv`,
`task6_summary.txt`. Per Quick Fact 0d, the original Phase-4 SEs were
fold-based and far too small. Redone with per-cohort **Hedges' g** (bias-
corrected standardized mean difference, R vs NR, on CLR abundance) and its
closed-form SE, combined via DerSimonian-Laird across all 4 cohorts
(prevalence filter: genus present in >=10% of samples within a cohort to be
tested in that cohort, same threshold used throughout the original pipeline).

| metric | value |
|---|---|
| Genera tested (>=2 cohorts) | 125 |
| FDR q < 0.05 | **0** (was 4) |
| Mean I² | **14.4%** (was 99.5%) |

**This changes the paper's meta-analysis narrative substantially.** With
correct SEs, heterogeneity across cohorts is actually *low* (14.4%, not
99.5%) — the earlier "near-universal extreme heterogeneity" conclusion was
itself an artifact of underestimated SEs inflating Cochran's Q. The corrected,
honest result: 0 genera reach FDR significance, full stop — a cleaner and more
defensible null than the original 4-genera/I²=99.5% result, and it no longer
needs the "hypothesis-generating only" hedge the original required.
**Recommend replacing the Phase-4 meta-analysis paragraph in the paper with
this corrected version.**

---

## 7. Cross-cohort transfer breakdown

File: `task7_cross_cohort_full_table.tsv` — all 18 Phase 0.5 LOCO+S2S rows
(source: `results/ml/cross_cohort/holdout_summary.tsv`).

| split | train | test | model | AUC | perm p |
|---|---|---|---|---|---|
| LOCO | 1+2 | 3 | ElasticNet | 0.3667 | 0.921 |
| LOCO | 1+2 | 3 | RandomForest | 0.3319 | 0.953 |
| LOCO | 1+3 | 2 | ElasticNet | 0.4107 | 0.825 |
| LOCO | 1+3 | 2 | RandomForest | 0.3920 | 0.867 |
| LOCO | 2+3 | 1 | ElasticNet | 0.3083 | 0.982 |
| LOCO | 2+3 | 1 | RandomForest | 0.2875 | 0.984 |
| S2S | 1 | 2 | ElasticNet | 0.4000 | 0.842 |
| S2S | 1 | 2 | RandomForest | 0.3373 | 0.944 |
| S2S | 2 | 1 | ElasticNet | 0.5500 | 0.338 |
| S2S | 2 | 1 | RandomForest | 0.3722 | 0.911 |
| S2S | 1 | 3 | ElasticNet | 0.3944 | 0.841 |
| S2S | 1 | 3 | RandomForest | 0.2986 | 0.984 |
| S2S | 3 | 1 | ElasticNet | 0.3750 | 0.898 |
| S2S | 3 | 1 | RandomForest | 0.3625 | 0.925 |
| S2S | 2 | 3 | ElasticNet | 0.4250 | 0.783 |
| S2S | 2 | 3 | RandomForest | 0.5014 | 0.477 |
| S2S | 3 | 2 | ElasticNet | 0.3253 | 0.968 |
| S2S | 3 | 2 | RandomForest | 0.4413 | 0.757 |

**Which held-out cohort accounts for most of the below-chance results?**
Below-chance (AUC<0.5) counts by held-out (test) cohort: cohort 1 = 5/6,
cohort 2 = **6/6**, cohort 3 = 5/6. Mean AUC by held-out cohort: cohort 1 =
0.376, cohort 2 = 0.384, cohort 3 = 0.386 — essentially identical. **No single
cohort dominates the below-chance pattern**; it is broadly distributed. The
only mild standout is cohort 2 (100% of its splits below chance) — notable
given cohort 2's identity (Peng 2020, GI cancer, not melanoma; see project
notes), consistent with it being the most biologically distinct cohort, but
the effect is small (5/6 vs 6/6) and does not change the overall conclusion
that transfer failure is general, not concentrated in one cohort.

---

## 8. Power over hypothetical effects

Files: `task8_power_curve.tsv`, `task8_c4_upper_bound.tsv`. Dirichlet
generator (Gamma/sum-normalize identity) on a single real cohort (C1,
Frankel n=39, top-30 most prevalent genera for tractability up to n=1000),
with a phase-2-style binary-search calibration of an injected CLR-space
signal shift per (n, target ΔR²) so the expected observed R² ≈
E0(2,n) + target ΔR².

**Power(α=0.05) for true ΔR2 at various n** (`task8_power_curve.tsv`):

| n | ΔR2=0.0025 | ΔR2=0.005 | ΔR2=0.01 | ΔR2=0.02 |
|---|---|---|---|---|
| 40 | 0.09 | 0.07 | 0.16 | 0.69 |
| 80 | 0.13 | 0.26 | 0.52 | 0.99 |
| 165 | 0.12 | 0.54 | 0.99 | 1.00 |
| 300 | 0.55 | 0.97 | 1.00 | 1.00 |
| 500 | 0.77 | 1.00 | 1.00 | 1.00 |
| 1000 | 1.00 | 1.00 | 1.00 | 1.00 |

(full table with calibrated injected-δ values in the TSV. The n=40 dip from
0.09→0.07 between ΔR²=0.0025 and 0.005 is Monte Carlo noise — only 100 reps
per cell at a small, noisy effect size; not a real non-monotonicity.)

**95% upper bound on C4's own response ΔR²** (Dirichlet bootstrap of the
*observed* Lee 2022 data — no injected effect, pure resampling uncertainty of
the effect actually measured, 500 reps):

| quantity | value |
|---|---|
| n (C4) | 165 |
| E0 | 0.006098 |
| observed ΔR2 (Task 4) | 0.000533 |
| bootstrap median ΔR2 | 0.002071 |
| **bootstrap 95th-percentile ΔR2 (upper bound)** | **0.007806** |
| n for 80% power at that bound | **~131** |

Even under a generous 95% upper bound on the true C4 response effect (using
the whole bootstrap distribution, not just the point estimate), 80% power
would need only ~131 same-protocol patients — well within reach of C4's own
n=165. The fact that C4 (n=165) still lands at p=0.277 (Task 4) despite this
is consistent with the true effect being close to its point estimate (near
zero), not with the study being underpowered for a real effect of this size.

---

## 9. batch_detector.py update (v2.0.0 — NOT yet published)

All six requested changes implemented in `scripts/batch_detector.py`:

1. **R2/E0/ΔR2 for both factors** — new `chance_expectation_r2()`, logged and
   in the JSON report (`response_E0`, `response_delta_R2`, `batch_E0`,
   `batch_delta_R2`).
2. **Joint model with within-cohort-permuted response** —
   `margin_permanova_two_factor()` (same Python implementation as Section 2,
   duplicated inline to keep the tool dependency-free, not imported from
   `scripts/revision/lib.py`).
3. **rho permutation-null percentile** — `rho_permutation_null()`.
4. **Power as a curve over hypothetical effect sizes** —
   `power_over_effect_grid()` is now the *default* (`--legacy-power` restores
   the old observed-effect-resampling behavior); output:
   `batch_detector_power_grid.tsv`.
5. **Zero-handling option** — `--zero-handling {pseudocount,multiplicative_replacement}`
   (`clr_multiplicative_replacement()`, no extra dependency).
6. **New chance-calibrated recommendation rule** — `_pooling_recommendation_v2()`
   (default; `--legacy-recommendation` restores the old rule). Thresholds:
   - **NO-GO** if response ΔR² <= 0, OR response p >= 0.05, OR rho > 4x, OR
     rho's permutation-null percentile >= 95.
   - **GO** if response ΔR² > 0, p < 0.05, and rho <= 0.25x.
   - **CAUTION** if response ΔR² > 0, p < 0.05, and 0.25x < rho <= 4x.
   The 0.25x/4x thresholds come directly from the paper's own Phase-2
   calibrated simulation, which found batch correction only reliably
   recovers a no-batch ceiling AUC when the biological signal is at least
   ~4x the batch effect (i.e. rho <= 0.25x); 4x is the mirror-image threshold.
   Full derivation is in the function's docstring.

**Rerun on n=283 (our data) and Duvallet, with runtime:**

| | our n=283 | Duvallet n=570 |
|---|---|---|
| response R² / E0 / ΔR² | 0.00395 / 0.00355 / +0.00040 | 0.00618 / 0.00176 / +0.00443 |
| response p | 0.266 | 0.001 |
| batch R² / E0 / ΔR² | 0.10098 / 0.01064 / +0.09034 | 0.24456 / 0.00527 / +0.23929 |
| batch p | 0.001 | 0.001 |
| joint margin R²(resp\|cohort) / p | 0.00335 / 0.384 | 0.00620 / 0.001 |
| rho (batch/response) | 25.6x | 39.55x |
| rho null percentile | 100.0 | 100.0 |
| **recommendation** | **NO-GO** | **NO-GO** |
| reason | response ΔR² positive but p=0.266 (not significant) — can't distinguish from noise | response IS significant, but rho sits at the 100th percentile of its own independent-shuffle null — structural batch dominance overrides a real signal |
| min n for 80% power (by target ΔR²) | 0.0025→566, 0.005→283, 0.01→142, 0.02→142 | 0.0025→570, 0.005→285, 0.01→285, 0.02→285 |
| **wall-clock runtime** (`/usr/bin/time -p`) | **217.2s real** (203.4s user, 5.0s sys) | **781.1s real** (748.1s user, 10.1s sys) |

Full JSON reports + PDF figures: `task9_n283/batch_detector_report.json` (+`_figure.pdf`),
`task9_duvallet/batch_detector_report.json` (+`_figure.pdf`).

This pair is the tool's headline demonstration: on our own data the call is
"can't tell if there's signal" (response not significant); on Duvallet the
call is "there IS real signal (p=0.001) but pooling is risky" (rho
verdict overrides). The v2 rule correctly produces two different *reasons*
for the same NO-GO verdict, which the v1 raw-ratio rule could not
distinguish.

---

## 10. Figure 1 (regenerated)

`dose_response.pdf` / `.png` — regenerated from `task1_chance_calibrated_permanova.tsv`
via `scripts/revision/task10_figure_dose_response.py`. Single-column width
(3.4 x 2.0 in), vector PDF, embedded fonts (`pdf.fonttype=42`).

- **Panel A**: observed response R² (filled circles) and cohort R² (open
  squares) at n=39/79/118/283, ordinal x-axis, dotted chance curves
  1/(n-1) [response] and (k-1)/(n-1) [cohort]; all 4 response and 3 cohort
  values labeled.
- **Panel B**: ΔR² for both factors, horizontal zero line; same labeling.

The regenerated figure visually makes the Section 1 finding unmissable: the
response curve in Panel A sits almost exactly on top of its own chance line,
while Panel B shows the response series pinned to ~0% throughout and the
cohort series sitting at +6 to +9.5 percentage points above chance at every
multi-cohort n.

---

## 11. Power check v2 — a bug in the C4 power claim, found and fixed

Prompted by a reviewer-style sanity check that the original claim ("95% upper
bound ΔR²=0.0078 → 80% power at n≈131") didn't sit consistently on the
plotted power curves. It didn't, because it wasn't computed the way it was
described. Two bugs, both now fixed; new script
`scripts/revision/task8b_c4_power_v2.py`; new data
`results/revision/task8_c4_power_v2.tsv`; new figure
`results/revision/figures/fig4_power_curves_v2.pdf` (see `figures/FIGURES.md`
for the full figure-level writeup — this section is the numeric summary).

**Bug 1 — wrong curve, silently substituted.** In
`scripts/revision/task8_power_hypothetical_effects.py` (lines 213-224), the
lookup for "n at 80% power for ΔR²=0.0078" picked the *nearest* of the four
calibrated target deltas `[0.0025, 0.005, 0.01, 0.02]` to 0.007806 — which is
**0.01**, not 0.0078 (0.0078 was never itself a calibrated target in that
run). It then interpolated n on the 0.01 curve (power 0.52 at n=80, 0.99 at
n=165) to get n=130.6, and reported that as "at ΔR²=0.0078." This is exactly
why the star sat on the ΔR²=0.01 curve in the original figure. The
substitution was recorded in a JSON field
(`closest_calibrated_delta_used_for_lookup: 0.01`) but never surfaced in the
headline number or the figure.

**Bug 2 — mismatched feature space.** That power curve was built from
**Cohort 1's top-30 genera** (`task8_power_hypothetical_effects.py` line
149-155), while the 95% upper bound itself was correctly computed from
**C4's own top-30 genera**. The lookup crossed feature spaces without saying so.

**The fix.** `task8b_c4_power_v2.py` recomputes the power curve entirely in
C4's own feature space. The aim was to match "the same prevalence
filtering as the Task 4 PERMANOVA": reading `scripts/revision/task4_c4_alone.py`
confirms Task 4 applied **no prevalence filter at all** — it calls
`clr_transform(raw.values)` directly on every genus column. So "the same
filtering" is no filtering, and this rerun uses C4's full, unfiltered genus
set (p=3005) for the Aitchison distance, exactly matching Task 4. (The 10
signal-injection columns are chosen as the most-prevalent genera purely so
the injected shift lands on a non-degenerate column — the distance
computation itself still uses every column.) Same injected-shift
binary-search calibration as Task 8. Targets now include 0.0078 directly —
no more nearest-neighbor substitution: ΔR² in {0.0025, 0.005, 0.0078, 0.01},
n in {80, 131, 165, 250, 400}, 200 reps/cell, α=0.05.

**Corrected power grid** (`task8_c4_power_v2.tsv`):

| n | ΔR²=0.0025 | ΔR²=0.005 | ΔR²=0.0078 | ΔR²=0.01 |
|---|---|---|---|---|
| 80 | 0.145 | 0.585 | 0.855 | 0.885 |
| 131 | 0.385 | 0.875 | 0.995 | 1.000 |
| 165 | 0.490 | 0.935 | **0.995** | 1.000 |
| 250 | 0.835 | 1.000 | 1.000 | 1.000 |
| 400 | 0.980 | 1.000 | 1.000 | 1.000 |

**Answers to the four questions asked:**
1. How was n=130.6 computed? Interpolated on the **ΔR²=0.01** curve (the
   nearest calibrated target to 0.0078), not an interpolation between two
   deltas and not a different generator — see Bug 1 above.
   `scripts/revision/task8_power_hypothetical_effects.py` lines 213-224.
2. Power at n=165 for ΔR²=0.0078, in C4's own feature space: **0.995**.
3. n reaching 80% power at ΔR²=0.0078: already reached at **n=80**, the
   smallest n tested (power=0.855 there).
4. Minimum detectable ΔR² (80% power) at n=165: **≈0.0042** — about half of
   the 0.0078 upper bound.

**Does "C4 was large enough to detect an effect at its own 95% upper bound"
still hold?** Yes — more strongly than originally claimed, not more weakly.
C4 had ~99.5% power for its own upper-bound effect (not the ~80% v1 implied),
and its 80%-power detection threshold (~0.0042) is about half the upper bound
itself. C4 was not underpowered for a real effect anywhere near what its own
data could plausibly contain — the null result (Task 4, p=0.277) is better
explained by the true effect sitting close to its point estimate (~0) than by
insufficient sample size. This is a *stronger* form of the paper's argument
than the (buggy) original version.

**Open, unresolved observation:** at matched (n, target ΔR²), C4's own
full-genus-space power is far higher than Cohort 1's top-30-genera power was
for the same nominal effect (e.g. at n=165, ΔR²=0.005: 0.935 here vs. 0.54 in
the original Task 8 curve). A plausible mechanism is that a higher-dimensional
base composition (p=3005 vs. p=30) reduces the sampling variance of the
Aitchison-distance PERMANOVA statistic across resamples once the injected
shift is recalibrated to hit the same mean R2 — but this has not been
independently verified with a dedicated ablation (e.g. rerunning C4's own
data restricted to its own top-30 genera). Flagged for anyone extending this
work, not asserted as settled.

---

## 12. Figure fixes — new numbers introduced

Full writeups are in `results/revision/figures/FIGURES.md`; only the numbers
that changed (not pure layout fixes) are repeated here, each with its source.

**Fig 7 corrected real-operating-point cell** — the outline was at the wrong
cell (signal f²=0.027, cohort f²=0.08); the actual real operating point is
signal f²=0.007 (calibrated to injected shift ≈0, i.e. the null model, chosen
to match the real n=118 pooled response R²— confirmed from
`results/ml/simulation/calibration_deltas.tsv`: target_r2=0.007 →
delta=4.657e-09, verified_r2=0.00851). Corrected per-method ΔAUC at (signal
f²=0.007, cohort f²=0.08), source `results/ml/simulation/grid_results.tsv`:

| method | delta_AUC |
|---|---|
| none | -0.0719 |
| mean_centering | -0.0590 |
| location_scale | -0.0364 |
| quantile_mapping | -0.0406 |
| cohort_covariate | -0.0658 |

All five are more negative than the previously-reported (wrong-cell) values —
correction looks even less helpful at the real operating point than the
mislabeled cell suggested.

---

## 13. Final checks before lock

Full writeup: `results/revision/final/FINAL_CHECKS.md`. Summary of every new
number, each traced to its source file.

**Fig 4 layout fix (v3)** — `results/revision/final/fig4_power_curves_v3.pdf`.
Same data as v2 (`task8_c4_power_v2.tsv`); annotation moved to empty
lower-right space with a leader line to the star, "C4 (n=165)" label's stray
leader line removed. No new numbers.

**batch_detector ρ threshold rationale (0.25x/4x) — NOT supported by the
grid.** Every ΔAUC>0 cell (nc=3, npc=40) from `results/ml/simulation/grid_results.tsv`,
percentile_norm excluded, spans signal:cohort f² ratios from 0.18 to infinity
with no clustering near 4 or 0.25 — e.g. at signal f²=0.10, four of five
methods post positive ΔAUC even at cohort f²=0.15 (ratio 0.67, i.e. signal
*smaller* than cohort effect). The only real pattern: positive cells occur
almost exclusively at signal f²=0.10 (the grid's highest absolute level),
regardless of ratio. The "~4x" in `batch_detector.py`'s
`_pooling_recommendation_v2()` docstring is very likely 0.10/0.027≈3.7 — a
signal-to-signal comparison across the grid's own two highest levels, not a
signal-to-cohort ratio. No R²-to-Cohen's-f² conversion exists anywhere in
`scripts/phase2_simulation.py` — "f²" there is a direct relabeling of
PERMANOVA R² (confirmed via `calibrate_delta(target_r2, ...)` and inline
comments), so `rho` and the grid's f² columns compare directly without
further conversion. Recommended paper wording (tool behavior unchanged, per
instructions): *"The tool's ρ thresholds (0.25x, 4x) are conservative
heuristic cutoffs; the calibrated simulation did not show a clean ρ-based
threshold for when correction helps — positive ΔAUC occurred across
signal:cohort ratios from well below 1 to infinity, with the only consistent
pattern being the grid's highest absolute signal level (f²=0.10)."*

**Dimensionality vs. cohort identity (optional check)** —
`results/revision/final/task8c_c4_power_top30.tsv` (C4, top-30 genera, n=165,
200 reps): power=0.625 at ΔR²=0.005, power=0.900 at ΔR²=0.0078 — compare to
C4's full-genus-space power (0.935 / 0.995, `task8_c4_power_v2.tsv`) and
Cohort 1's top-30 power (0.54 / ~0.79 interpolated, `task8_power_curve.tsv`).
**Dimensionality is the dominant driver** of the original power gap (shrinking
C4's own feature space from p=3005 to p=30 reproduces most of the gap), but a
smaller real cohort-specific component also exists (C4-top30 still outperforms
C1-top30 at matched n and target ΔR²).

**Number verification (Section 4 of FINAL_CHECKS.md)** — 10 of 12 claims
matched exactly against their source files (Dirichlet bootstrap CI, PERMDISP
response/cohort at n=39-283 and at phylum/genus/species, C1 EN LOOCV AUC,
sign-inversion test, tested-on-Riaz range, median TMB, percentile-norm R²
reduction, full 6-method batch-correction table — all confirmed, files listed
in FINAL_CHECKS.md). **Two do not hold as currently worded:**
- Nested-CV gradient-boosting "82-85% of folds" (`results/ml/nested_cv/nested_cv_best_params.tsv`):
  the full (depth=2, n_est=100, lr=0.05) combination occurs in only 34/118=28.8%
  of folds. 82 and 85 are real but are raw **fold counts** for two different
  single-hyperparameter marginals (lr=0.05: 82/118 folds=69.5%; n_est=100:
  85/118 folds=72.0%), not a joint percentage.
- "Hugo→Liu transfer AUC 0.58-0.62" (`results/ml/tumor/tumor_loso_results.tsv`,
  cross-checked against `scripts/tumor_ml.py` lines 68-70): the numbers are
  right but the direction is backwards — `held_out=hugo2016` means trained on
  Liu2019+Riaz2017 pooled, **tested on Hugo2016**, not "Hugo→Liu." No
  single-cohort-to-single-cohort tumor result exists anywhere in the pipeline.
- "Top Liu features include TTN, MUC16, ADGRV1" conflates two different
  rankings: TTN/MUC16 top raw mutation *frequency* (`liu2019_features.tsv`,
  computed directly: 63.9%/59.7%), but XGBoost feature *importance*
  (`logs/liu_single_cohort.log`) top-10 is RAC1, DNAH5, OBSCN, ADGRV1, NEB,
  PKHD1L1, CSMD2, FAT4, MAP2K1, KIT — TTN/MUC16 aren't in it, and RAC1 is
  ranked #1, not merely "the only canonical driver buried in the list."

**References (Section 5)** — all 8 checked against CrossRef (one via web
search after a CrossRef rate-limit) and all match exactly on
title/volume/issue/pages/DOI. One update: Orletskaia &
Olekhnovich's bioRxiv preprint (doi:10.1101/2025.05.07.652660) has since been
peer-reviewed and published in *Computational and Structural Biotechnology
Journal* (doi:10.34133/csbj.0065); that version can be cited instead of,
or alongside, the preprint.

---

## 14. Review-response tasks

Full writeup: `results/revision/review/REVIEW_TASKS.md`. Summary of every new
number, each traced to its source file.

**1. Cohort inclusion accounting** — `results/revision/review/cohort_screening.tsv`,
source `metadata/response_sheet.xlsx` (624 rows, 11 studies, sums exactly to
the compilation total). Of 11 studies, this project processed 4 (Frankel
2017/C1, Peng 2020/C2, Matson 2018/C3, Lee 2022/C4, n=283); the reason the
other 7 (Gopalakrishnan 2018, Gunjur 2024, Heshiki 2020, Liu 2022, McCulloch
2022, Spencer 2021, Tsakmaklis 2023; 341 samples) weren't included is **not
documented anywhere in project files** — a project-wide search found zero
references to any of the 7 by name. **C2's accession (PRJNA615114) was never
recorded in any project file** (only discovered missing during this review);
confirmed correct via a live ENA API lookup on `SRR11413606` during this review,
but the paper should note this provenance gap. Two small unreconciled
sample-count discrepancies: Matson_2019 is 38 in the compilation vs. 39 in
this project's C3; Lee_2022 is 164 vs. 165 in C4.

**2. Per-cohort PERMANOVA** — `results/revision/review/percohort_permanova_C1_C4.tsv`.
All four cohorts individually null: C1 ΔR²=+0.00043 (p=0.419), C2
ΔR²=+0.0023 (p=0.275), C3 ΔR²=-0.0025 (p=0.635), C4 ΔR²=+0.00053 (p=0.277).

**3. Strict response definition (CR/PR vs PD)** — **not obtained for any
cohort.** C1 (Frankel Table 2, confirmed to have per-patient RECIST via PMC)
blocked by a missing patient-ID-to-SRR crosswalk; C2 blocked by an AACR
paywall (previously documented); C3 blocked by two local supplement files
(`metadata/Matson_Supplement.pdf`, `Matson_Supplemental_Tables.zip`) that
turned out to be failed downloads (HTML error pages saved with the wrong
extension, confirmed via `file`/`unzip`) — a previously-undetected data
problem this review surfaced; C4 (Lee 2022, confirmed to have RECIST v1.1 in
its supplement) not attempted, same crosswalk problem as C1 would apply. No
PERMANOVA reruns were possible without usable labels.

**4a. Depth and classification rate** — `results/revision/review/depth_classification_summary.tsv`
(278/283 samples; 5 C4 samples missing intermediate files). Median read
pairs post-fastp: C1=41.1M, C2=25.3M, C3=35.7M, C4=20.6M. Median %
Kraken2-classified excluding human: C1=52.5%, C2=63.3%, C3=75.0%, C4=55.1%.

**4b. Rarefaction** — required re-parsing raw Kraken2 reports for integer
counts (existing `X_genus_raw.tsv` files store *percentages*, not counts, so
cannot be rarefied directly — confirmed by reading `scripts/build_matrix_3cohort.py`).
Rarefied all 283 samples to 2,098,257 reads (the true per-sample minimum; 0
samples dropped). Source: `results/revision/review/rarefied_permanova.tsv`.
Response: R²=0.003552, E0=0.003546, ΔR²=**+0.000005**, p=0.336. Cohort:
R²=0.081041, E0=0.010638, ΔR²=+0.070403, p=0.001. Rarefaction changes
essentially nothing — response ΔR² lands even closer to exactly zero than
the unrarefied result; depth was never the confound.

**5. Table VII p-values** — `results/ml/nested_cv/permutation_test_nested.tsv`.
File's own values (n_ge/n_perms): EN=0.990, RF=0.240, GB=0.100. Standard
(n_ge+1)/(n_perms+1) correction (used elsewhere in this project, including
`scripts/revision/lib.py`): EN=0.990, RF=0.248, GB=0.109 — recommend the
corrected values for consistency.

**6. Runtime hardware** — Apple M3 Pro (arm64), 18 GB RAM, macOS 26.2
(Darwin 25.2.0) — the same persistent session environment that ran the
`batch_detector.py` timings (217.2s / 781.1s) reported in Section 11.

**7. Fig 7 v3** — `results/revision/figures/fig7_simulation_grid_v3.pdf`.
Axes relabeled "signal/cohort target R²" (no R2-to-f2 conversion exists in
this codebase, per Section 13). `percentile_norm` added as a 6th panel with
an "(assumes same signal direction in every cohort)" footnote. Found and
fixed on inspection: percentile_norm's extreme cells (up to delta_auc≈+0.95)
were saturating the shared color scale and washing out the other 5 panels to
near-white; fixed by keeping vmin/vmax computed from the original 5 methods
only (±0.0911, unchanged from v2) with `extend="both"` colorbar caps so
percentile_norm's cells still render (clipped) rather than distorting
everything else. At the real operating point specifically, percentile_norm's
ΔAUC (-0.0746) is in-range with the other five methods, not an outlier there.

---

## 15. Three small reruns

Full writeup: `results/revision/review2/REVIEW2_TASKS.md`.

**1. PERMDISP for C4 site** — `results/revision/review2/permdisp_c4_site.tsv`.
n=165, 5 sites: F=2.2805, p=0.120 — dispersion homogeneous (unlike the
batch/cohort PERMDISP results elsewhere in this project, which are all
significantly heterogeneous). Site PERMANOVA R² (Task 4: 9.40%, p=0.001)
reflects real centroid separation, not a dispersion artifact.

**2. Milder rarefaction (5,000,000 reads)** —
`results/revision/review2/rarefied_5M_permanova.tsv`,
`dropped_below_5M.tsv`. 29 samples dropped, all from C4 (matches the "about
30, mostly C4" expectation). Final n=254: response R²=0.004382, E0=0.003953,
ΔR²=+0.000429, p=0.186; cohort R²=0.095608, E0=0.011858, ΔR²=+0.083751,
p=0.001. Same conclusion as the deeper (2.1M) rarefaction — robust across
two very different depths.

**3. One-sample discrepancies** — `results/revision/review2/permanova_n281_dropped2.tsv`.
Identified via direct set comparison against `metadata/response_sheet.xlsx`:
C3's extra sample is `SRR6000943` (label R, sourced from this project's own
`metadata/response_labels_PRJNA399742_extended.tsv`, not from the Orletskaia
compilation); C4's extra sample is `ERR10290768` (label NR, sourced from this
project's own `metadata/lee2022_labels.tsv`). Both trace to this project's
own independently-sourced label files, not to the Orletskaia curation, and
nothing suggests either is mislabeled — they simply fall outside that
compilation's own (undocumented, from here) inclusion criteria. Dropping
both and rerunning at n=281: response R²=0.003868, E0=0.003571,
ΔR²=+0.000296, p=0.305; cohort R²=0.105186, E0=0.010714, ΔR²=+0.094471,
p=0.001 — essentially unchanged from the full n=283 result.

**4. Data availability, six unprocessed studies + Gopalakrishnan 2018** —
script `scripts/revision/review2_task4_ena_lookup.py`; outputs
`results/revision/review2/ena_per_run_response_sheet_studies.tsv` and
`ena_project_totals.tsv`; every URL in `REVIEW2_TASKS.md` §4. All six have
**public** raw reads (none controlled-access or undeposited). Whole-BioProject
shotgun size (ENA FASTQ.gz): Gunjur PRJEB49516 107 runs / 348.8 GB (NovaSeq
6000, shotgun; combination nivolumab+ipilimumab, trial CA209-538, no FMT);
Heshiki PRJNA494824 71 runs / 304.0 GB (HiSeq 1500; chemo ± immunotherapy, not
an anti-PD-1 cohort); Liu PRJNA866654 14 runs / 69.8 GB (NSCLC; paper says
"pending review", now public since 2022-10-01); McCulloch PRJNA762360 94
shotgun runs / 397.5 GB + 80 16S runs (anti-PD-1 alone or +peg-IFN in 14 pts);
Spencer PRJNA770295 309 shotgun runs / 557.4 GB + 498 16S runs (observational,
no human FMT); Tsakmaklis PRJNA1011235 29 runs / 92.8 GB (anti-PD-1 15,
anti-PD-1+anti-CTLA-4 13, anti-CTLA-4 1; no FMT). **Gopalakrishnan 2018:** the
paper's statement puts fecal WGS at PRJEB22893 on ENA (25 runs, 56.5 GB, public
since 2017-11-06) and only human exome (EGAS00001002698) on EGA, so the shotgun
data is not controlled-access. **Flag:** the compilation's 22
"Gopalakrishnan_2018" runs resolve to Spencer's PRJNA770295 (46.4 GB), not
PRJEB22893; patient overlap with the 25 PRJEB22893 runs and with the 134
Spencer runs could not be checked (no crosswalk). Not confirmed: Liu's
agents/chemo, whether Spencer patients had combination ICB.
