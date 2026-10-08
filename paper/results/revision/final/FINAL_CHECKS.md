# Final checks before the paper is locked

Rules followed: no existing result file overwritten; new outputs live in
`results/revision/final/`; every new number is also appended to
`results/revision/SUMMARY.md` Section 13 with its source.

---

## 1. Fig 4 power figure — layout fix

`results/revision/final/fig4_power_curves_v3.{pdf,png}` — same data as
`fig4_power_curves_v2` (`results/revision/task8_c4_power_v2.tsv`), two layout
changes:
- "power=0.995 at C4's own n and 95% UB effect" moved to empty space in the
  lower right, with a thin leader line to the star at (n=165, power=0.995).
- The "C4 (n=165)" label's leader line (the stray-looking horizontal segment)
  removed; it's now plain offset text clear of the dotted vertical line.

PNG opened and checked after saving: no collisions, no clipping, no
overlapping text. See `results/revision/figures/FIGURES.md` for the full
figure-log entry.

---

## 2. batch_detector recommendation thresholds — rationale check

**Verdict up front: the 0.25x / 4x cutoffs are NOT supported by the grid as
described. The docstring's derivation is wrong — it appears to conflate two
different comparisons.** No code was changed — this is a
paper-wording fix only.

### (a) Exact current docstring text (`scripts/batch_detector.py`, `_pooling_recommendation_v2`, lines ~759-764)

> "The rho thresholds (0.25x / 4x) come directly from the companion Phase-2
> calibrated simulation, which found batch correction only reliably recovers
> a no-batch ceiling AUC when the biological signal is at least ~4x the batch
> effect size (sig_f2 >= ~4x batch_f2 in that grid) — i.e. rho <= 0.25. The
> mirror-image threshold (rho >= 4) marks where batch dominates by the same
> margin in the other direction."

### (b) Every cell with ΔAUC > 0 (nc=3, npc=40), with its signal:cohort f² ratio

Source: `results/ml/simulation/grid_results.tsv`, filtered to
`n_cohorts==3 & n_per_cohort==40`, `percentile_norm` excluded (already flagged
elsewhere as a generator artifact). All 29 positive-ΔAUC cells:

| method | signal f² | cohort f² | signal:cohort ratio | ΔAUC |
|---|---|---|---|---|
| cohort_covariate | 0.027 | 0.04 | 0.68 | +0.0014 |
| cohort_covariate | 0.100 | 0.00 | inf | +0.0314 |
| cohort_covariate | 0.100 | 0.04 | 2.50 | +0.0615 |
| cohort_covariate | 0.100 | 0.08 | **1.25** | +0.0617 |
| cohort_covariate | 0.100 | 0.15 | **0.67** | +0.0056 |
| location_scale | 0.027 | 0.00 | inf | +0.0035 |
| location_scale | 0.100 | 0.00 | inf | +0.0418 |
| location_scale | 0.100 | 0.04 | 2.50 | +0.0497 |
| location_scale | 0.100 | 0.08 | **1.25** | +0.0354 |
| location_scale | 0.100 | 0.15 | **0.67** | +0.0300 |
| mean_centering | 0.027 | 0.04 | 0.68 | +0.0090 |
| mean_centering | 0.100 | 0.00 | inf | +0.0301 |
| mean_centering | 0.100 | 0.04 | 2.50 | +0.0444 |
| mean_centering | 0.100 | 0.08 | **1.25** | +0.0444 |
| mean_centering | 0.100 | 0.15 | **0.67** | +0.0513 |
| none | 0.100 | 0.00 | inf | +0.0439 |
| none | 0.100 | 0.04 | 2.50 | +0.0542 |
| none | 0.100 | 0.08 | **1.25** | +0.0232 |
| none | 0.100 | 0.15 | **0.67** | +0.0025 |
| quantile_mapping | 0.027 | 0.00 | inf | +0.0078 |
| quantile_mapping | 0.027 | 0.15 | 0.18 | +0.0192 |
| quantile_mapping | 0.050 | 0.00 | inf | +0.0822 |
| quantile_mapping | 0.050 | 0.04 | 1.25 | +0.0193 |
| quantile_mapping | 0.050 | 0.08 | 0.625 | +0.0146 |
| quantile_mapping | 0.050 | 0.15 | 0.33 | +0.0713 |
| quantile_mapping | 0.100 | 0.00 | inf | +0.0356 |
| quantile_mapping | 0.100 | 0.04 | 2.50 | +0.0289 |
| quantile_mapping | 0.100 | 0.08 | **1.25** | +0.0265 |
| quantile_mapping | 0.100 | 0.15 | **0.67** | +0.0628 |

**This directly contradicts the docstring's claim.** At signal f²=0.10 (the
grid's highest signal level), four of five methods still post positive ΔAUC
at cohort f²=0.15 — a ratio of **0.67, i.e. signal SMALLER than the cohort
effect**, the opposite of "signal >= 4x cohort." There is no ratio threshold
in the data at all: positive ΔAUC cells span ratios from 0.18 to infinity,
with no clustering near 4 or 0.25. The one clean pattern in the grid is that
positive cells occur almost exclusively at **signal f²=0.10** (the top
absolute signal level), regardless of the cohort effect size or their ratio.

Confirmed: **0.10 / 0.027 ≈ 3.7 ≈ "~4x"** — this is a
signal-to-signal comparison across the grid's own two highest calibrated
levels (0.10 vs. 0.027, the real per-cohort operating point from the old
Cohort-1-alone framing), not a signal-to-cohort ratio. The docstring's "~4x"
almost certainly is that number, mislabeled as a rho threshold.

### (c) f² vs. R² — is there a conversion?

**No.** Checked `scripts/phase2_simulation.py` end to end: `sig_f2` and
`batch_f2` are calibrated directly against PERMANOVA R² as the target
(`calibrate_delta(target_r2, ...)` — the parameter is literally named
`target_r2`), and the code comments confirm this identification directly,
e.g. line 76: `REAL_SIG_F2_N118 = 0.007   # PERMANOVA response R² pooled
n=118`. There is no `R2/(1-R2)` (or any other) Cohen's-f²-from-R²
transformation anywhere in the script. In this codebase, "f²" is a label
applied directly to PERMANOVA R² values, not Cohen's f² in the textbook
sense. Since every R² value in the grid is small (max 0.10), the informal
usage and the textbook definition are numerically close (at R²=0.10, true
Cohen's f²=0.111, an ~11% difference; at R²=0.027, f²=0.0277, ~2.6% off) — so
this labeling looseness does not change the substance of the finding in (b),
but `rho` (`batch_detector.py`'s literal `R2_batch/R2_response`) and the
grid's `sig_f2`/`batch_f2` columns can be compared directly as ratios without
any further conversion, because none was ever applied on the grid side either.

### (d) Are the 0.25x / 4x cutoffs supported by the grid?

**No.** The paper's sentence should be replaced. Two options:

- Plain version, no false precision: *"The tool's pooling recommendation
  applies conservative, heuristic ρ thresholds (0.25x and 4x) chosen to be
  cautious relative to the calibrated simulation, which found no reliable
  ρ-based cutoff for when batch correction restores the single-cohort AUC
  ceiling; instead, correction only reliably helped at the grid's highest
  absolute signal level (f²=0.10), independent of the cohort effect size."*
- Slightly shorter variant: *"These are
  heuristic cutoffs, conservative relative to the simulation grid, which did
  not show a clean ρ-based threshold for when correction helps — positive
  ΔAUC cells spanned signal:cohort ratios from well below 1 to infinity, with
  the only consistent pattern being the grid's highest absolute signal level
  (f²=0.10)."*

Either is accurate. The original sentence ("reliably recovered the
no-cohort-effect ceiling only when the biological effect was about four times
larger than the cohort effect") should not be used — it is not what the grid
shows.

---

## 3. Optional — feature-space dimensionality vs. power (≤1 hour)

Rerun: `scripts/revision/task8c_c4_power_top30.py`, same method as
`task8b_c4_power_v2.py`, restricted to C4's own top-30 most prevalent genera
(same 165 C4 samples, only the feature dimensionality changes: p=30 vs.
p=3005). n=165, ΔR² in {0.005, 0.0078}, 200 reps. Output:
`results/revision/final/task8c_c4_power_top30.tsv`.

| source | feature space | n | ΔR²=0.005 power | ΔR²=0.0078 power |
|---|---|---|---|---|
| `task8_power_curve.tsv` (original) | Cohort 1, top-30 | 165 | 0.54 | ~0.79 (interpolated; 0.0078 wasn't directly computed there) |
| `task8c_c4_power_top30.tsv` (new) | **C4, top-30** | 165 | **0.625** | **0.900** |
| `task8_c4_power_v2.tsv` | **C4, full genus set (p=3005)** | 165 | **0.935** | **0.995** |

**Result: dimensionality is the dominant driver, but not the entire story.**
Holding the cohort fixed (C4) and only shrinking the feature space from
p=3005 to p=30 drops power substantially (0.935→0.625 at ΔR²=0.005;
0.995→0.900 at ΔR²=0.0078) — most of the original C1-vs-C4 gap is reproduced
by dimensionality alone. But C4's own top-30 power (0.625 / 0.900) is still
noticeably higher than Cohort 1's top-30 power (0.54 / ~0.79) at the same n
and target ΔR² — so a real, smaller cohort-specific component (differing
covariance/dispersion structure between C1's and C4's top-30 genera) also
contributes. Both are real effects; dimensionality is the larger one.

---

## 4. Number verification against result files

| # | Claim | File(s) | Value found | Verdict |
|---|---|---|---|---|
| 1 | C1 Dirichlet bootstrap 95% CI, response R²: [0.0181, 0.0295] | `results/ml/batch_diagnostics/bootstrap_cis.tsv` (row `PERMANOVA_response_R2`, n=39) | ci_lower=0.018109, ci_upper=0.029493 | **Matches** ([0.0181, 0.0295]) |
| 2 | PERMDISP response p, n=39/79/118/283: 0.341/0.192/0.779/0.366 | `results/ml/batch_diagnostics/permdisp_results.tsv` | n39: 0.3413; n79: 0.1922; n118: 0.7788; n283: 0.366 | **Matches** |
| 3 | PERMDISP cohort: n=79 F=9.59 p=0.005; n=118 F=5.63 p=0.006; n=283 p=0.001 | same file, `factor=batch` rows | n79: F=9.5885, p=0.005; n118: F=5.6285, p=0.006; n283: p=0.001 (F=9.7644) | **Matches** |
| 4 | PERMDISP phylum/genus/species (n=118): response p 0.969/0.762/0.267; cohort p 0.008/0.002/0.015 | `results/ml/phase3c/taxonomic_sensitivity_summary.tsv` | response: 0.969/0.7618/0.2673; batch: 0.008/0.002/0.015 | **Matches** |
| 5 | Nested LOOCV gradient boosting: depth 2, 100 trees, lr 0.05 in 82-85% of folds (n=118) | `results/ml/nested_cv/nested_cv_best_params.tsv` | Exact triple (depth=2, n_est=100, lr=0.05): **34/118 = 28.8%** of folds. Marginals: lr=0.05 alone in **82/118 = 69.5%** of folds; n_est=100 alone in **85/118 = 72.0%** of folds; depth=2 alone in 83/118=70.3%. | **NOT supported as written.** The numbers 82 and 85 are real — but they are **raw fold counts** for two different single-hyperparameter marginals (lr=0.05: 82 folds; n_estimators=100: 85 folds), not a joint percentage. The exact combination (all three together) occurs in under 29% of folds. This reads like "82 and 85 folds" got rewritten as "82-85%." Recommend: *"Inner CV most often selected learning_rate=0.05 (82/118 folds, 69.5%) and n_estimators=100 (85/118 folds, 72.0%); the full combination (depth=2, n_estimators=100, learning_rate=0.05) was selected in 34/118 folds (28.8%)."* |
| 6 | C1 elastic net LOOCV AUC 0.449, p=0.560 | `results/ml/cross_cohort/holdout_summary.tsv`, `REF` row, train=test="1 (LOOCV)" | AUC=0.4486, permutation_p=0.56 | **Matches** |
| 7 | Sign inversion: 48 genera, 64.6% vs 58.8% null, p=0.213 | `results/ml/phase3b/inversion_permutation_test.tsv` | obs_genera_analyzed=48, obs_inversion_rate=0.6458, null_mean=0.5878, p=0.213 | **Matches** |
| 8a | Tumor pooled AUCs 0.45-0.52 | `logs/tumor_ml.log` (not a structured TSV — the pooled LOOCV loop in `scripts/tumor_ml.py` only prints these, never saves them to a file) | ElasticNet=0.5210, RandomForest=0.4698, XGBoost=0.4482 | **Matches range**, but flag: only traceable via a log file, not a saved result table |
| 8b | XGB AUC 0.448 p=0.745 | `results/ml/tumor/tumor_permutation_results.tsv` | observed_auc=0.448167, p_value=0.745 | **Matches** |
| 8c | Liu EN AUC 0.513 p=0.290 | `results/ml/tumor/single_cohort/liu2019_permutation.tsv` | observed_auc=0.512564, p_value=0.29 | **Matches** |
| 8d | Hugo→Liu 0.58-0.62 | `results/ml/tumor/tumor_loso_results.tsv` | See below — **direction is mislabeled** | **Not supported as written** |
| 8e | tested-on-Riaz 0.34-0.45 | same file, `held_out=riaz2017` | EN=0.4359, RF=0.4469, XGB=0.3391 | **Matches** |
| 8f | median TMB 4.7/6.5/13.1 | `results/ml/tumor/combined_tumor_features.tsv`, grouped median TMB by study (computed directly from this file; not a pre-saved statistic) | riaz2017=4.658, liu2019=6.474, hugo2016=13.079 | **Matches** (riaz/liu/hugo order) |
| 8g | Top Liu features include TTN, MUC16, ADGRV1; RAC1 the only canonical driver in top 10 | `logs/liu_single_cohort.log` (XGBoost importance) vs. `results/ml/tumor/liu2019_features.tsv` (raw mutation frequency, computed directly from the file) | See below — **conflates two different rankings** | **Not supported as written** |
| 9 | Percentile normalization reduced response R² by 81% (0.0068→0.0013) | `results/ml/batch_correction_shootout/method_comparison.tsv` | uncorrected response_R2=0.0068; percentile_norm=0.001261 → 81.5% reduction | **Matches** |
| 10 | Batch-correction table (n=118 LOOCV): response R², EN AUC, RF AUC, all 6 methods | same file | See table below | **Matches** |

**#8d detail — LOSO direction is mislabeled.** `scripts/tumor_ml.py` lines
68-70: `train_mask = studies != held_out; test_mask = studies == held_out` —
"held_out" is the **test** cohort, trained on the **other two pooled**. So
`held_out=hugo2016` (EN=0.6103, RF=0.6191, XGB=0.5765 → range ~0.58-0.62,
matching the claimed numbers) means **trained on Liu2019+Riaz2017 pooled,
tested on Hugo2016** — not "Hugo→Liu" (which would mean trained on Hugo
alone, tested on Liu alone; no such single-cohort-to-single-cohort result
exists anywhere in the tumor pipeline — only this 2-cohorts-pooled-vs-1-held-out
design was run). The numbers are real and correctly sourced; the **direction
in the label is backwards**. Recommend: *"(Liu2019+Riaz2017)→Hugo2016 AUC
0.58-0.62"* or simply *"tested-on-Hugo 0.58-0.62"* (matching how
tested-on-Riaz is already phrased correctly).

**#8g detail — two different "top features" lists got merged.** By XGBoost
feature importance (`logs/liu_single_cohort.log`, full-data fit), the actual
top 10 for Liu 2019 are: **RAC1** (0.0484), DNAH5 (0.0419), OBSCN (0.0393),
**ADGRV1** (0.0353), NEB (0.0345), PKHD1L1 (0.0327), CSMD2 (0.0306), FAT4
(0.0303), MAP2K1 (0.0293), KIT (0.0290). **TTN and MUC16 are not in this list
at all.** Separately, by raw mutation frequency (computed directly from
`results/ml/tumor/liu2019_features.tsv`, mean of each `mut_*` binary column),
TTN (63.9% of patients) and MUC16 (59.7%) are indeed the two most frequently
mutated genes in Liu 2019 — consistent with them being famous
large/passenger-mutation genes — but that is a **different ranking** (mutation
prevalence, not predictive importance) and RAC1/ADGRV1 do not appear near the
top of it. The paper's sentence conflates the two rankings. Recommend:
*"By raw mutation frequency, TTN (63.9%) and MUC16 (59.7%) — well-known
large, frequently-mutated passenger genes — dominate Liu 2019; by XGBoost
feature importance, however, the top predictor is RAC1, a canonical MAPK
pathway driver, followed by large/passenger genes (DNAH5, OBSCN, ADGRV1, NEB,
PKHD1L1, CSMD2, FAT4); RAC1 is the only canonical driver among the top 10 by
importance."* This preserves the intended scientific point (the model mostly
rediscovers TMB-correlated passenger mutations rather than real driver
biology) while fixing the specific gene attributions.

**Batch-correction table (n=118, LOOCV)** — source
`results/ml/batch_correction_shootout/method_comparison.tsv`:

| method | response R² | ElasticNet AUC | RandomForest AUC |
|---|---|---|---|
| uncorrected | 0.0068 | 0.3511 | 0.5421 |
| mean_centering | 0.0060 | 0.3294 | 0.5321 |
| location_scale | 0.0060 | 0.3634 | 0.5127 |
| quantile_mapping | 0.0057 | 0.3466 | 0.4193 |
| percentile_norm | 0.0013 | 0.4193 | 0.4647 |
| cohort_covariate | 0.0068 | 0.3537 | 0.3439 |

All six rows match the paper's table exactly.

---

## 5. References

Checked against CrossRef (`api.crossref.org`) for all but one; the
Orletskaia & Olekhnovich preprint's CrossRef lookup was rate-limited (HTTP
429), so that one was confirmed via web search (bioRxiv page + PubMed +
publisher record) instead.

| Reference | Result |
|---|---|
| Peng et al. 2020, Cancer Immunol. Res. 8(10):1251-1261 | **Matches exactly.** Title: "The Gut Microbiome Is Associated with Clinical Response to Anti–PD-1/PD-L1 Immunotherapy in Gastrointestinal Cancer." DOI: 10.1158/2326-6066.CIR-19-1014. |
| Orletskaia & Olekhnovich 2025, bioRxiv doi:10.1101/2025.05.07.652660 | **DOI and title match exactly**: "Ecological and functional stratification of the stool microbiome predicts response to immune checkpoint inhibitors across cancer types." **Update worth noting**: this work has since been peer-reviewed and published in *Computational and Structural Biotechnology Journal* (doi:10.34133/csbj.0065; PMID 42146906) — you may want to cite the peer-reviewed CSBJ version instead of, or alongside, the bioRxiv preprint before final submission. |
| Anderson 2001, Austral Ecol. 26(1):32-46 | **Matches exactly.** "A new method for non-parametric multivariate analysis of variance." DOI: 10.1046/j.1442-9993.2001.01070.x. |
| McArdle & Anderson 2001, Ecology 82(1):290-297 | **Matches exactly.** "Fitting multivariate models to community data: a comment on distance-based redundancy analysis." DOI: 10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2. |
| DerSimonian & Laird 1986, Control. Clin. Trials 7(3):177-188 | **Matches exactly.** "Meta-analysis in clinical trials." DOI: 10.1016/0197-2456(86)90046-2. |
| Parker, Günter & Bedo 2007, BMC Bioinformatics 8:326 | **Matches exactly** (volume 8, article 326). Title: "Stratification bias in low signal microarray studies." DOI: 10.1186/1471-2105-8-326. |
| Airola et al. 2011, Comput. Stat. Data Anal. 55(4):1828-1844 | **Matches exactly.** "An experimental comparison of cross-validation techniques for estimating the area under the ROC curve." DOI: 10.1016/j.csda.2010.11.018. |
| Martín-Fernández et al. 2003, Math. Geol. 35(3):253-278 | **Matches exactly.** Title: "Dealing with Zeros and Missing Values in Compositional Data Sets Using Nonparametric Imputation." DOI: 10.1023/A:1023866030544. |

All eight references check out on volume/issue/pages/DOI. Only actionable
item: the now-published CSBJ version of Orletskaia &
Olekhnovich can be cited alongside or instead of the bioRxiv preprint.
