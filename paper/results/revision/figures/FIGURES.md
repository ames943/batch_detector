# arXiv-version figures

All figures generated from result files already on disk in `results/revision/`
and `results/ml/` — no analysis was rerun. Scripts: `scripts/revision/figures/`.
Shared style: `scripts/revision/figures/style.py` (colorblind-safe palette,
`pdf.fonttype=42`, all text >=7pt, PDF+300dpi PNG for every figure). Every PNG
was opened and inspected after saving; label collisions/clipping found on the
first pass are listed under "fixed" below and are not present in the final files.

---

## fig1_dose_response.{pdf,png} — single/double column (7.1 in)

**Source:** `results/revision/task1_chance_calibrated_permanova.tsv`

**Fixes applied** (vs. the original `results/revision/dose_response.pdf`):
legend moved from inside Panel A to a shared row below both panels (was
overlapping the "10.6" data label); the response and cohort chance curves are
now drawn with distinct line styles (dotted / dash-dot) and each has its own
legend entry ("Chance (response), 1/(n-1)" / "Chance (cohort), (k-1)/(n-1)"),
previously they shared one ambiguous "Chance E0" entry; all Panel B labels
moved to sit above their points with widened y-limits (was: "-0.18" touching
the x-axis, "+0.03" clipped at the bottom).

**Numbers on the figure:**
| n | response R2 | E0 | delta_R2 | cohort R2 | cohort delta_R2 |
|---|---|---|---|---|---|
| 39 | 2.67% | 2.63% | +0.04% | n/a (1 cohort) | n/a |
| 79 | 1.10% | 1.28% | -0.18% | 8.5% | +7.2% |
| 118 | 0.68% | 0.85% | -0.18% | 7.7% | +6.0% |
| 283 | 0.39% | 0.35% | +0.03% | 10.6% | +9.5% |

**Shows:** response R2 tracks its own chance line almost exactly at every n
(and dips slightly below it at n=79/118), while cohort R2 stays 6-9.5
percentage points above its own chance line throughout.

**Unfixed problems:** none.

---

## fig2_ordination.{pdf,png} — double column (7.1 in), 3 panels

**Source:** `results/ml/n283_4cohort/X_genus_clr.tsv` +
`response_labels_n283.tsv` (panels A/B); `results/ml/lee2022/X_genus_clr.tsv` +
`metadata/lee2022_labels.tsv` + `metadata/lee2022_sites.tsv` (panel C).
PCoA (Gower double-centering + eigendecomposition) is a deterministic
projection of these existing CLR matrices for visualization only — no new
statistic was computed. All annotated R2/E0/delta_R2/p values are copied
verbatim from `task1_chance_calibrated_permanova.tsv` (panels A/B) and
`task4_c4_alone_permanova.tsv` (panel C).

**Numbers on the figure:**
- Panel A/B axes: PCo1 = 12.8% variance, PCo2 = 6.6% variance (n=283)
- Panel A text box: Cohort R2=10.56%, E0=1.06%, delta_R2=+9.50%, p=0.001
- Panel B text box: Response R2=0.388%, E0=0.355%, delta_R2=+0.033%, p=0.284
- Panel C axes: PCo1 = 9.2% variance, PCo2 = 7.5% variance (n=165)
- Panel C text box: Site R2=9.40%, delta_R2=+6.96%, p=0.001; Response
  R2=0.66%, delta_R2=+0.05%, p=0.277

**Shows:** the same n=283 ordination visibly separates by cohort (Panel A,
especially cohort4/PRIMM-NL forming its own cluster) but shows no visible
separation by response at all (Panel B) — the identical two axes, just
recolored. Panel C shows C4 alone separates somewhat by recruiting site
(Barcelona stands apart) with no visible response structure.

**Fixes applied:** legends given opaque white backgrounds (`frameon=True,
facecolor="white"`) — they were rendering directly on top of data points in
the crowded upper-right region of each panel with no visual separation.

**Unfixed problems:** none.

### Update — fig2_ordination_v2.{pdf,png} (supersedes the above for the paper)

Per follow-up review, three more issues fixed:
1. Cohort legend now reads **C1/C2/C3/C4** (was `cohort1`-`cohort4`).
2. All legends moved **outside the axes**, below each panel — the opaque
   white backgrounds added above still sat on top of the data (notably the
   C4/cohort4 cluster in panel A), just now legibly.
3. All stats text boxes moved **outside the axes** too, below each panel's
   legend — panel C's box previously covered the lower half of that panel's
   data.

Same source files and same numbers as above (unchanged: PCo1/PCo2 % variance,
R2/E0/delta_R2/p for cohort, response, and site) — only the layout changed.
Figure is now taller (7.5 x 3.6 in) to fit the two-row footer under each panel.

**Unfixed problems:** none.

---

## fig3_chance_calibrated_summary.{pdf,png} — single column (3.5 in)

**Source:** `results/revision/task1_chance_calibrated_permanova.tsv`
(microbiome n=79/118/283, tumor n=223, Duvallet n=570) and
`results/revision/task4_c4_alone_permanova.tsv` (C4 site/response).

**Numbers on the figure (delta_R2%, p):**
| dataset | response delta_R2 | response p | cohort/site delta_R2 | cohort/site p |
|---|---|---|---|---|
| Microbiome n=79 | -0.18% | 0.736 | +7.2%* | 0.001 |
| Microbiome n=118 | -0.18% | 0.874 | +6.0%* | 0.001 |
| Microbiome n=283 | +0.03% | 0.284 | +9.5%* | 0.001 |
| C4 alone n=165 (site) | +0.05% | 0.277 | +7.0%* | 0.001 |
| Tumor n=223 | +0.75%* | 0.003 | +13.4%* | 0.001 |
| CRC (Duvallet) n=570 | +0.44%* | 0.001 | +23.9%* | 0.001 |
(* = p<0.05)

**Shows:** across every dataset in the paper, cohort/site delta_R2 is large
and significant; response delta_R2 is significant only for tumor and Duvallet
(both real but small), and indistinguishable from chance for every microbiome
dataset generated by this project.

**Fixes applied:** response-bar p-value labels were anchored at the bar's own
(near-zero, sometimes negative) x position, landing under the y-axis category
labels; re-anchored to always sit just right of x=0 regardless of sign.

**Unfixed problems:** none.

### Update — fig3_chance_calibrated_summary_v2.{pdf,png} (supersedes the above)

Legend moved from `loc="lower right"` (overlapping the long CRC/Duvallet bar,
+23.9%) to `loc="upper right"` (empty space above the short microbiome-n=79
bar). Same data, same numbers as above — layout only.

**Unfixed problems:** none.

---

## fig4_power_curves.{pdf,png} — single column (3.5 in)

**Source:** `results/revision/task8_power_curve.tsv`,
`results/revision/task8_c4_upper_bound.tsv`.

**Numbers on the figure:** full 6(n) x 4(delta_R2) power grid from
`task8_power_curve.tsv` (e.g. n=165: power=0.12/0.54/0.99/1.00 for
delta_R2=0.0025/0.005/0.01/0.02); C4 95% upper-bound marker at
delta_R2=0.0078, n=131, power=0.80 (from `task8_c4_upper_bound.tsv`:
`dirichlet_boot_95pct_upper_bound_delta_R2`=0.007806,
`n_for_80pct_power_at_bound`=130.6); vertical reference line at n=165 (C4's
actual size).

**Shows:** C4 (n=165) would already have >=80% power for a true delta_R2 of
about 0.0078 or larger — its own 95% bootstrap upper bound on the effect it
actually has. The fact that C4 still lands at p=0.277 (Task 4) is more
consistent with the true effect being near its point estimate (essentially
zero) than with the study being underpowered for an effect of that size.

**Note:** no calibrated power curve exists on disk for
delta_R2=0.0078 exactly (only 0.0025/0.005/0.01/0.02 were run in Task 8) —
only the single (n=131, power=0.80) marker is drawn, not an interpolated
curve. Printed explicitly in the script's log output.

**Unfixed problems (at the time): none** — but see below, the whole figure
turned out to rest on a bug.

### Update — fig4_power_curves_v2.{pdf,png} (supersedes the above; the v1 star was wrong)

**The bug** (found in `scripts/revision/task8_power_hypothetical_effects.py`,
lines 213-224, function `main`): the "n for 80% power at C4's 95% UB" lookup
picked the *nearest calibrated target delta* to 0.007806 out of
`[0.0025, 0.005, 0.01, 0.02]` — that nearest value is **0.01**
(`abs(0.01-0.0078)=0.0022` vs `abs(0.005-0.0078)=0.0028`), NOT 0.0078 itself,
which was never one of the calibrated targets. It then linearly interpolated
n on the **delta_R2=0.01** curve between n=80 (power=0.52) and n=165
(power=0.99), got n=130.6, and reported it as "at delta_R2=0.0078." The v1
star sits on the 0.01 curve, exactly as suspected. This was not an
interpolation-between-curves, and not a different generator — it was a
silent nearest-neighbor substitution disclosed only in a JSON field
(`closest_calibrated_delta_used_for_lookup: 0.01`) that the figure/headline
text didn't surface.

**Second bug**: that power curve (`task8_power_curve.tsv`, v1's Part A) was
built entirely from **Cohort 1's top-30 genera**
(`scripts/revision/task8_power_hypothetical_effects.py` line 149-155,
`base_c1 = raw118.loc[mask, top]...`), while the 95% upper bound itself
(Part B) was correctly computed in **C4's own top-30 genera**. The power
lookup therefore mixed two different feature spaces.

**The fix**: `scripts/revision/task8b_c4_power_v2.py` recomputes the whole
power curve in C4's own feature space — specifically its **full, unfiltered
genus set** (p=3005), matching `scripts/revision/task4_c4_alone.py` exactly
(confirmed by reading that script: it calls `clr_transform(raw.values)`
directly with no prevalence/variance filter at all — "the same filtering as
Task 4" turned out to mean *no filtering*). Same injected-shift binary-search
calibration method as Task 8. Targets **now include 0.0078 directly**
(no more nearest-neighbor substitution): delta_R2 in
{0.0025, 0.005, 0.0078, 0.01}, n in {80, 131, 165, 250, 400}, 200 reps/cell.
Output: `results/revision/task8_c4_power_v2.tsv`.

**Result — the corrected numbers are dramatically different, in the
opposite direction from what v1 implied:**
| n | ΔR²=0.0025 | ΔR²=0.005 | ΔR²=0.0078 | ΔR²=0.01 |
|---|---|---|---|---|
| 80 | 0.145 | 0.585 | 0.855 | 0.885 |
| 131 | 0.385 | 0.875 | 0.995 | 1.000 |
| 165 | 0.490 | 0.935 | **0.995** | 1.000 |
| 250 | 0.835 | 1.000 | 1.000 | 1.000 |
| 400 | 0.980 | 1.000 | 1.000 | 1.000 |

**Answers to the specific questions asked:**
- Power at n=165 for ΔR²=0.0078: **0.995** (v1's mislabeled figure said 0.80).
- n reaching 80% power at ΔR²=0.0078: already reached by **n=80**, the
  smallest n tested (power=0.855 there) — nowhere near needing n≈131 or 165.
- Minimum detectable ΔR² (80% power) at n=165: **≈0.0042** (interpolated
  between the computed n=165 points at ΔR²=0.0025→power=0.49 and
  ΔR²=0.005→power=0.935) — about half of the 0.0078 upper bound.

**Does the conclusion still hold?** Yes, and more strongly than v1 claimed.
C4 (n=165) had ~99.5% power to detect an effect at its own 95% upper bound —
not a marginal ~80%. Its minimum detectable effect at 80% power (~0.0042) is
roughly half the upper bound itself. C4 was **not** underpowered for an
effect anywhere near what its own data could plausibly contain; the null
result (Task 4: p=0.277) is therefore better explained by the true effect
being close to its point estimate (~0) than by insufficient sample size.

**Why v1's C1-based curve showed so much lower power at the same (n, delta_R2)**
than C4's own feature space: likely the much higher dimensionality of C4's
full genus set (p=3005 vs. the top-30 subset) — with the injected shift
recalibrated to hit the same mean R2 either way, a higher-dimensional
composition appears to reduce the sampling variance of the Aitchison-distance
PERMANOVA statistic across resamples, increasing power for a matched mean
effect. This is a plausible mechanism, not independently verified by a
dedicated ablation — flagged here rather than asserted as certain.

New figure: adds the ΔR²=0.0078 curve (now actually computed, so drawn per
instructions); the star sits at (n=165, power=0.995), an ACTUAL COMPUTED grid
point (v1's star was an interpolated point on the wrong curve); the
"C4 (n=165)" label is now offset above-right of the dotted line instead of
sitting on top of it.

**Unfixed problems:** none (at the time) — one cosmetic layout issue found
next.

### Update — fig4_power_curves_v3.{pdf,png} (layout-only; supersedes v2's layout, same data)

Two layout fixes, made after review of v2: (1) the "power=0.995 at C4's
own n and 95% UB effect" annotation was up near the saturated curves at the
top of the plot; moved to empty space in the lower right with a thin leader
line back to the star. (2) The "C4 (n=165)" label had its own leader line,
which read as a stray horizontal segment disconnected from anything;
replaced with plain offset text (no line) clear of the dotted vertical line.
Legend also moved to upper-left (empty space) since lower-right is now
occupied by the relocated annotation. Same `task8_c4_power_v2.tsv` data —
no numbers changed.

**Unfixed problems:** none.

### Update — fig4_power_curves_v4.{pdf,png} (layout-only; supersedes v3's layout, same data)

`results/revision/final/fig4_power_curves_v4.{pdf,png}`. Same
`task8_c4_power_v2.tsv` data as v2/v3 — no numbers changed. Six fixes:
1. Legend moved fully outside the axes, below the x-axis label, in 2 columns
   (was inside the plot, covering the orange ΔR²=0.005 point at n=80 and
   crossing the green ΔR²=0.0078 line).
2. x-axis ticks replaced with plain numbers at the 5 tested sample sizes (80,
   131, 165, 250, 400); log scale kept; minor tick labels removed.
3. The long leader line from "power=0.995 at C4's own n and 95% UB effect"
   (which crossed the orange and blue curves) removed entirely; replaced with
   a short "0.995" label immediately to the right of the star, no line
   needed (the paper caption explains the star).
4. "C4 (n=165)" label kept above the plot, clear of the dotted line
   (unchanged from v3).
5. 80%-power dashed line and label kept; label repositioned to the far left
   (n≈85, just above the line) where the blue ΔR²=0.0025 curve is still low
   (~0.15-0.4), so it no longer risks touching that curve as n grows.
6. Figure made taller (3.5 x 3.9 in, was 3.5 x 3.2 in) to fit the two-row
   legend below while keeping single-column width.

**Unfixed problems:** none.

---

## fig5_cv_vs_null.{pdf,png} — double column (7.1 in), 4x3 grid

**Source:** `results/revision/task5_{c1,n118,n283,c4}_fold_aucs.tsv`,
`..._permutation_aucs.tsv`, `..._permutation_summary.tsv`.

**Numbers on the figure (observed mean AUC, permutation p):**
| dataset | EN | RF | XGB |
|---|---|---|---|
| C1 (n=39) | 0.390, p=0.752 | 0.431, p=0.653 | 0.392, p=0.762 |
| n=118 | 0.376, p=0.950 | 0.546, p=0.248 | 0.541, p=0.248 |
| n=283 | 0.508, p=0.455 | 0.537, p=0.228 | 0.533, p=0.198 |
| C4 (n=165) | 0.571, p=0.109 | 0.502, p=0.485 | 0.531, p=0.267 |

**Shows:** for every dataset x model cell, the 50 individual per-fold AUCs
(blue strip) and their observed mean (black diamond) sit inside the spread of
the 100-permutation null distribution of mean AUC (gray violin) — visually
confirming none of the 12 results are significant at alpha=0.05.

**Fixes applied:** the shared legend and per-column model titles
("ElasticNet"/"RandomForest"/"XGBoost") were overlapping at the top of the
figure; re-laid out with explicit vertical spacing (`fig.text` for the title,
`fig.legend` pinned below it, more headroom via `subplots_adjust(top=0.88)`).

**Unfixed problems:** none.

### Update — fig5_cv_vs_null_v2.{pdf,png} (supersedes the above)

Three fixes: (1) removed the suptitle entirely (the paper caption states the
protocol, so it was redundant, and it's also where a factual error was — see
next point); (2) the removed suptitle had said "repeated 5x10-fold CV" — the
actual protocol (`scripts/revision/task5_repeated_cv.py`) is **10 repeats of
stratified 5-fold CV** (50 outer folds total), not "5x10-fold"; (3) column
titles renamed **"Elastic net" / "Random forest" / "Gradient boosting"**
(were "ElasticNet"/"RandomForest"/"XGBoost"). Same data, same numbers as
above.

**Unfixed problems:** none.

---

## fig6_volcano.{pdf,png} — single column (3.5 in)

**Source:** `results/revision/task6_meta_analysis_hedges_g.tsv`.

**Numbers on the figure:** "125 genera, 0 FDR-significant, mean I2 = 14.4%";
q=0.05 threshold line at -log10(q)=1.301; the 5 smallest-q genera labeled —
Dialister (g=-0.262, q=0.957), Pseudomonas (g=+0.216, q=0.957),
Ligilactobacillus (g=+0.252, q=0.957), Phascolarctobacterium (g=+0.262,
q=0.957), Wujia (g=+0.329, q=0.870). Maximum observed -log10(q) across all
125 genera is 0.060 — no point comes anywhere near the significance line.

**Shows:** even the "best" (smallest-q) genus from the corrected meta-analysis
is nowhere close to FDR significance; the null is a wide, clean margin, not a
borderline call.

**Fixes applied:** the initial y-axis spanned all the way to the data's actual
near-zero range PLUS the empty gap up to q=0.05, wasting most of the panel;
capped the y-axis just above the threshold line. The top-5 labels initially
overlapped each other (4 of the 5 genera sit within a narrow x-range); replaced
with a staggered vertical label stack + thin leader lines to each point, and
widened the x-axis so the longest label ("Phascolarctobacterium") no longer
clipped at the right edge.

**Unfixed problems:** none (at the time) — on closer inspection the 4
right-side leader lines still read as crossing; fixed next.

### Update — fig6_volcano_v2.{pdf,png} (layout-only; same data)

`results/revision/final/fig6_volcano_v2.{pdf,png}`. Same
`task6_meta_analysis_hedges_g.tsv` data — no numbers changed. Three of the
four right-side points (Pseudomonas, Ligilactobacillus, Phascolarctobacterium)
sit within 0.046 units of each other in x, three of the four share the exact
same q-value/y-height — at that proximity, any pair of individually-offset
diagonal leader lines with overlapping x-ranges reads as crossing to the eye
even when the segments don't technically intersect. Fix: anchor all 4 labels
at a single fixed x to the right of the cluster, stacked top-to-bottom in the
same left-to-right order as their points (rightmost point <-> topmost label,
Dialister kept separate on the left as before) — spokes from one common
vertical line to monotonically-ordered points cannot cross by construction
(same principle as non-crossing dendrogram tip labels). Kept in the lower
half (labels span y=0.10-0.38, points sit at y=0.02-0.06, well below the
q=0.05 line at y=1.30).

Base marker size increased 10->20 (with the 5 highlighted points slightly
larger still, at 22, open circles) so the 125-point cluster near y=0 reads as
a scatter rather than a flat smear — axis scale and limits otherwise
unchanged, still no broken axis; the point (nothing comes near q=0.05)
is unchanged and still visually obvious.

**Unfixed problems:** none.

---

## fig7_simulation_grid.{pdf,png} — double column (7.1 in), 5 panels

**Source:** `results/ml/simulation/grid_results.tsv`, filtered to
n_cohorts=3, n_per_cohort=40.

**Numbers on the figure:** delta_AUC heatmap (diverging, centered at 0,
vmin=-0.091/vmax=+0.091) for methods none / mean-centering / location-scale /
quantile mapping / cohort covariate, across the 4x4 grid of signal
f^2 in {0.007, 0.027, 0.05, 0.10} x cohort f^2 in {0.00, 0.04, 0.08, 0.15}.
Real-operating-point cell (signal f^2=0.027, cohort f^2=0.08, outlined in
black in every panel) delta_AUC by method: none=-0.069, mean_centering=-0.030,
location_scale=-0.040, quantile_mapping=-0.040, cohort_covariate=-0.003.

**Shows:** at the real operating point, every correction method is at best
roughly neutral (cohort_covariate) to clearly negative (none) relative to the
single-cohort ceiling — none of the 5 methods recovers the ceiling AUC there,
consistent with the paper's "correction-helps zone is narrow and the real
data doesn't sit in it" finding.

**percentile_norm excluded as planned**: its simulated gains are a known
generator artifact (the simulation assumes a consistent cross-cohort signal
direction, which the real-data Phase 0.5 cross-cohort holdout showed does not
hold) — logged explicitly by the script, not silently dropped.

**Unfixed problems:** none.

### Update — fig7_simulation_grid_v2.{pdf,png} (supersedes the above)

The outlined "real operating point" cell was wrong: it was at (signal
f2=0.027, cohort f2=0.08), but per `results/ml/simulation/calibration_deltas.tsv`
and `scripts/phase2_simulation.py` (`REAL_SIG_F2_N118=0.007`, comment
`"sig_f2=0.007 row = NULL MODEL (delta=0 from calibration)"`), **signal
f2=0.007 is the actual real operating point** — it was calibrated to an
injected shift of delta=4.657e-09 (numerically zero, i.e. no injected signal
at all), chosen specifically because its simulated R2 (0.00851) matches the
REAL measured n=118 pooled response R2 (~0.007-0.0085, itself indistinguishable
from chance per Task 1). Signal f2=0.027 is a *different* calibrated point
(delta=2.987, verified R2=0.0265) matching Cohort-1-ALONE's (n=39) response
R2, not the pooled n=118 point this figure is meant to mark. Outline moved to
(0.007, 0.08); suptitle removed (paper caption covers it).

Corrected real-operating-point cell (signal f2=0.007, cohort f2=0.08)
delta_AUC by method: none=-0.072, mean_centering=-0.059, location_scale=-0.036,
quantile_mapping=-0.041, cohort_covariate=-0.066 — all methods more negative
than at the old (wrong) cell, i.e. correction looks even less helpful at the
actual real operating point than the mislabeled cell suggested.

### Update — fig7_simulation_grid_v3.{pdf,png} (adds a panel; same data)

Two changes: (1) axes relabeled "signal target $R^2$" / "cohort
target $R^2$" (were "signal $f^2$" / "cohort $f^2$") — no R2-to-Cohen's-f2
conversion is ever applied anywhere in this codebase (confirmed by reading
`scripts/phase2_simulation.py`; see `SUMMARY.md` Section 13), so this is a
pure relabeling, not a data change. (2) `percentile_norm` added back as a
sixth panel, with a figure-level footnote noting it "assumes the same signal
direction in every cohort" (kept as a footnote rather than a per-panel
subtitle, to keep panel titles short at 6-across).

**Fix found on inspection**: naively including percentile_norm's cells (up
to delta_AUC≈+0.95 elsewhere in its grid — a known generator artifact) in
the shared color scale saturated it and washed out all five other panels to
near-uniform white. Fixed by computing vmin/vmax from the original 5 methods
only (±0.0911, unchanged from v2); percentile_norm's own cells still render
(clipped to the colorbar's end color), with `extend="both"` triangle caps
added to the colorbar so the clipping is visually flagged rather than silent.
At the real operating point specifically, percentile_norm's ΔAUC (-0.0746)
is in-range with the other five methods — its extreme cells live elsewhere
in the grid (higher signal levels, where its signal-direction-consistency
assumption is exercised more).

**Unfixed problems:** none.

**Unfixed problems:** none.
