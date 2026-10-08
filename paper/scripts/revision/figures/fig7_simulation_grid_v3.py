#!/usr/bin/env python3
"""
Fig 7 (v3) — Phase-2 simulation grid: Delta_AUC vs. single-cohort ceiling, nc=3, npc=40.
Source: results/ml/simulation/grid_results.tsv,
        results/ml/simulation/calibration_deltas.tsv (no recompute).

v3 changes vs. v2:
  1. Axes relabeled "signal target R2" / "cohort target R2" (were "signal f2" /
     "cohort f2") -- the paper now calls these target R2 values (see
     SUMMARY.md Section 13: no R2-to-Cohen's-f2 conversion is ever applied in
     this codebase; "f2" was always a direct relabeling of PERMANOVA R2).
  2. percentile_norm added back as a sixth panel, with a note in its title:
     "(assumes same signal direction in every cohort)" -- its own generator
     artifact, now shown rather than omitted.
Outline kept on the real operating point (signal target R2=0.007, cohort
target R2=0.08) -- see v2 docstring / SUMMARY.md for why 0.007, not 0.027, is
the correct cell.
"""
import sys
sys.path.insert(0, "scripts/revision/figures")
import numpy as np
import pandas as pd
from style import apply, save

apply()
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

df = pd.read_csv("results/ml/simulation/grid_results.tsv", sep="\t")
sub = df[(df["n_cohorts"] == 3) & (df["n_per_cohort"] == 40)].copy()

METHODS = ["none", "mean_centering", "location_scale", "quantile_mapping",
           "cohort_covariate", "percentile_norm"]
TITLES = {"none": "none", "mean_centering": "mean-\ncentering", "location_scale": "location-\nscale",
          "quantile_mapping": "quantile\nmapping", "cohort_covariate": "cohort\ncovariate",
          "percentile_norm": "percentile\nnorm*"}

sig_vals = sorted(sub["sig_f2"].unique(), reverse=True)  # high->low, top->bottom
bat_vals = sorted(sub["batch_f2"].unique())               # low->high, left->right

# Color scale range set from the 5 original methods only (not percentile_norm):
# percentile_norm's known generator-artifact cells (delta_auc up to ~+0.95,
# see calibration notes) would otherwise saturate the shared scale and wash
# out every other panel to near-uniform white. percentile_norm's own cells
# still render -- they just clip to the colorbar's end color, which itself
# communicates "far outside the other methods' range."
CORE_METHODS = ["none", "mean_centering", "location_scale", "quantile_mapping", "cohort_covariate"]
vmax = sub[sub["method"].isin(CORE_METHODS)]["delta_auc"].abs().max()
vmin = -vmax

REAL_SIG, REAL_BAT = 0.007, 0.08

fig, axes = plt.subplots(1, len(METHODS), figsize=(7.5, 2.5), sharey=True)
fig.subplots_adjust(left=0.085, right=0.89, top=0.80, bottom=0.34, wspace=0.14)

im = None
for ax, method in zip(axes, METHODS):
    dfm = sub[sub["method"] == method]
    piv = dfm.pivot(index="sig_f2", columns="batch_f2", values="delta_auc").reindex(
        index=sig_vals, columns=bat_vals)
    im = ax.imshow(piv.values, aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(bat_vals)))
    ax.set_xticklabels([f"{v:.2f}" for v in bat_vals], fontsize=6.5, rotation=45)
    ax.set_title(TITLES[method], fontsize=6.8)
    ax.set_xlabel("cohort\ntarget $R^2$", fontsize=6.5)

    ri = sig_vals.index(REAL_SIG)
    ci = bat_vals.index(REAL_BAT)
    rect = mpatches.Rectangle((ci - 0.5, ri - 0.5), 1, 1, fill=False,
                               edgecolor="black", linewidth=1.6, zorder=5)
    ax.add_patch(rect)

axes[0].set_yticks(range(len(sig_vals)))
axes[0].set_yticklabels([f"{v:.3f}" for v in sig_vals], fontsize=6.5)
axes[0].set_ylabel("signal target $R^2$", fontsize=6.8)

cbar = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.015, aspect=18, extend="both")
cbar.set_label(r"$\Delta$AUC vs. single-cohort ceiling", fontsize=7)
cbar.ax.tick_params(labelsize=7)

# Footnote for the percentile_norm caveat (keeps the panel title short).
fig.text(0.5, 0.995, "*percentile normalization assumes the same signal "
          "direction in every cohort", fontsize=6.4, ha="center", va="top",
          style="italic", color="#555555")

save(fig, "results/revision/figures/fig7_simulation_grid_v3")

calib = pd.read_csv("results/ml/simulation/calibration_deltas.tsv", sep="\t")
calib_007 = calib[(calib["mode"] == "signal") & (calib["target_r2"] == 0.007)].iloc[0]
calib_027 = calib[(calib["mode"] == "signal") & (calib["target_r2"] == 0.027)].iloc[0]
print(f"\nsignal target R2=0.007 calibration: injected delta={calib_007['delta']:.3e} "
      f"(numerically zero), verified_r2={calib_007['verified_r2']:.5f}  -- NULL MODEL, "
      f"matches real n=118 pooled response R2.")
print(f"signal target R2=0.027 calibration (NOT the outlined cell): injected delta={calib_027['delta']:.4f}, "
      f"verified_r2={calib_027['verified_r2']:.5f} -- matches Cohort-1-alone (n=39) response R2.")

real_cell = sub[(sub["sig_f2"] == REAL_SIG) & (sub["batch_f2"] == REAL_BAT) & (sub["method"].isin(METHODS))]
print(f"\nReal operating point cell (signal target R2={REAL_SIG}, cohort target R2={REAL_BAT}) "
      "delta_auc by method (now including percentile_norm):")
print(real_cell[["method", "delta_auc", "mean_auc", "ceiling_auc"]].to_string(index=False))
print(f"Color scale: vmin={vmin:.4f}  vmax={vmax:.4f}")
