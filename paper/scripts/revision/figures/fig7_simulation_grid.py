#!/usr/bin/env python3
"""
Fig 7 — Phase-2 simulation grid: Delta_AUC vs. single-cohort ceiling, nc=3, npc=40.
Source: results/ml/simulation/grid_results.tsv (no recompute).
percentile_norm is skipped: its simulated gains are a known generator artifact
(the simulation's consistent-signal-direction assumption across cohorts, which
Phase 0.5's real-data cross-cohort holdout showed does not hold) -- see project
notes / Phase 2 follow-up.
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

METHODS = ["none", "mean_centering", "location_scale", "quantile_mapping", "cohort_covariate"]
TITLES = {"none": "none", "mean_centering": "mean-\ncentering", "location_scale": "location-\nscale",
          "quantile_mapping": "quantile\nmapping", "cohort_covariate": "cohort\ncovariate"}

sig_vals = sorted(sub["sig_f2"].unique(), reverse=True)  # high->low, top->bottom
bat_vals = sorted(sub["batch_f2"].unique())               # low->high, left->right

vmax = sub[sub["method"].isin(METHODS)]["delta_auc"].abs().max()
vmin = -vmax

REAL_SIG, REAL_BAT = 0.027, 0.08

fig, axes = plt.subplots(1, len(METHODS), figsize=(7.1, 2.3), sharey=True)
fig.subplots_adjust(left=0.09, right=0.90, top=0.80, bottom=0.30, wspace=0.12)

im = None
for ax, method in zip(axes, METHODS):
    dfm = sub[sub["method"] == method]
    piv = dfm.pivot(index="sig_f2", columns="batch_f2", values="delta_auc").reindex(
        index=sig_vals, columns=bat_vals)
    im = ax.imshow(piv.values, aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(bat_vals)))
    ax.set_xticklabels([f"{v:.2f}" for v in bat_vals], fontsize=7, rotation=45)
    ax.set_title(TITLES[method], fontsize=7)
    ax.set_xlabel("cohort $f^2$", fontsize=7)

    ri = sig_vals.index(REAL_SIG)
    ci = bat_vals.index(REAL_BAT)
    rect = mpatches.Rectangle((ci - 0.5, ri - 0.5), 1, 1, fill=False,
                               edgecolor="black", linewidth=1.6, zorder=5)
    ax.add_patch(rect)

axes[0].set_yticks(range(len(sig_vals)))
axes[0].set_yticklabels([f"{v:.3f}" for v in sig_vals], fontsize=7)
axes[0].set_ylabel("signal $f^2$", fontsize=7)

cbar = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02, aspect=18)
cbar.set_label(r"$\Delta$AUC vs. single-cohort ceiling", fontsize=7)
cbar.ax.tick_params(labelsize=7)

fig.suptitle("Simulation grid (nc=3, n/cohort=40): black outline = real operating point "
             f"(signal $f^2$={REAL_SIG}, cohort $f^2$={REAL_BAT})", fontsize=7.2, y=0.98)

save(fig, "results/revision/figures/fig7_simulation_grid")

real_cell = sub[(sub["sig_f2"] == REAL_SIG) & (sub["batch_f2"] == REAL_BAT) & (sub["method"].isin(METHODS))]
print("percentile_norm excluded (known generator artifact -- see docstring).")
print("Real operating point cell (sig_f2=0.027, batch_f2=0.08) delta_auc by method:")
print(real_cell[["method", "delta_auc", "mean_auc", "ceiling_auc"]].to_string(index=False))
print(f"Color scale: vmin={vmin:.4f}  vmax={vmax:.4f}")
