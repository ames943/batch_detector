#!/usr/bin/env python3
"""
Fig 7 (v2) — Phase-2 simulation grid: Delta_AUC vs. single-cohort ceiling, nc=3, npc=40.
Source: results/ml/simulation/grid_results.tsv,
        results/ml/simulation/calibration_deltas.tsv (no recompute).
percentile_norm is skipped: its simulated gains are a known generator artifact
(the simulation's consistent-signal-direction assumption across cohorts, which
Phase 0.5's real-data cross-cohort holdout showed does not hold) -- see project
notes / Phase 2 follow-up.

v2: outline moved from (signal f2=0.027, cohort f2=0.08) to
(signal f2=0.007, cohort f2=0.08) -- the actual real operating point. Per
results/ml/simulation/calibration_deltas.tsv and scripts/phase2_simulation.py
(REAL_SIG_F2_N118=0.007, comment "sig_f2=0.007 row = NULL MODEL (delta=0 from
calibration)"): signal f2=0.007 was calibrated to an injected shift of
delta=4.66e-09 (numerically zero -- no injected signal at all), chosen because
its simulated R2 (0.00851, verified_r2 column) matches the REAL measured n=118
response R2 (~0.007-0.0085). signal f2=0.027 is a DIFFERENT, larger calibrated
point (delta=2.987, verified_r2=0.0265) matching Cohort-1-alone's (n=39) response
R2 -- not the pooled n=118 operating point this figure is meant to mark.
Suptitle removed (paper caption covers it).
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

REAL_SIG, REAL_BAT = 0.007, 0.08

fig, axes = plt.subplots(1, len(METHODS), figsize=(7.1, 2.3), sharey=True)
fig.subplots_adjust(left=0.09, right=0.90, top=0.86, bottom=0.30, wspace=0.12)

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

save(fig, "results/revision/figures/fig7_simulation_grid_v2")

calib = pd.read_csv("results/ml/simulation/calibration_deltas.tsv", sep="\t")
calib_007 = calib[(calib["mode"] == "signal") & (calib["target_r2"] == 0.007)].iloc[0]
calib_027 = calib[(calib["mode"] == "signal") & (calib["target_r2"] == 0.027)].iloc[0]
print("percentile_norm excluded (known generator artifact -- see docstring).")
print(f"\nsignal f2=0.007 calibration: injected delta={calib_007['delta']:.3e} "
      f"(numerically zero), verified_r2={calib_007['verified_r2']:.5f}  -- NULL MODEL, "
      f"matches real n=118 pooled response R2.")
print(f"signal f2=0.027 calibration (NOT the outlined cell): injected delta={calib_027['delta']:.4f}, "
      f"verified_r2={calib_027['verified_r2']:.5f} -- matches Cohort-1-alone (n=39) response R2.")

real_cell = sub[(sub["sig_f2"] == REAL_SIG) & (sub["batch_f2"] == REAL_BAT) & (sub["method"].isin(METHODS))]
print(f"\nReal operating point cell (sig_f2={REAL_SIG}, batch_f2={REAL_BAT}) delta_auc by method:")
print(real_cell[["method", "delta_auc", "mean_auc", "ceiling_auc"]].to_string(index=False))
print(f"Color scale: vmin={vmin:.4f}  vmax={vmax:.4f}")
