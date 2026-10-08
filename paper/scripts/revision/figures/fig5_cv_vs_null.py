#!/usr/bin/env python3
"""
Fig 5 — Repeated-CV per-fold AUC distributions vs. the permutation null.
Source: results/revision/task5_{dataset}_fold_aucs.tsv,
        task5_{dataset}_permutation_aucs.tsv, task5_{dataset}_permutation_summary.tsv
(no recompute).
"""
import sys
sys.path.insert(0, "scripts/revision/figures")
import numpy as np
import pandas as pd
from style import apply, save, panel_letter, RESPONSE, CHANCE

apply()
import matplotlib.pyplot as plt

DATASETS = [("c1", "C1 (n=39)"), ("n118", "n=118"), ("n283", "n=283"), ("c4", "C4 (n=165)")]
MODELS = [("en", "ElasticNet"), ("rf", "RandomForest"), ("xgb", "XGBoost")]

fig, axes = plt.subplots(len(DATASETS), len(MODELS), figsize=(7.1, 8.3), sharey=True)
fig.subplots_adjust(left=0.07, right=0.98, top=0.88, bottom=0.05, hspace=0.55, wspace=0.15)

rng = np.random.default_rng(0)

for ri, (dkey, dlabel) in enumerate(DATASETS):
    fold_df = pd.read_csv(f"results/revision/task5_{dkey}_fold_aucs.tsv", sep="\t")
    perm_df = pd.read_csv(f"results/revision/task5_{dkey}_permutation_aucs.tsv", sep="\t")
    summ_df = pd.read_csv(f"results/revision/task5_{dkey}_permutation_summary.tsv", sep="\t").set_index("model")

    for ci, (mkey, mlabel) in enumerate(MODELS):
        ax = axes[ri, ci]
        fold_vals = fold_df[f"auc_{mkey}"].values
        null_vals = perm_df[mkey].values
        obs_mean = float(summ_df.loc[mkey, "observed_mean_auc"])
        p_val = float(summ_df.loc[mkey, "p_value"])

        # strip plot of the 50 per-fold AUCs at x=0
        jitter = (rng.random(len(fold_vals)) - 0.5) * 0.28
        ax.scatter(jitter, fold_vals, s=5, color=RESPONSE, alpha=0.6, zorder=3,
                   label="per-fold AUC" if (ri == 0 and ci == 0) else None)
        ax.scatter([0], [obs_mean], marker="D", s=28, color="black", zorder=5,
                   label="observed mean" if (ri == 0 and ci == 0) else None)

        # violin of the 100 permutation-null MEAN AUCs at x=1
        vp = ax.violinplot(null_vals, positions=[1], widths=0.6, showextrema=False)
        for b in vp["bodies"]:
            b.set_facecolor(CHANCE)
            b.set_alpha(0.55)
            b.set_edgecolor(CHANCE)

        ax.axhline(0.5, color="gray", linestyle=":", linewidth=0.8, zorder=1)
        ax.set_xlim(-0.5, 1.5)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["folds", "null"], fontsize=7)
        ax.set_ylim(0.05, 1.0)

        ax.annotate(f"p={p_val:.3f}", (0.5, 0.98), xycoords="axes fraction",
                    ha="center", va="top", fontsize=7,
                    fontweight="bold" if p_val < 0.05 else "normal")

        if ri == 0:
            ax.set_title(mlabel, fontsize=7.5, pad=6)
        if ci == 0:
            ax.set_ylabel(dlabel, fontsize=7.2)
        ax.tick_params(axis="y", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)

fig.text(0.5, 0.99, "Per-fold AUC (repeated 5x10-fold CV) vs. permutation-null mean AUC (N=100)",
          fontsize=8, ha="center", va="top")
fig.legend(loc="upper center", ncol=2, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.955))

save(fig, "results/revision/figures/fig5_cv_vs_null")

print("Observed mean AUC / p by dataset x model:")
for dkey, dlabel in DATASETS:
    summ_df = pd.read_csv(f"results/revision/task5_{dkey}_permutation_summary.tsv", sep="\t")
    for _, row in summ_df.iterrows():
        print(f"  {dlabel:>10} {row['model']:>4}: obs={row['observed_mean_auc']:.4f}  p={row['p_value']:.4f}")
