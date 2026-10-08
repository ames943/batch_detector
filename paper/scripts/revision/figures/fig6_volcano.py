#!/usr/bin/env python3
"""
Fig 6 — Volcano plot of the Hedges'-g meta-analysis.
Source: results/revision/task6_meta_analysis_hedges_g.tsv, task6_summary.txt (no recompute).
"""
import sys
sys.path.insert(0, "scripts/revision/figures")
import numpy as np
import pandas as pd
from style import apply, save, RESPONSE, CHANCE

apply()
import matplotlib.pyplot as plt

df = pd.read_csv("results/revision/task6_meta_analysis_hedges_g.tsv", sep="\t")
n_genera = len(df)
n_sig = int(df["fdr_sig"].sum())
mean_i2 = df["I2"].mean()

x = df["pooled_hedges_g"].values
y = -np.log10(df["q_value"].values)

fig, ax = plt.subplots(figsize=(3.5, 3.4))
fig.subplots_adjust(left=0.16, right=0.96, top=0.95, bottom=0.14)

ax.scatter(x, y, s=10, color=RESPONSE, alpha=0.65, edgecolor="none", zorder=3)

thresh = -np.log10(0.05)
ax.axhline(thresh, color=CHANCE, linestyle="--", linewidth=1.0, zorder=2)
ax.text(x.max() * 0.98, thresh + 0.04, "q=0.05", fontsize=7, color=CHANCE, ha="right", va="bottom")

# All 5 top hits sit far below the significance line and close together in x
# (max -log10(q) ~ 0.06) -- stagger label heights + leader lines so they don't
# overlap each other, and cap the y-axis just above the threshold line instead
# of wasting most of the panel on empty space between the data and q=0.05.
top5 = df.nsmallest(5, "q_value").sort_values("pooled_hedges_g").reset_index(drop=True)
label_ys = np.linspace(0.45, 0.95, len(top5))  # stacked label row heights (data units)
for (_, row), ly in zip(top5.iterrows(), label_ys):
    xi, yi = row["pooled_hedges_g"], -np.log10(row["q_value"])
    ha = "left" if xi >= 0 else "right"
    lx = xi + (0.03 if xi >= 0 else -0.03)
    ax.scatter([xi], [yi], s=16, facecolor="none", edgecolor="black", linewidth=0.8, zorder=4)
    ax.annotate(row["genus"], xy=(xi, yi), xytext=(lx, ly),
                textcoords="data", ha=ha, va="center", fontsize=7, style="italic",
                arrowprops=dict(arrowstyle="-", color="black", linewidth=0.5, alpha=0.7))

ax.set_xlabel("Pooled Hedges' $g$ (R vs NR, random-effects)")
ax.set_ylabel(r"$-\log_{10}(q)$")
ax.set_ylim(0, thresh * 1.25)
ax.set_xlim(x.min() - 0.15, x.max() + 0.45)
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, color="#dddddd", linewidth=0.4, zorder=0)
ax.set_axisbelow(True)

ax.text(0.02, 0.97, f"{n_genera} genera, {n_sig} FDR-significant\nmean $I^2$ = {mean_i2:.1f}%",
        transform=ax.transAxes, fontsize=7, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#bbbbbb", alpha=0.92))

save(fig, "results/revision/figures/fig6_volcano")

print(f"n_genera={n_genera}  n_sig={n_sig}  mean_I2={mean_i2:.2f}")
print("Top 5 by smallest q:")
print(top5[["genus", "pooled_hedges_g", "q_value", "I2"]].to_string(index=False))
print(f"max -log10(q) observed = {y.max():.3f}  (threshold line at {thresh:.3f}); "
      f"{'no point crosses q=0.05' if y.max() < thresh else 'WARNING: a point crosses q=0.05'}")
