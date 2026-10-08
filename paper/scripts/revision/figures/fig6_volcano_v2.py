#!/usr/bin/env python3
"""
Fig 6 (v2) — Volcano plot of the Hedges'-g meta-analysis.
Source: results/revision/task6_meta_analysis_hedges_g.tsv, task6_summary.txt
(UNCHANGED -- layout-only fix, no data or numbers changed from v1).

v2 fixes vs. v1:
  1. The 4 right-side labels (Pseudomonas, Ligilactobacillus,
     Phascolarctobacterium, Wujia) were stacked from y=0.45 to 0.95 (upper
     half), which crossed leader lines. Re-laid-out strictly ordered by each
     point's x position, in the LOWER half of the plot, close to their
     points, with short leader lines -- both x and y increase together in
     the same point order, so lines cannot cross by construction.
  2. Base marker size increased (10 -> 20) so the near-zero cluster of 125
     points reads as a scatter, not a flat smear. Axis scale/limits
     unchanged -- still no broken axis, same y-range showing the real gap to
     q=0.05.
  3. Stats box and "q=0.05" label repositioned slightly so they don't overlap
     each other (stats box lowered a touch; q=0.05 label kept at the far
     right of the dashed line as before -- confirmed clear of the box).
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

ax.scatter(x, y, s=20, color=RESPONSE, alpha=0.55, edgecolor="none", zorder=3)

thresh = -np.log10(0.05)
ax.axhline(thresh, color=CHANCE, linestyle="--", linewidth=1.0, zorder=2)
ax.text(x.max() * 0.98, thresh + 0.04, "q=0.05", fontsize=7, color=CHANCE, ha="right", va="bottom")

top5 = df.nsmallest(5, "q_value").copy()
dialister = top5[top5["pooled_hedges_g"] < 0]
right4 = top5[top5["pooled_hedges_g"] >= 0].sort_values("pooled_hedges_g").reset_index(drop=True)

# The 4 right-side points sit almost on top of each other in x (0.216-0.329;
# 3 of them share the exact same q-value/y-height). Diagonal fans from
# individually-offset label positions still read as crossing at this
# proximity -- any two lines with overlapping x-ranges but different slopes
# look crossed even when they technically aren't. The robust fix (same idea
# as non-crossing dendrogram tip labels): anchor all 4 labels at a SINGLE
# fixed x to the right of the cluster, stacked top-to-bottom in the SAME
# left-to-right order as their points (rightmost point <-> topmost label).
# Spokes from a common vertical line to monotonically-ordered points cannot
# cross by construction. Kept in the lower half (max label y=0.55).
LABEL_X = 0.42
label_ys = [0.38, 0.28, 0.18, 0.10]  # top->bottom = right->left point order
for (_, row), ly in zip(right4.iloc[::-1].iterrows(), label_ys):
    xi, yi = row["pooled_hedges_g"], -np.log10(row["q_value"])
    ax.scatter([xi], [yi], s=22, facecolor="none", edgecolor="black", linewidth=0.9, zorder=4)
    ax.annotate(row["genus"], xy=(xi, yi), xytext=(LABEL_X, ly),
                textcoords="data", ha="left", va="center", fontsize=7, style="italic",
                arrowprops=dict(arrowstyle="-", color="black", linewidth=0.5, alpha=0.75))

# Dialister: isolated on the left, own short label, also lower half.
for _, row in dialister.iterrows():
    xi, yi = row["pooled_hedges_g"], -np.log10(row["q_value"])
    ax.scatter([xi], [yi], s=22, facecolor="none", edgecolor="black", linewidth=0.9, zorder=4)
    ax.annotate(row["genus"], xy=(xi, yi), xytext=(xi - 0.02, 0.30),
                textcoords="data", ha="right", va="center", fontsize=7, style="italic",
                arrowprops=dict(arrowstyle="-", color="black", linewidth=0.5, alpha=0.75))

ax.set_xlabel("Pooled Hedges' $g$ (R vs NR, random-effects)")
ax.set_ylabel(r"$-\log_{10}(q)$")
ax.set_ylim(0, thresh * 1.25)
ax.set_xlim(x.min() - 0.15, LABEL_X + 0.70)
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, color="#dddddd", linewidth=0.4, zorder=0)
ax.set_axisbelow(True)

ax.text(0.02, 0.97, f"{n_genera} genera, {n_sig} FDR-significant\nmean $I^2$ = {mean_i2:.1f}%",
        transform=ax.transAxes, fontsize=7, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#bbbbbb", alpha=0.92))

save(fig, "results/revision/final/fig6_volcano_v2")

print(f"n_genera={n_genera}  n_sig={n_sig}  mean_I2={mean_i2:.2f}")
print("Top 5 by smallest q:")
print(top5.sort_values("q_value")[["genus", "pooled_hedges_g", "q_value", "I2"]].to_string(index=False))
print(f"max -log10(q) observed = {y.max():.3f}  (threshold line at {thresh:.3f}); "
      f"{'no point crosses q=0.05' if y.max() < thresh else 'WARNING: a point crosses q=0.05'}")
print("No numbers changed from v1 -- layout-only fix.")
