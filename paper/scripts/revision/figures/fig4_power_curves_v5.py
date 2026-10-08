#!/usr/bin/env python3
"""
Fig 4 (v5) — Power curves, recomputed in C4's OWN feature space.
Source: results/revision/task8_c4_power_v2.tsv (UNCHANGED -- layout-only fix,
no data or numbers changed from v2/v3/v4).

v5 fixes vs. v4:
  1. The "0.995" text label next to the star removed entirely (the paper
     caption already explains the star).
  2. "80% power" label moved to the right end of the dashed line (near
     n=400, just above it), clear of every curve there (all four curves sit
     at ~0.98-1.0 at n=400, well above the label's y=0.815).
"""
import sys
sys.path.insert(0, "scripts/revision/figures")
import pandas as pd
from matplotlib.ticker import NullLocator, FixedLocator, FixedFormatter
from style import apply, save, OKABE_ITO

apply()
import matplotlib.pyplot as plt

curve = pd.read_csv("results/revision/task8_c4_power_v2.tsv", sep="\t")

deltas = sorted(curve["target_delta_R2"].unique())
colors = {d: OKABE_ITO[i] for i, d in enumerate(deltas)}
markers = {d: m for d, m in zip(deltas, ["o", "s", "^", "D"])}
labels = {0.0025: "$\\Delta R^2$=0.0025", 0.005: "$\\Delta R^2$=0.005",
          0.0078: "$\\Delta R^2$=0.0078 (C4 95% UB)", 0.01: "$\\Delta R^2$=0.01"}

fig, ax = plt.subplots(figsize=(3.5, 3.9))
fig.subplots_adjust(left=0.16, right=0.97, top=0.88, bottom=0.30)

for d in deltas:
    sub = curve[curve["target_delta_R2"] == d].sort_values("n")
    lw = 1.6 if abs(d - 0.0078) < 1e-9 else 1.1
    ax.plot(sub["n"], sub["power"], marker=markers[d], markersize=3.8,
            color=colors[d], linewidth=lw, label=labels[d],
            zorder=4 if abs(d - 0.0078) < 1e-9 else 3)

ax.axhline(0.80, color="gray", linestyle="--", linewidth=1.0, zorder=1)
# Label at the right end of the dashed line, near n=400, just above it --
# every curve is already up near 0.98-1.0 by n=400, well clear of y=0.815.
ax.text(400, 0.815, "80% power", fontsize=7, color="gray", ha="right", va="bottom")

ax.axvline(165, color="black", linestyle=":", linewidth=1.0, zorder=1)
ax.annotate("C4\n(n=165)", xy=(165, 1.08), xytext=(180, 1.08),
            fontsize=7, ha="left", va="center", annotation_clip=False)

# Star at an ACTUALLY COMPUTED grid point: n=165, delta_R2=0.0078
row = curve[(curve["n"] == 165) & (curve["target_delta_R2"] == 0.0078)].iloc[0]
ax.scatter([165], [row["power"]], marker="*", s=110, color="black", zorder=6,
           edgecolor="white", linewidth=0.5)
# No text label -- the paper caption explains the star.

ax.set_xscale("log")
ax.set_xlim(60, 500)
ax.set_ylim(0, 1.0)
tick_ns = [80, 131, 165, 250, 400]
ax.xaxis.set_major_locator(FixedLocator(tick_ns))
ax.xaxis.set_major_formatter(FixedFormatter([str(n) for n in tick_ns]))
ax.xaxis.set_minor_locator(NullLocator())
ax.tick_params(axis="x", labelsize=6.5, rotation=0)
ax.set_xlabel("$n$ (log scale)")
ax.set_ylabel("Power (fraction $p$<0.05)")
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, color="#dddddd", linewidth=0.4, zorder=0)
ax.set_axisbelow(True)

# Legend fully outside the axes, below the x-axis label, 2 columns.
ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.22),
          ncol=2, frameon=False, handlelength=1.6, columnspacing=1.2,
          labelspacing=0.4)

save(fig, "results/revision/final/fig4_power_curves_v5")

print(curve.to_string(index=False))
print(f"\nStar: n=165, delta_R2=0.0078, power={row['power']:.4f} (actual computed grid point, unlabeled)")
print("No numbers changed from v2/v3/v4 -- layout-only fix.")
