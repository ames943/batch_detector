#!/usr/bin/env python3
"""
Fig 4 (v2) — Power curves, recomputed in C4's OWN feature space.
Source: results/revision/task8_c4_power_v2.tsv (no recompute; see
scripts/revision/task8b_c4_power_v2.py for how these were generated, and
results/revision/SUMMARY.md "Power check v2" for the full writeup of the bug
this fixes).

Unlike the v1 figure, a delta_R2=0.0078 curve IS drawn here -- it was
actually computed this time (v1 only had 0.0025/0.005/0.01/0.02). The star
marking C4's 95% upper-bound effect is placed at (n=165, power=0.995), an
ACTUAL COMPUTED grid point (v1 placed it at an interpolated n on the wrong
curve).
"""
import sys
sys.path.insert(0, "scripts/revision/figures")
import pandas as pd
from style import apply, save, OKABE_ITO

apply()
import matplotlib.pyplot as plt

curve = pd.read_csv("results/revision/task8_c4_power_v2.tsv", sep="\t")

deltas = sorted(curve["target_delta_R2"].unique())
colors = {d: OKABE_ITO[i] for i, d in enumerate(deltas)}
markers = {d: m for d, m in zip(deltas, ["o", "s", "^", "D"])}
labels = {0.0025: "$\\Delta R^2$=0.0025", 0.005: "$\\Delta R^2$=0.005",
          0.0078: "$\\Delta R^2$=0.0078 (C4 95% UB)", 0.01: "$\\Delta R^2$=0.01"}

fig, ax = plt.subplots(figsize=(3.5, 3.2))
fig.subplots_adjust(left=0.16, right=0.97, top=0.88, bottom=0.15)

for d in deltas:
    sub = curve[curve["target_delta_R2"] == d].sort_values("n")
    lw = 1.6 if abs(d - 0.0078) < 1e-9 else 1.1
    ax.plot(sub["n"], sub["power"], marker=markers[d], markersize=3.8,
            color=colors[d], linewidth=lw, label=labels[d], zorder=4 if abs(d-0.0078) < 1e-9 else 3)

ax.axhline(0.80, color="gray", linestyle="--", linewidth=1.0, zorder=1)
ax.text(390, 0.815, "80% power", fontsize=7, color="gray", ha="right")

ax.axvline(165, color="black", linestyle=":", linewidth=1.0, zorder=1)
# Offset to the right of the dotted line, at the top (was sitting directly on
# top of the line, and a lower placement collided with the legend).
ax.annotate("C4\n(n=165)", xy=(165, 1.08), xytext=(230, 1.08),
            fontsize=7, ha="left", va="center", annotation_clip=False,
            arrowprops=dict(arrowstyle="-", color="black", linewidth=0.6))

# Star at an ACTUALLY COMPUTED grid point: n=165, delta_R2=0.0078
row = curve[(curve["n"] == 165) & (curve["target_delta_R2"] == 0.0078)].iloc[0]
ax.scatter([165], [row["power"]], marker="*", s=110, color="black", zorder=6,
           edgecolor="white", linewidth=0.5)
ax.annotate(f"power={row['power']:.3f}\nat C4's own n and\n95% UB effect",
            (165, row["power"]), textcoords="offset points", xytext=(-95, -8),
            fontsize=7, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color="black", linewidth=0.6))

ax.set_xscale("log")
ax.set_xlim(60, 500)
ax.set_ylim(0, 1.0)
ax.set_xlabel("$n$ (log scale)")
ax.set_ylabel("Power (fraction $p$<0.05)")
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, color="#dddddd", linewidth=0.4, zorder=0)
ax.set_axisbelow(True)
ax.legend(fontsize=7, loc="lower right", frameon=False, handlelength=1.6)

save(fig, "results/revision/figures/fig4_power_curves_v2")

print(curve.to_string(index=False))
print(f"\nStar: n=165, delta_R2=0.0078, power={row['power']:.4f} (actual computed grid point)")
