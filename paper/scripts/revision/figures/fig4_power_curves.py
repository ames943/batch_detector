#!/usr/bin/env python3
"""
Fig 4 — Power curves over hypothetical effect sizes, with the C4 95% upper bound.
Source: results/revision/task8_power_curve.tsv, task8_c4_upper_bound.tsv (no recompute).
No calibrated curve exists on disk for delta_R2=0.0078 (only 0.0025/0.005/0.01/0.02
were run) -- only the single (n,power) marker for that bound is
drawn, not an invented curve.
"""
import sys
sys.path.insert(0, "scripts/revision/figures")
import pandas as pd
from style import apply, save, OKABE_ITO

apply()
import matplotlib.pyplot as plt

curve = pd.read_csv("results/revision/task8_power_curve.tsv", sep="\t")
bound = pd.read_csv("results/revision/task8_c4_upper_bound.tsv", sep="\t").iloc[0]

deltas = sorted(curve["target_delta_R2"].unique())
colors = {d: OKABE_ITO[i] for i, d in enumerate(deltas)}
markers = {d: m for d, m in zip(deltas, ["o", "s", "^", "D"])}

fig, ax = plt.subplots(figsize=(3.5, 3.2))
fig.subplots_adjust(left=0.16, right=0.97, top=0.96, bottom=0.15)

for d in deltas:
    sub = curve[curve["target_delta_R2"] == d].sort_values("n")
    ax.plot(sub["n"], sub["power"], marker=markers[d], markersize=3.6,
            color=colors[d], linewidth=1.1, label=f"$\\Delta R^2$={d}")

ax.axhline(0.80, color="gray", linestyle="--", linewidth=1.0, zorder=1)
ax.text(1050, 0.815, "80% power", fontsize=7, color="gray", ha="right")

ax.axvline(165, color="black", linestyle=":", linewidth=1.0, zorder=1)
ax.text(165, 1.03, "C4 (n=165)", fontsize=7, ha="center", va="bottom")

n_bound = float(bound["n_for_80pct_power_at_bound"])
ax.scatter([n_bound], [0.80], marker="*", s=90, color="black", zorder=5,
           label=f"C4 95% UB $\\Delta R^2$={bound['dirichlet_boot_95pct_upper_bound_delta_R2']:.4f}")
ax.annotate(f"n≈{n_bound:.0f}\nat $\\Delta R^2$=0.0078", (n_bound, 0.80),
            textcoords="offset points", xytext=(8, -28), fontsize=7, ha="left",
            arrowprops=dict(arrowstyle="-", color="black", linewidth=0.6))

ax.set_xscale("log")
ax.set_xlim(30, 1300)
ax.set_ylim(0, 1.08)
ax.set_xlabel("$n$ (log scale)")
ax.set_ylabel("Power (fraction $p$<0.05)")
ax.spines[["top", "right"]].set_visible(False)
ax.yaxis.grid(True, color="#dddddd", linewidth=0.4, zorder=0)
ax.set_axisbelow(True)
ax.legend(fontsize=7, loc="lower right", frameon=False, handlelength=1.6)

save(fig, "results/revision/figures/fig4_power_curves")

print("No calibrated power curve exists on disk for delta_R2=0.0078 (C4 95% UB) --")
print("only the (n,power)=(%.1f, 0.80) marker is drawn." % n_bound)
print(curve.to_string(index=False))
print(bound.to_string())
