#!/usr/bin/env python3
"""
Fig 1 (fixed) — chance-calibrated dose-response.
Source: results/revision/task1_chance_calibrated_permanova.tsv (no recompute).
Fixes vs. the first draft: legend moved below the panels (was overlapping the
"10.6" data label); both chance curves now separately styled and legend-labeled;
Panel B labels all placed above their points with widened y-limits (was
colliding with the x-axis / clipped at the bottom).
"""
import sys
sys.path.insert(0, "scripts/revision/figures")
import pandas as pd
from matplotlib.ticker import MultipleLocator
from style import apply, save, RESPONSE, COHORT, CHANCE

apply()
import matplotlib.pyplot as plt

df = pd.read_csv("results/revision/task1_chance_calibrated_permanova.tsv", sep="\t")
NS = [39, 79, 118, 283]
POS = {n: i for i, n in enumerate(NS)}
NAME = {39: "genus_n39", 79: "genus_n79", 118: "genus_n118", 283: "genus_n283"}


def get(n, factor, col):
    sub = df[(df["dataset"] == NAME[n]) & (df["factor"] == factor)]
    return None if sub.empty else float(sub.iloc[0][col])


resp_R2 = {n: get(n, "response", "R2") for n in NS}
resp_E0 = {n: get(n, "response", "E0") for n in NS}
resp_dR2 = {n: get(n, "response", "delta_R2") for n in NS}
coh_R2 = {n: get(n, "cohort", "R2") for n in NS if get(n, "cohort", "R2") is not None}
coh_E0 = {n: get(n, "cohort", "E0") for n in NS if get(n, "cohort", "R2") is not None}
coh_dR2 = {n: get(n, "cohort", "delta_R2") for n in NS if get(n, "cohort", "R2") is not None}

XT = [f"{n}\n{c}" for n, c in zip(NS, ["C1", "+C2", "+C3", "+C4"])]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.1, 2.6))
fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.30, wspace=0.30)

# ══════════════════════════════════ Panel A ══════════════════════════════════
rx = [POS[n] for n in NS]
ry = [resp_R2[n] * 100 for n in NS]
chance_r = [resp_E0[n] * 100 for n in NS]
bx = [POS[n] for n in NS if n in coh_R2]
by = [coh_R2[n] * 100 for n in NS if n in coh_R2]
chance_b = [coh_E0[n] * 100 for n in NS if n in coh_E0]

lcr, = axA.plot(rx, chance_r, color=CHANCE, linestyle=":", linewidth=1.3, zorder=2,
                label="Chance (response), $1/(n-1)$")
lcb, = axA.plot(bx, chance_b, color=CHANCE, linestyle="-.", linewidth=1.3, zorder=2,
                label="Chance (cohort), $(k-1)/(n-1)$")
lr, = axA.plot(rx, ry, color=RESPONSE, marker="o", markersize=3.6, markerfacecolor=RESPONSE,
               markeredgecolor=RESPONSE, zorder=4, clip_on=False, label="Response $R^2$")
lb, = axA.plot(bx, by, color=COHORT, marker="s", markersize=3.8,
               linestyle=(0, (3.5, 2)), markerfacecolor="white",
               markeredgecolor=COHORT, markeredgewidth=1.0, zorder=4,
               clip_on=False, label="Cohort $R^2$")

for x, y in zip(rx, ry):
    axA.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 5),
                 ha="center", va="bottom", fontsize=7, color=RESPONSE)
for x, y in zip(bx, by):
    axA.annotate(f"{y:.1f}", (x, y), textcoords="offset points", xytext=(0, -9),
                 ha="center", va="top", fontsize=7, color=COHORT)

axA.set_ylabel("PERMANOVA $R^2$ (%)", labelpad=2)
axA.set_ylim(-0.6, 13)
axA.yaxis.set_major_locator(MultipleLocator(4))
from style import panel_letter
panel_letter(axA, "A")

# ══════════════════════════════════ Panel B ══════════════════════════════════
rdx = [POS[n] for n in NS]
rdy = [resp_dR2[n] * 100 for n in NS]
bdx = [POS[n] for n in NS if n in coh_dR2]
bdy = [coh_dR2[n] * 100 for n in NS if n in coh_dR2]

axB.axhline(0, color="black", linewidth=0.7, zorder=2)
axB.plot(rdx, rdy, color=RESPONSE, marker="o", markersize=3.6, markerfacecolor=RESPONSE,
         markeredgecolor=RESPONSE, zorder=4, clip_on=False)
axB.plot(bdx, bdy, color=COHORT, marker="s", markersize=3.8,
         linestyle=(0, (3.5, 2)), markerfacecolor="white",
         markeredgecolor=COHORT, markeredgewidth=1.0, zorder=4, clip_on=False)

# All response labels placed ABOVE their point (values all near zero; keeps
# them off the x-axis and away from the cohort labels, which sit much higher).
for x, y in zip(rdx, rdy):
    axB.annotate(f"{y:+.2f}", (x, y), textcoords="offset points", xytext=(0, 8),
                 ha="center", va="bottom", fontsize=7, color=RESPONSE)
for x, y in zip(bdx, bdy):
    axB.annotate(f"{y:+.1f}", (x, y), textcoords="offset points", xytext=(0, 6),
                 ha="center", va="bottom", fontsize=7, color=COHORT)

axB.set_ylabel(r"$\Delta R^2 = R^2 - E_0$ (%)", labelpad=2)
axB.set_ylim(-1.5, 11.5)
axB.yaxis.set_major_locator(MultipleLocator(3))
panel_letter(axB, "B")

for ax in (axA, axB):
    ax.set_xlim(-0.35, 3.35)
    ax.set_xticks(range(4))
    ax.set_xticklabels(XT)
    ax.tick_params(axis="x", pad=1.5)
    ax.set_xlabel("cumulative $n$", labelpad=1)
    ax.yaxis.grid(True, color="#dddddd", linewidth=0.4, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=2)

# Fig-level legend BELOW both panels (was overlapping data inside Panel A)
handles = [lr, lb, lcr, lcb]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
           bbox_to_anchor=(0.5, -0.02), fontsize=7, handlelength=1.8,
           columnspacing=1.2, handletextpad=0.4)
fig.subplots_adjust(bottom=0.40)

save(fig, "results/revision/figures/fig1_dose_response")

print("\nValues plotted (source: task1_chance_calibrated_permanova.tsv):")
for n in NS:
    print(f"  n={n}: response R2={resp_R2[n]*100:.3f}%  E0={resp_E0[n]*100:.3f}%  "
          f"dR2={resp_dR2[n]*100:+.3f}%   cohort R2={coh_R2.get(n)}  dR2={coh_dR2.get(n)}")
