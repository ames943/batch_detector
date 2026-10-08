#!/usr/bin/env python3
"""
Task 10 — Regenerated Figure 1: chance-calibrated dose-response.

Panel A: observed response R2 (filled circles) and cohort R2 (open squares)
at n=39/79/118/283, with dotted chance curves 1/(n-1) [response] and
(k-1)/(n-1) [cohort]; every value labeled.
Panel B: delta_R2 for both factors, horizontal zero line.

Single-column width, vector PDF, embedded fonts (pdf.fonttype=42).
Reads results/revision/task1_chance_calibrated_permanova.tsv — nothing hardcoded.
"""
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

df = pd.read_csv("results/revision/task1_chance_calibrated_permanova.tsv", sep="\t")
NS = [39, 79, 118, 283]
POS = {n: i for i, n in enumerate(NS)}
NAME = {39: "genus_n39", 79: "genus_n79", 118: "genus_n118", 283: "genus_n283"}
N_COHORTS = {39: 1, 79: 2, 118: 3, 283: 4}


def get(n, factor, col):
    sub = df[(df["dataset"] == NAME[n]) & (df["factor"] == factor)]
    if sub.empty:
        return None
    return float(sub.iloc[0][col])


resp_R2 = {n: get(n, "response", "R2") for n in NS}
resp_E0 = {n: get(n, "response", "E0") for n in NS}
resp_dR2 = {n: get(n, "response", "delta_R2") for n in NS}
coh_R2 = {n: get(n, "cohort", "R2") for n in NS if get(n, "cohort", "R2") is not None}
coh_E0 = {n: get(n, "cohort", "E0") for n in NS if get(n, "cohort", "R2") is not None}
coh_dR2 = {n: get(n, "cohort", "delta_R2") for n in NS if get(n, "cohort", "R2") is not None}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "axes.labelsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "lines.linewidth": 1.1,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

INK = "#000000"       # response
BATCH_INK = "#6b6b6b"  # cohort/batch
CHANCE_INK = "#b0b0b0"
GRID = "#dddddd"
XT = [f"{n}\n{c}" for n, c in zip(NS, ["C1", "+C2", "+C3", "+C4"])]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(3.4, 2.0))
fig.subplots_adjust(left=0.13, right=0.985, top=0.88, bottom=0.30, wspace=0.42)

# ══════════════════════════════════════ Panel A: R2 (%) + chance curves ═════
rx = [POS[n] for n in NS]
ry = [resp_R2[n] * 100 for n in NS]
chance_r = [resp_E0[n] * 100 for n in NS]
bx = [POS[n] for n in NS if n in coh_R2]
by = [coh_R2[n] * 100 for n in NS if n in coh_R2]
chance_b = [coh_E0[n] * 100 for n in NS if n in coh_E0]

axA.plot(rx, chance_r, color=CHANCE_INK, linestyle=":", linewidth=1.0, zorder=2)
axA.plot(bx, chance_b, color=CHANCE_INK, linestyle=":", linewidth=1.0, zorder=2)

lr, = axA.plot(rx, ry, color=INK, marker="o", markersize=3.2, markerfacecolor=INK,
               markeredgecolor=INK, zorder=4, clip_on=False, label="Response $R^2$")
lb, = axA.plot(bx, by, color=BATCH_INK, marker="s", markersize=3.4,
               linestyle=(0, (3.5, 2)), markerfacecolor="white",
               markeredgecolor=BATCH_INK, markeredgewidth=0.9, zorder=4,
               clip_on=False, label="Cohort $R^2$")
lc, = axA.plot([], [], color=CHANCE_INK, linestyle=":", linewidth=1.0, label="Chance $E_0$")

axA.legend(handles=[lr, lb, lc], loc="upper right", frameon=False,
           handlelength=1.6, handletextpad=0.4, labelspacing=0.25,
           borderaxespad=0.0, fontsize=6.3)

for x, y in zip(rx, ry):
    axA.annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 4),
                 ha="center", va="bottom", fontsize=6.3, color=INK)
for x, y in zip(bx, by):
    axA.annotate(f"{y:.1f}", (x, y), textcoords="offset points", xytext=(0, -8),
                 ha="center", va="top", fontsize=6.3, color=BATCH_INK)

axA.set_ylabel("PERMANOVA $R^2$ (%)", labelpad=2)
axA.set_ylim(-0.6, 13)
axA.yaxis.set_major_locator(MultipleLocator(4))
axA.set_title("A", loc="left", fontsize=9, fontweight="bold", pad=3)

# ══════════════════════════════════════ Panel B: delta_R2 (%) ══════════════
rdx = [POS[n] for n in NS]
rdy = [resp_dR2[n] * 100 for n in NS]
bdx = [POS[n] for n in NS if n in coh_dR2]
bdy = [coh_dR2[n] * 100 for n in NS if n in coh_dR2]

axB.axhline(0, color="black", linewidth=0.6, zorder=2)
axB.plot(rdx, rdy, color=INK, marker="o", markersize=3.2, markerfacecolor=INK,
         markeredgecolor=INK, zorder=4, clip_on=False)
axB.plot(bdx, bdy, color=BATCH_INK, marker="s", markersize=3.4,
         linestyle=(0, (3.5, 2)), markerfacecolor="white",
         markeredgecolor=BATCH_INK, markeredgewidth=0.9, zorder=4, clip_on=False)

resp_label_dy = [16, -16, 16, -16]  # alternate above/below to avoid collisions near y=0
for (x, y), dy in zip(zip(rdx, rdy), resp_label_dy):
    va = "bottom" if dy > 0 else "top"
    axB.annotate(f"{y:+.2f}", (x, y), textcoords="offset points", xytext=(0, dy),
                 ha="center", va=va, fontsize=6.3, color=INK)
for x, y in zip(bdx, bdy):
    axB.annotate(f"{y:+.1f}", (x, y), textcoords="offset points", xytext=(0, 5),
                 ha="center", va="bottom", fontsize=6.3, color=BATCH_INK)

axB.set_ylabel(r"$\Delta R^2 = R^2 - E_0$ (%)", labelpad=2)
axB.set_ylim(-3.5, 11.5)
axB.yaxis.set_major_locator(MultipleLocator(3))
axB.set_title("B", loc="left", fontsize=9, fontweight="bold", pad=3)

for ax in (axA, axB):
    ax.set_xlim(-0.35, 3.35)
    ax.set_xticks(range(4))
    ax.set_xticklabels(XT)
    ax.tick_params(axis="x", pad=1.5)
    ax.set_xlabel("cumulative $n$", labelpad=1)
    ax.yaxis.grid(True, color=GRID, linewidth=0.4, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=2)

import os
os.makedirs("results/revision", exist_ok=True)
fig.savefig("results/revision/dose_response.pdf")
fig.savefig("results/revision/dose_response.png", dpi=600)
plt.close(fig)

print("Saved results/revision/dose_response.{pdf,png}")
print(df[df["level"] == "genus"].to_string(index=False))
