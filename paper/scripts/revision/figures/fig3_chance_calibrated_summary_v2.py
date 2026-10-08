#!/usr/bin/env python3
"""
Fig 3 — Chance-calibrated Delta_R2 summary across datasets, single column.
Source: results/revision/task1_chance_calibrated_permanova.tsv,
        results/revision/task4_c4_alone_permanova.tsv (no recompute).
"""
import sys
sys.path.insert(0, "scripts/revision/figures")
import pandas as pd
import numpy as np
from style import apply, save, panel_letter, RESPONSE, COHORT

apply()
import matplotlib.pyplot as plt

task1 = pd.read_csv("results/revision/task1_chance_calibrated_permanova.tsv", sep="\t")
task4 = pd.read_csv("results/revision/task4_c4_alone_permanova.tsv", sep="\t")


def get1(dataset, factor):
    sub = task1[(task1["dataset"] == dataset) & (task1["factor"] == factor)].iloc[0]
    return float(sub["delta_R2"]), float(sub["p_value"])


def get4(factor):
    sub = task4[task4["factor"] == factor].iloc[0]
    return float(sub["delta_R2"]), float(sub["p_value"])


groups = [
    ("Microbiome\nn=79",  *get1("genus_n79", "response"),  *get1("genus_n79", "cohort")),
    ("Microbiome\nn=118", *get1("genus_n118", "response"), *get1("genus_n118", "cohort")),
    ("Microbiome\nn=283", *get1("genus_n283", "response"), *get1("genus_n283", "cohort")),
    ("C4 alone\nn=165\n(site)", *get4("response"), *get4("site")),
    ("Tumor\nn=223",      *get1("tumor_n223", "response"), *get1("tumor_n223", "cohort")),
    ("CRC (Duvallet)\nn=570", *get1("duvallet_n570", "response"), *get1("duvallet_n570", "cohort")),
]
# each tuple: (label, resp_dR2, resp_p, coh_dR2, coh_p)

labels = [g[0] for g in groups]
resp_d = np.array([g[1] * 100 for g in groups])
resp_p = np.array([g[2] for g in groups])
coh_d = np.array([g[3] * 100 for g in groups])
coh_p = np.array([g[4] for g in groups])

y = np.arange(len(groups))
h = 0.36

fig, ax = plt.subplots(figsize=(3.5, 3.6))
fig.subplots_adjust(left=0.30, right=0.95, top=0.94, bottom=0.10)

ax.axvline(0, color="black", linewidth=0.7, zorder=2)
bars_r = ax.barh(y + h / 2, resp_d, height=h, color=RESPONSE, label="Response", zorder=3)
bars_c = ax.barh(y - h / 2, coh_d, height=h, color=COHORT,
                  label=("Cohort" if True else "Site"), zorder=3)

xmax = max(coh_d.max(), resp_d.max()) * 1.35
xmin = min(0, resp_d.min()) * 1.6
ax.set_xlim(xmin, xmax)

for yi, (d, p) in zip(y, zip(resp_d, resp_p)):
    star = "*" if p < 0.05 else ""
    # Response bars are all tiny (near-zero delta_R2); always anchor the label
    # just right of x=0 regardless of sign, so it never lands under the
    # y-axis category labels on the left.
    ax.annotate(f"p={p:.3f}{star}", (max(d, 0), yi + h / 2), textcoords="offset points",
                xytext=(5, 0), ha="left", va="center", fontsize=7, color=RESPONSE)

for yi, (d, p) in zip(y, zip(coh_d, coh_p)):
    star = "*" if p < 0.05 else ""
    ha = "left" if d >= 0 else "right"
    ax.annotate(f"{d:+.1f}%{star}", (d, yi - h / 2), textcoords="offset points",
                xytext=(3 if d >= 0 else -3, 0), ha=ha, va="center", fontsize=7, color=COHORT)

ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=7)
ax.set_xlabel(r"$\Delta R^2 = R^2 - E_0$ (%)")
ax.invert_yaxis()
ax.spines[["top", "right"]].set_visible(False)
ax.xaxis.grid(True, color="#dddddd", linewidth=0.4, zorder=0)
ax.set_axisbelow(True)
ax.legend(handles=[bars_r, bars_c], labels=["Response", "Cohort / Site"],
          loc="upper right", fontsize=7, frameon=False)

save(fig, "results/revision/figures/fig3_chance_calibrated_summary_v2")

print("Values (dataset, response dR2%, response p, cohort/site dR2%, cohort/site p):")
for g in groups:
    print(f"  {g[0].replace(chr(10),' ')}: resp_dR2={g[1]*100:+.3f}% p={g[2]:.3f}  "
          f"coh_dR2={g[3]*100:+.3f}% p={g[4]:.3f}")
