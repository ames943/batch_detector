#!/usr/bin/env python3
"""
Fig 2 — Aitchison-distance PCoA ordination.
Source data (no re-analysis, only a deterministic PCoA projection of the
existing CLR matrices for visualization):
  results/ml/n283_4cohort/X_genus_clr.tsv (+ response_labels_n283.tsv)
  results/ml/lee2022/X_genus_clr.tsv (+ metadata/lee2022_labels.tsv, lee2022_sites.tsv)
Annotated R2/E0/delta_R2/p values are pulled verbatim from
results/revision/task1_chance_calibrated_permanova.tsv and task4_c4_alone_permanova.tsv
-- NOT recomputed.
"""
import sys
sys.path.insert(0, "scripts/revision/figures")
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from style import apply, save, panel_letter, OKABE_ITO, RESPONSE, COHORT

apply()
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


def pcoa(D):
    n = D.shape[0]
    A = -0.5 * (D ** 2)
    J = np.eye(n) - np.ones((n, n)) / n
    G = J @ A @ J
    eigvals, eigvecs = np.linalg.eigh(G)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals > 1e-10
    total = eigvals[pos].sum()
    pct = eigvals / total * 100
    coords = eigvecs * np.sqrt(np.clip(eigvals, 0, None))
    return coords, pct


def lookup(df, dataset, factor, col):
    sub = df[(df["dataset"] == dataset) & (df["factor"] == factor)]
    return float(sub.iloc[0][col])


task1 = pd.read_csv("results/revision/task1_chance_calibrated_permanova.tsv", sep="\t")
task4 = pd.read_csv("results/revision/task4_c4_alone_permanova.tsv", sep="\t")

# ── n=283 data ────────────────────────────────────────────────────────────────
clr283 = pd.read_csv("results/ml/n283_4cohort/X_genus_clr.tsv", sep="\t", index_col="run_accession")
lab283 = pd.read_csv("results/ml/n283_4cohort/response_labels_n283.tsv", sep="\t").set_index("run_accession").reindex(clr283.index)
D283 = squareform(pdist(clr283.values, metric="euclidean"))
coords283, pct283 = pcoa(D283)

# ── C4 alone data ────────────────────────────────────────────────────────────
clr_c4 = pd.read_csv("results/ml/lee2022/X_genus_clr.tsv", sep="\t", index_col="run_accession")
resp_c4 = pd.read_csv("metadata/lee2022_labels.tsv", sep="\t").set_index("run_accession").reindex(clr_c4.index)["response"]
site_c4 = pd.read_csv("metadata/lee2022_sites.tsv", sep="\t").set_index("run_accession").reindex(clr_c4.index)["site"]
keep = resp_c4.notna() & site_c4.notna()
clr_c4, resp_c4, site_c4 = clr_c4.loc[keep], resp_c4.loc[keep], site_c4.loc[keep]
D_c4 = squareform(pdist(clr_c4.values, metric="euclidean"))
coords_c4, pct_c4 = pcoa(D_c4)

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(7.1, 2.7))
fig.subplots_adjust(left=0.07, right=0.985, top=0.87, bottom=0.30, wspace=0.32)


def scatter_panel(ax, coords, pct, color_labels, color_map, shape_labels, title_letter,
                   legend_title, textbox, marker="o"):
    x, y = coords[:, 0], coords[:, 1]
    for cat in pd.unique(color_labels):
        m_cat = color_labels == cat
        # response(-like) shape: 'R' filled, else open
        m_r = m_cat & (shape_labels == shape_labels.unique()[0]) if hasattr(shape_labels, "unique") else m_cat
        ax.scatter(x[m_cat & (shape_labels == "R")], y[m_cat & (shape_labels == "R")],
                   s=14, marker="o", facecolor=color_map[cat], edgecolor=color_map[cat],
                   linewidth=0.5, zorder=3, label=None)
        ax.scatter(x[m_cat & (shape_labels == "NR")], y[m_cat & (shape_labels == "NR")],
                   s=16, marker="o", facecolor="none", edgecolor=color_map[cat],
                   linewidth=0.9, zorder=3, label=None)
    ax.set_xlabel(f"PCo1 ({pct[0]:.1f}%)", labelpad=2)
    ax.set_ylabel(f"PCo2 ({pct[1]:.1f}%)", labelpad=2)
    panel_letter(ax, title_letter)
    ax.text(0.02, 0.02, textbox, transform=ax.transAxes, fontsize=7,
            va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#bbbbbb", alpha=0.9))
    ax.tick_params(length=2)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


# Panel A: n=283, color=cohort, shape=response
cohort_cats = sorted(lab283["cohort"].unique())
cohort_color_map = {c: OKABE_ITO[i] for i, c in enumerate(cohort_cats)}
r2_coh = lookup(task1, "genus_n283", "cohort", "R2")
e0_coh = lookup(task1, "genus_n283", "cohort", "E0")
d_coh = lookup(task1, "genus_n283", "cohort", "delta_R2")
p_coh = lookup(task1, "genus_n283", "cohort", "p_value")
txtA = f"Cohort: $R^2$={r2_coh*100:.2f}%\n$E_0$={e0_coh*100:.2f}%  $\\Delta R^2$={d_coh*100:+.2f}%\np={p_coh:.3f}"
scatter_panel(axA, coords283, pct283, lab283["cohort"].values, cohort_color_map,
              lab283["response"].values, "A", "Cohort", txtA)
handles_a = [mlines.Line2D([], [], marker="o", linestyle="", color=cohort_color_map[c],
                            markersize=4.5, label=c) for c in cohort_cats]
handles_a += [mlines.Line2D([], [], marker="o", linestyle="", markerfacecolor="black",
                             markeredgecolor="black", markersize=4.5, label="R"),
              mlines.Line2D([], [], marker="o", linestyle="", markerfacecolor="none",
                             markeredgecolor="black", markersize=4.5, label="NR")]
axA.legend(handles=handles_a, fontsize=7, loc="upper right", frameon=True,
           facecolor="white", edgecolor="none", framealpha=0.85,
           handletextpad=0.3, labelspacing=0.25, borderaxespad=0.1, ncol=1)

# Panel B: n=283, color=response only (same coordinates)
resp_color_map = {"R": RESPONSE, "NR": COHORT}
r2_resp = lookup(task1, "genus_n283", "response", "R2")
e0_resp = lookup(task1, "genus_n283", "response", "E0")
d_resp = lookup(task1, "genus_n283", "response", "delta_R2")
p_resp = lookup(task1, "genus_n283", "response", "p_value")
txtB = f"Response: $R^2$={r2_resp*100:.3f}%\n$E_0$={e0_resp*100:.3f}%  $\\Delta R^2$={d_resp*100:+.3f}%\np={p_resp:.3f}"
scatter_panel(axB, coords283, pct283, lab283["response"].values, resp_color_map,
              lab283["response"].values, "B", "Response", txtB)
handles_b = [mlines.Line2D([], [], marker="o", linestyle="", markerfacecolor=RESPONSE,
                            markeredgecolor=RESPONSE, markersize=4.5, label="R"),
             mlines.Line2D([], [], marker="o", linestyle="", markerfacecolor="none",
                            markeredgecolor=COHORT, markersize=4.5, label="NR")]
axB.legend(handles=handles_b, fontsize=7, loc="upper right", frameon=True,
           facecolor="white", edgecolor="none", framealpha=0.85,
           handletextpad=0.3, labelspacing=0.25, borderaxespad=0.1)

# Panel C: C4 alone, color=site, shape=response
site_cats = sorted(site_c4.unique())
site_color_map = {s: OKABE_ITO[i] for i, s in enumerate(site_cats)}
r2_site = float(task4[task4["factor"] == "site"]["R2"].iloc[0])
e0_site = float(task4[task4["factor"] == "site"]["E0"].iloc[0])
d_site = float(task4[task4["factor"] == "site"]["delta_R2"].iloc[0])
p_site = float(task4[task4["factor"] == "site"]["p_value"].iloc[0])
r2_resp_c4 = float(task4[task4["factor"] == "response"]["R2"].iloc[0])
d_resp_c4 = float(task4[task4["factor"] == "response"]["delta_R2"].iloc[0])
p_resp_c4 = float(task4[task4["factor"] == "response"]["p_value"].iloc[0])
txtC = (f"Site: $R^2$={r2_site*100:.2f}%, $\\Delta R^2$={d_site*100:+.2f}%\n"
        f"   p={p_site:.3f}\n"
        f"Response: $R^2$={r2_resp_c4*100:.2f}%, $\\Delta R^2$={d_resp_c4*100:+.2f}%\n"
        f"   p={p_resp_c4:.3f}")
scatter_panel(axC, coords_c4, pct_c4, site_c4.values, site_color_map,
              resp_c4.values, "C", "Site", txtC)
handles_c = [mlines.Line2D([], [], marker="o", linestyle="", color=site_color_map[s],
                            markersize=4.2, label=s) for s in site_cats]
axC.legend(handles=handles_c, fontsize=7, loc="upper right", frameon=True,
           facecolor="white", edgecolor="none", framealpha=0.85,
           handletextpad=0.3, labelspacing=0.22, borderaxespad=0.1)

save(fig, "results/revision/figures/fig2_ordination")

print(f"Panel A/B (n=283): PCo1 pct={pct283[0]:.2f}  PCo2 pct={pct283[1]:.2f}")
print(f"Panel C (C4, n={len(clr_c4)}): PCo1 pct={pct_c4[0]:.2f}  PCo2 pct={pct_c4[1]:.2f}")
print(f"Cohort R2={r2_coh:.5f} E0={e0_coh:.5f} dR2={d_coh:.5f} p={p_coh:.3f}")
print(f"Response(n283) R2={r2_resp:.5f} E0={e0_resp:.5f} dR2={d_resp:.5f} p={p_resp:.3f}")
print(f"Site R2={r2_site:.5f} E0={e0_site:.5f} dR2={d_site:.5f} p={p_site:.3f}")
print(f"Response(C4) R2={r2_resp_c4:.5f} dR2={d_resp_c4:.5f} p={p_resp_c4:.3f}")
