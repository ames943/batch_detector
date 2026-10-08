#!/usr/bin/env python3
"""
Task 6 — Meta-analysis redo with proper effect sizes.

Quick-fact 0d established that the original Phase 4 meta-analysis computed
per-cohort "SE" as std(LOOCV fold coefficients)/sqrt(n) -- LOOCV folds share
n-1 of n training samples, so they are highly correlated, and dividing by
sqrt(n) on top of that produces SEs far too small (a false-precision problem,
not a sampling-uncertainty estimate). This script redoes the meta-analysis
with a standard, analytically-correct effect size: per-cohort Hedges' g
(bias-corrected standardized mean difference, R vs NR) on CLR abundance,
with its closed-form SE, combined across cohorts via DerSimonian-Laird.
"""
import sys, time
sys.path.insert(0, "scripts/revision")
import numpy as np
import pandas as pd
from lib import hedges_g, dsl_meta, bh_fdr

OUT = "results/revision"
PREVALENCE_FRAC = 0.10

# ── Load per-cohort CLR (row-wise op -> identical whether computed on subset
#    or full matrix) and RAW (for prevalence) ─────────────────────────────────
clr118 = pd.read_csv("results/ml/n118_3cohort/X_genus_clr.tsv", sep="\t", index_col="run_accession")
raw118 = pd.read_csv("results/ml/n118_3cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
labels118 = pd.read_csv("metadata/response_labels_3cohort.tsv", sep="\t").set_index("run_accession").reindex(clr118.index)

clr_lee = pd.read_csv("results/ml/lee2022/X_genus_clr.tsv", sep="\t", index_col="run_accession")
raw_lee = pd.read_csv("results/ml/lee2022/X_genus_raw.tsv", sep="\t", index_col="run_accession")
resp_lee = pd.read_csv("metadata/lee2022_labels.tsv", sep="\t").set_index("run_accession").reindex(clr_lee.index)["response"]

cohorts = {}
for cname, prefix in [("cohort1", "SRR5930"), ("cohort2", "SRR11413"), ("cohort3", "SRR6000")]:
    mask = clr118.index.str.startswith(prefix)
    cohorts[cname] = dict(clr=clr118.loc[mask], raw=raw118.loc[mask], resp=labels118.loc[mask, "response"])
keep_lee = resp_lee.notna()
cohorts["cohort4"] = dict(clr=clr_lee.loc[keep_lee], raw=raw_lee.loc[keep_lee], resp=resp_lee.loc[keep_lee])

for c, d in cohorts.items():
    print(f"{c}: n={len(d['clr'])}, genera={d['clr'].shape[1]}, "
          f"R={int((d['resp']=='R').sum())} NR={int((d['resp']=='NR').sum())}")

# ── Per-cohort prevalence-filtered genus set + per-genus Hedges' g ──────────
per_cohort_effects = {}   # cname -> {genus: (g, se)}
for cname, d in cohorts.items():
    n = len(d["raw"])
    min_prev = max(2, int(np.ceil(PREVALENCE_FRAC * n)))
    prevalence = (d["raw"] > 0).sum(axis=0)
    prevalent_genera = prevalence[prevalence >= min_prev].index.tolist()
    y = d["resp"].values
    effects = {}
    for genus in prevalent_genera:
        vals = d["clr"][genus].values
        x_r = vals[y == "R"]
        x_nr = vals[y == "NR"]
        g, se = hedges_g(x_r, x_nr)
        if np.isfinite(g) and np.isfinite(se) and se > 0:
            effects[genus] = (g, se)
    per_cohort_effects[cname] = effects
    print(f"[{time.strftime('%H:%M:%S')}] {cname}: {len(prevalent_genera)} prevalent genera, "
          f"{len(effects)} with valid Hedges' g", flush=True)

# ── Save per-cohort effects for transparency ────────────────────────────────
all_genera = sorted(set().union(*[set(e.keys()) for e in per_cohort_effects.values()]))
pc_rows = []
for genus in all_genera:
    row = {"genus": genus}
    for cname in cohorts:
        ge = per_cohort_effects[cname].get(genus)
        row[f"hedges_g_{cname}"] = round(ge[0], 6) if ge else None
        row[f"se_{cname}"] = round(ge[1], 6) if ge else None
    pc_rows.append(row)
pd.DataFrame(pc_rows).to_csv(f"{OUT}/task6_per_cohort_hedges_g.tsv", sep="\t", index=False)

# ── DerSimonian-Laird meta-analysis across genera present in >=2 cohorts ────
meta_rows = []
for genus in all_genera:
    eff, var, coh_contrib = [], [], []
    for cname, effects in per_cohort_effects.items():
        if genus in effects:
            g, se = effects[genus]
            eff.append(g); var.append(se ** 2); coh_contrib.append(cname)
    if len(eff) < 2:
        continue
    eff = np.array(eff); var = np.array(var)
    res = dsl_meta(eff, var)
    dir_consistent = bool(np.all(eff > 0) or np.all(eff < 0))
    meta_rows.append(dict(
        genus=genus, n_cohorts=len(eff), cohorts=",".join(coh_contrib),
        pooled_hedges_g=round(res["theta_re"], 6), SE=round(res["se_re"], 6),
        ci_lo_95=round(res["ci_lo"], 6), ci_hi_95=round(res["ci_hi"], 6),
        tau2=round(res["tau2"], 6), I2=round(res["I2"], 1),
        Q=round(res["Q"], 4), Q_p=round(res["Q_p"], 4),
        z_stat=round(res["z"], 4), p_raw=round(res["p_z"], 6),
        direction_consistent=dir_consistent,
        direction="R+" if res["theta_re"] > 0 else "NR+",
    ))

meta_df = pd.DataFrame(meta_rows)
meta_df["q_value"] = bh_fdr(meta_df["p_raw"].values).round(6)
meta_df["fdr_sig"] = meta_df["q_value"] < 0.05
meta_df = meta_df.sort_values("p_raw").reset_index(drop=True)
meta_df.to_csv(f"{OUT}/task6_meta_analysis_hedges_g.tsv", sep="\t", index=False)

n_sig = int(meta_df["fdr_sig"].sum())
mean_I2 = round(float(meta_df["I2"].mean()), 1)
sig_genera = meta_df.loc[meta_df["fdr_sig"], "genus"].tolist()

print("\n" + "=" * 65)
print("TASK 6 — Hedges' g meta-analysis (proper SEs; supersedes Phase 4)")
print("=" * 65)
print(f"  Genera tested (>=2 cohorts): {len(meta_df)}")
print(f"  FDR q < 0.05               : {n_sig}")
print(f"  Mean I^2                   : {mean_I2}%")
print(f"  FDR-significant genera     : {sig_genera}")
print(meta_df.head(10).to_string(index=False))

with open(f"{OUT}/task6_summary.txt", "w") as f:
    f.write(f"Genera tested (>=2 cohorts): {len(meta_df)}\n")
    f.write(f"FDR q < 0.05: {n_sig}\n")
    f.write(f"Mean I2: {mean_I2}%\n")
    f.write(f"FDR-significant genera: {sig_genera}\n")

print(f"\n[{time.strftime('%H:%M:%S')}] Task 6 complete.")
