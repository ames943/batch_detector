#!/usr/bin/env python3
"""
Task 1 — Chance-calibrated PERMANOVA across every dataset in the paper.

For each dataset: R2, E0=(g-1)/(n-1), delta_R2=R2-E0, permutation p, for both
response and cohort/batch (where applicable). Also: (a) empirical sanity check
that mean permuted R2 ~= E0 for one dataset, (b) permutation null for
rho = R2_cohort / R2_response (999 reps, independent shuffles of each factor).

Run from cancer_project/:  python3 scripts/revision/task1_chance_calibrated_permanova.py
"""
import sys, time
sys.path.insert(0, "scripts/revision")
import numpy as np
import pandas as pd
from lib import (clr_transform, aitchison_dist, permanova, E0,
                  rho_permutation_null, cohort_of_n283)

OUT = "results/revision"
N_PERMS = 999
SEED = 42

rows = []          # main R2/E0/dR2/p table
rho_rows = []       # rho permutation null table
sanity_rows = []    # empirical mean-permuted-R2 vs E0 sanity check


def add_row(dataset, level, n, factor, g, r2, p, note=""):
    e0 = E0(g, n)
    rows.append(dict(dataset=dataset, level=level, n=n, factor=factor, groups=g,
                      R2=round(r2, 6), E0=round(e0, 6), delta_R2=round(r2 - e0, 6),
                      p_value=p, note=note))


def permuted_r2_distribution(D, grp, n_perms, seed):
    """Full permuted-R2 distribution (not just F) for the sanity check."""
    from lib import _ss_total, _ss_within
    d2 = D ** 2
    ST = _ss_total(d2)
    rng = np.random.default_rng(seed)
    r2s = np.empty(n_perms)
    for i in range(n_perms):
        g_p = rng.permutation(grp)
        sw = _ss_within(d2, g_p)
        r2s[i] = (ST - sw) / ST if ST > 0 else 0.0
    return r2s


def run_microbiome(name, clr_df, response, cohort=None, level="genus"):
    n = len(clr_df)
    D = aitchison_dist(clr_df.values)

    res = permanova(D, response.values, n_perms=N_PERMS, seed=SEED)
    add_row(name, level, n, "response", 2, res["R2"], res["p_value"])

    if cohort is not None:
        g = cohort.nunique()
        resb = permanova(D, cohort.values, n_perms=N_PERMS, seed=SEED + 1)
        add_row(name, level, n, "cohort", g, resb["R2"], resb["p_value"])

        rr = rho_permutation_null(D, cohort.values, response.values, n_perms=N_PERMS, seed=SEED + 2)
        rho_rows.append(dict(dataset=name, level=level, n=n, **rr))

    return D


print(f"[{time.strftime('%H:%M:%S')}] Loading matrices ...", flush=True)

raw118 = pd.read_csv("results/ml/n118_3cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
labels118 = pd.read_csv("metadata/response_labels_3cohort.tsv", sep="\t").set_index("run_accession")
labels118 = labels118.reindex(raw118.index)

raw283 = pd.read_csv("results/ml/n283_4cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
labels283 = pd.read_csv("results/ml/n283_4cohort/response_labels_n283.tsv", sep="\t").set_index("run_accession")
labels283 = labels283.reindex(raw283.index)

# ── n=39, 79, 118 genus (subsets of the 118 matrix) ──────────────────────────
for name, prefixes in [("genus_n39", ["SRR5930"]),
                        ("genus_n79", ["SRR5930", "SRR11413"]),
                        ("genus_n118", ["SRR5930", "SRR11413", "SRR6000"])]:
    mask = raw118.index.to_series().apply(lambda s: any(s.startswith(p) for p in prefixes))
    raw_sub = raw118.loc[mask]
    resp_sub = labels118.loc[mask, "response"]
    coh_sub = labels118.loc[mask, "cohort"] if len(prefixes) > 1 else None
    clr_sub = pd.DataFrame(clr_transform(raw_sub.values), index=raw_sub.index, columns=raw_sub.columns)
    print(f"[{time.strftime('%H:%M:%S')}] {name}: n={len(raw_sub)}", flush=True)
    D = run_microbiome(name, clr_sub, resp_sub, coh_sub, level="genus")
    if name == "genus_n118":
        D118 = D  # keep for sanity check below

# ── n=283 genus ───────────────────────────────────────────────────────────────
cohort283 = pd.Series([cohort_of_n283(s) for s in raw283.index], index=raw283.index)
clr283 = pd.DataFrame(clr_transform(raw283.values), index=raw283.index, columns=raw283.columns)
print(f"[{time.strftime('%H:%M:%S')}] genus_n283: n={len(raw283)}", flush=True)
run_microbiome("genus_n283", clr283, labels283["response"], cohort283, level="genus")

# ── phylum / species at n=118 (already-CLR matrices from Phase 3c) ──────────
for level, path in [("phylum", "results/ml/phase3c/X_phylum_clr.tsv"),
                     ("species", "results/ml/phase3c/X_species_clr.tsv")]:
    clr_lvl = pd.read_csv(path, sep="\t", index_col="run_accession")
    clr_lvl = clr_lvl.reindex(raw118.index)  # same 118 samples, same order
    resp_lvl = labels118["response"]
    coh_lvl = labels118["cohort"]
    print(f"[{time.strftime('%H:%M:%S')}] {level}_n118: n={len(clr_lvl)}, features={clr_lvl.shape[1]}", flush=True)
    run_microbiome(f"{level}_n118", clr_lvl, resp_lvl, coh_lvl, level=level)

# ── Tumor n=223 (Euclidean on standardized features, per original tumor_batch.py) ──
print(f"[{time.strftime('%H:%M:%S')}] tumor_n223 ...", flush=True)
from sklearn.preprocessing import StandardScaler
tdf = pd.read_csv("results/ml/tumor/combined_tumor_features.tsv", sep="\t")
tdf = tdf.dropna(subset=["response", "TMB"])
feat_cols = ["TMB", "n_mutations"] + [c for c in tdf.columns if c.startswith("mut_")]
Xt = tdf[feat_cols].fillna(0).values.astype(float)
Xt_scaled = StandardScaler().fit_transform(Xt)
Dt = squareform_helper = None
from scipy.spatial.distance import pdist as _pdist, squareform as _squareform
Dt = _squareform(_pdist(Xt_scaled, metric="euclidean"))
y_t = tdf["response"].values
study_t = tdf["study"].values
res_t_resp = permanova(Dt, y_t, n_perms=N_PERMS, seed=SEED)
add_row("tumor_n223", "gene_flags", len(tdf), "response", 2, res_t_resp["R2"], res_t_resp["p_value"])
res_t_coh = permanova(Dt, study_t, n_perms=N_PERMS, seed=SEED + 1)
add_row("tumor_n223", "gene_flags", len(tdf), "cohort", pd.Series(study_t).nunique(), res_t_coh["R2"], res_t_coh["p_value"])
rr_t = rho_permutation_null(Dt, study_t, y_t, n_perms=N_PERMS, seed=SEED + 2)
rho_rows.append(dict(dataset="tumor_n223", level="gene_flags", n=len(tdf), **rr_t))

# ── Duvallet n=570 (raw genus counts -> CLR, per external validation run) ──────
print(f"[{time.strftime('%H:%M:%S')}] duvallet_n570 ...", flush=True)
dfm = pd.read_csv("results/batch_detector_external_validation/input_feature_matrix.tsv", sep="\t", index_col="sample_id")
dlab = pd.read_csv("results/batch_detector_external_validation/input_labels.tsv", sep="\t", index_col="sample_id")
dlab = dlab.reindex(dfm.index)
feat_var = dfm.var(axis=0)
dfm = dfm.loc[:, feat_var > 0]
clr_d = pd.DataFrame(clr_transform(dfm.values, pseudocount=1e-6), index=dfm.index, columns=dfm.columns)
Dd = aitchison_dist(clr_d.values)
res_d_resp = permanova(Dd, dlab["response"].values, n_perms=N_PERMS, seed=SEED)
add_row("duvallet_n570", "genus", len(dfm), "response", 2, res_d_resp["R2"], res_d_resp["p_value"])
res_d_coh = permanova(Dd, dlab["study"].values, n_perms=N_PERMS, seed=SEED + 1)
add_row("duvallet_n570", "genus", len(dfm), "cohort", dlab["study"].nunique(), res_d_coh["R2"], res_d_coh["p_value"])
rr_d = rho_permutation_null(Dd, dlab["study"].values, dlab["response"].values, n_perms=N_PERMS, seed=SEED + 2)
rho_rows.append(dict(dataset="duvallet_n570", level="genus", n=len(dfm), **rr_d))

# ── Sanity check: mean permuted R2 ~ E0, for genus_n118 response and cohort ──
print(f"[{time.strftime('%H:%M:%S')}] Sanity check: mean permuted R2 vs E0 (genus_n118) ...", flush=True)
resp118 = labels118["response"].values
coh118 = labels118["cohort"].values
r2s_resp = permuted_r2_distribution(D118, resp118, n_perms=999, seed=123)
r2s_coh = permuted_r2_distribution(D118, coh118, n_perms=999, seed=124)
sanity_rows.append(dict(dataset="genus_n118", factor="response", g=2, n=118,
                         mean_permuted_R2=float(r2s_resp.mean()), E0=E0(2, 118),
                         std_permuted_R2=float(r2s_resp.std())))
sanity_rows.append(dict(dataset="genus_n118", factor="cohort", g=3, n=118,
                         mean_permuted_R2=float(r2s_coh.mean()), E0=E0(3, 118),
                         std_permuted_R2=float(r2s_coh.std())))

# ── Save ──────────────────────────────────────────────────────────────────────
df_main = pd.DataFrame(rows)
df_main.to_csv(f"{OUT}/task1_chance_calibrated_permanova.tsv", sep="\t", index=False)
print(df_main.to_string(index=False))

df_rho = pd.DataFrame(rho_rows)
df_rho.to_csv(f"{OUT}/task1_rho_permutation_null.tsv", sep="\t", index=False)
print("\n", df_rho.to_string(index=False))

df_sanity = pd.DataFrame(sanity_rows)
df_sanity.to_csv(f"{OUT}/task1_sanity_check_E0.tsv", sep="\t", index=False)
print("\n", df_sanity.to_string(index=False))

print(f"\n[{time.strftime('%H:%M:%S')}] Task 1 complete.", flush=True)
