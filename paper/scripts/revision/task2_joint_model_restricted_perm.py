#!/usr/bin/env python3
"""
Task 2 — Joint (Type-III "by=margin") two-factor PERMANOVA with restricted
permutation, Python substitute for:
    adonis2(dist ~ cohort + response, by="margin",
            permutations=how(nperm=999, blocks=cohort))
(No R/vegan available in this environment — see scripts/revision/lib.py docstring.)

Run on: n=79, n=118, n=283 (genus), tumor n=223, Duvallet n=570.
"""
import sys, time
sys.path.insert(0, "scripts/revision")
import numpy as np
import pandas as pd
from lib import clr_transform, aitchison_dist, margin_permanova_two_factor, cohort_of_n283
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler

OUT = "results/revision"
N_PERMS = 999
SEED = 42
rows = []

t0 = time.time()

raw118 = pd.read_csv("results/ml/n118_3cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
labels118 = pd.read_csv("metadata/response_labels_3cohort.tsv", sep="\t").set_index("run_accession").reindex(raw118.index)

for name, prefixes in [("genus_n79", ["SRR5930", "SRR11413"]),
                        ("genus_n118", ["SRR5930", "SRR11413", "SRR6000"])]:
    mask = raw118.index.to_series().apply(lambda s: any(s.startswith(p) for p in prefixes))
    raw_sub = raw118.loc[mask]
    resp_sub = labels118.loc[mask, "response"].values
    coh_sub = labels118.loc[mask, "cohort"].values
    clr_sub = clr_transform(raw_sub.values)
    D = aitchison_dist(clr_sub)
    print(f"[{time.strftime('%H:%M:%S')}] {name} (n={len(raw_sub)}) joint model ...", flush=True)
    res = margin_permanova_two_factor(D, coh_sub, resp_sub, n_perms=N_PERMS, seed=SEED)
    res["dataset"] = name
    rows.append(res)
    print(f"    R2_resp|coh={res['R2_response_margin']:.5f} p={res['p_response_margin']:.4f}  "
          f"R2_coh|resp={res['R2_cohort_margin']:.5f} p={res['p_cohort_margin']:.4f}  "
          f"({time.time()-t0:.1f}s elapsed)", flush=True)

raw283 = pd.read_csv("results/ml/n283_4cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
labels283 = pd.read_csv("results/ml/n283_4cohort/response_labels_n283.tsv", sep="\t").set_index("run_accession").reindex(raw283.index)
cohort283 = np.array([cohort_of_n283(s) for s in raw283.index])
clr283 = clr_transform(raw283.values)
D283 = aitchison_dist(clr283)
print(f"[{time.strftime('%H:%M:%S')}] genus_n283 (n=283) joint model ...", flush=True)
res283 = margin_permanova_two_factor(D283, cohort283, labels283["response"].values, n_perms=N_PERMS, seed=SEED)
res283["dataset"] = "genus_n283"
rows.append(res283)
print(f"    R2_resp|coh={res283['R2_response_margin']:.5f} p={res283['p_response_margin']:.4f}  "
      f"R2_coh|resp={res283['R2_cohort_margin']:.5f} p={res283['p_cohort_margin']:.4f}  "
      f"({time.time()-t0:.1f}s elapsed)", flush=True)

tdf = pd.read_csv("results/ml/tumor/combined_tumor_features.tsv", sep="\t").dropna(subset=["response", "TMB"])
feat_cols = ["TMB", "n_mutations"] + [c for c in tdf.columns if c.startswith("mut_")]
Xt = StandardScaler().fit_transform(tdf[feat_cols].fillna(0).values.astype(float))
Dt = squareform(pdist(Xt, metric="euclidean"))
print(f"[{time.strftime('%H:%M:%S')}] tumor_n223 joint model ...", flush=True)
res_t = margin_permanova_two_factor(Dt, tdf["study"].values, tdf["response"].values, n_perms=N_PERMS, seed=SEED)
res_t["dataset"] = "tumor_n223"
rows.append(res_t)
print(f"    R2_resp|coh={res_t['R2_response_margin']:.5f} p={res_t['p_response_margin']:.4f}  "
      f"R2_coh|resp={res_t['R2_cohort_margin']:.5f} p={res_t['p_cohort_margin']:.4f}  "
      f"({time.time()-t0:.1f}s elapsed)", flush=True)

dfm = pd.read_csv("results/batch_detector_external_validation/input_feature_matrix.tsv", sep="\t", index_col="sample_id")
dlab = pd.read_csv("results/batch_detector_external_validation/input_labels.tsv", sep="\t", index_col="sample_id").reindex(dfm.index)
feat_var = dfm.var(axis=0)
dfm = dfm.loc[:, feat_var > 0]
clr_d = clr_transform(dfm.values, pseudocount=1e-6)
Dd = aitchison_dist(clr_d)
print(f"[{time.strftime('%H:%M:%S')}] duvallet_n570 joint model ...", flush=True)
res_d = margin_permanova_two_factor(Dd, dlab["study"].values, dlab["response"].values, n_perms=N_PERMS, seed=SEED)
res_d["dataset"] = "duvallet_n570"
rows.append(res_d)
print(f"    R2_resp|coh={res_d['R2_response_margin']:.5f} p={res_d['p_response_margin']:.4f}  "
      f"R2_coh|resp={res_d['R2_cohort_margin']:.5f} p={res_d['p_cohort_margin']:.4f}  "
      f"({time.time()-t0:.1f}s elapsed)", flush=True)

df = pd.DataFrame(rows)
cols = ["dataset", "n", "R2_response_margin", "p_response_margin",
        "R2_cohort_margin", "p_cohort_margin",
        "response_perm_restricted_within_cohort", "cohort_perm_restricted", "n_perms"]
df = df[cols]
df.to_csv(f"{OUT}/task2_joint_model_margin_permanova.tsv", sep="\t", index=False)
print("\n" + df.to_string(index=False))
print(f"\n[{time.strftime('%H:%M:%S')}] Task 2 complete. Total time: {time.time()-t0:.1f}s", flush=True)
