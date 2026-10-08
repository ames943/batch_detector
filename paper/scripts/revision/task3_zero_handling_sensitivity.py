#!/usr/bin/env python3
"""
Task 3 — Zero-handling sensitivity. Redo genus n=39/79/118/283 PERMANOVA under
three zero-replacement strategies: pseudocount eps=0.5, eps=1, multiplicative
replacement (skbio). Checks whether cohort delta_R2 stays large and response
delta_R2 stays ~0 regardless of how zeros are handled.
"""
import sys, time
sys.path.insert(0, "scripts/revision")
import numpy as np
import pandas as pd
from lib import (clr_pseudocount, clr_multiplicative_replacement, aitchison_dist,
                  permanova, E0, cohort_of_n283)

OUT = "results/revision"
N_PERMS = 999
SEED = 42
rows = []

raw118 = pd.read_csv("results/ml/n118_3cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
labels118 = pd.read_csv("metadata/response_labels_3cohort.tsv", sep="\t").set_index("run_accession").reindex(raw118.index)
raw283 = pd.read_csv("results/ml/n283_4cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
labels283 = pd.read_csv("results/ml/n283_4cohort/response_labels_n283.tsv", sep="\t").set_index("run_accession").reindex(raw283.index)
cohort283 = pd.Series([cohort_of_n283(s) for s in raw283.index], index=raw283.index)

METHODS = {
    "pseudocount_0.5": lambda X: clr_pseudocount(X, 0.5),
    "pseudocount_1.0": lambda X: clr_pseudocount(X, 1.0),
    "multiplicative_replacement": lambda X: clr_multiplicative_replacement(X, delta=1e-5),
}

datasets = []
for name, prefixes in [("n39", ["SRR5930"]), ("n79", ["SRR5930", "SRR11413"]),
                        ("n118", ["SRR5930", "SRR11413", "SRR6000"])]:
    mask = raw118.index.to_series().apply(lambda s: any(s.startswith(p) for p in prefixes))
    datasets.append((f"genus_{name}", raw118.loc[mask], labels118.loc[mask, "response"],
                      labels118.loc[mask, "cohort"] if len(prefixes) > 1 else None))
datasets.append(("genus_n283", raw283, labels283["response"], cohort283))

t0 = time.time()
for dname, raw_sub, resp_sub, coh_sub in datasets:
    n = len(raw_sub)
    for mname, fn in METHODS.items():
        clr_vals = fn(raw_sub.values.astype(float))
        D = aitchison_dist(clr_vals)
        res_r = permanova(D, resp_sub.values, n_perms=N_PERMS, seed=SEED)
        e0_r = E0(2, n)
        rows.append(dict(dataset=dname, zero_method=mname, n=n, factor="response", groups=2,
                          R2=round(res_r["R2"], 6), E0=round(e0_r, 6),
                          delta_R2=round(res_r["R2"] - e0_r, 6), p_value=res_r["p_value"]))
        if coh_sub is not None:
            g = coh_sub.nunique()
            res_c = permanova(D, coh_sub.values, n_perms=N_PERMS, seed=SEED + 1)
            e0_c = E0(g, n)
            rows.append(dict(dataset=dname, zero_method=mname, n=n, factor="cohort", groups=g,
                              R2=round(res_c["R2"], 6), E0=round(e0_c, 6),
                              delta_R2=round(res_c["R2"] - e0_c, 6), p_value=res_c["p_value"]))
        print(f"[{time.strftime('%H:%M:%S')}] {dname} / {mname} done ({time.time()-t0:.1f}s elapsed)", flush=True)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/task3_zero_handling_sensitivity.tsv", sep="\t", index=False)
print("\n" + df.to_string(index=False))

# Summary check
print("\n--- Summary: response delta_R2 range and cohort delta_R2 range per zero-handling method ---")
summ = df.groupby(["dataset", "factor"])["delta_R2"].agg(["min", "max"]).reset_index()
print(summ.to_string(index=False))
summ.to_csv(f"{OUT}/task3_summary_range.tsv", sep="\t", index=False)
print(f"\n[{time.strftime('%H:%M:%S')}] Task 3 complete.", flush=True)
