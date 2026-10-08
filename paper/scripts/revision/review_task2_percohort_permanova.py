#!/usr/bin/env python3
"""Review Task 2 — response PERMANOVA run separately within C2 and C3 (C1, C4
already done in Task 4). Same settings: genus level, Aitchison distance,
CLR pseudocount=1e-6, 999 permutations."""
import sys
sys.path.insert(0, "scripts/revision")
import pandas as pd
from lib import clr_transform, aitchison_dist, permanova, E0

OUT = "results/revision/review"
N_PERMS = 999
SEED = 42

raw118 = pd.read_csv("results/ml/n118_3cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
labels118 = pd.read_csv("metadata/response_labels_3cohort.tsv", sep="\t").set_index("run_accession").reindex(raw118.index)

rows = []
for name, prefix in [("C2", "SRR11413"), ("C3", "SRR6000")]:
    mask = raw118.index.str.startswith(prefix)
    raw_sub = raw118.loc[mask]
    resp_sub = labels118.loc[mask, "response"]
    n = len(raw_sub)
    clr_sub = clr_transform(raw_sub.values)
    D = aitchison_dist(clr_sub)
    res = permanova(D, resp_sub.values, n_perms=N_PERMS, seed=SEED)
    e0 = E0(2, n)
    rows.append(dict(dataset=name, n=n, factor="response", groups=2,
                      R2=round(res["R2"], 6), E0=round(e0, 6),
                      delta_R2=round(res["R2"] - e0, 6), p_value=res["p_value"]))
    print(f"{name}: n={n}  R={int((resp_sub=='R').sum())} NR={int((resp_sub=='NR').sum())}  "
          f"R2={res['R2']:.5f}  E0={e0:.5f}  dR2={res['R2']-e0:+.5f}  p={res['p_value']:.3f}")

# Pull C1 and C4 from already-computed files (Task 1 / Task 4) for the combined table
task1 = pd.read_csv("results/revision/task1_chance_calibrated_permanova.tsv", sep="\t")
c1 = task1[(task1["dataset"] == "genus_n39") & (task1["factor"] == "response")].iloc[0]
task4 = pd.read_csv("results/revision/task4_c4_alone_permanova.tsv", sep="\t")
c4 = task4[task4["factor"] == "response"].iloc[0]
rows.insert(0, dict(dataset="C1", n=int(c1["n"]), factor="response", groups=2,
                     R2=c1["R2"], E0=c1["E0"], delta_R2=c1["delta_R2"], p_value=c1["p_value"]))
rows.append(dict(dataset="C4", n=int(c4["n"]), factor="response", groups=2,
                  R2=c4["R2"], E0=c4["E0"], delta_R2=c4["delta_R2"], p_value=c4["p_value"]))

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/percohort_permanova_C1_C4.tsv", sep="\t", index=False)
print("\nAll 4 cohorts:")
print(df.to_string(index=False))
