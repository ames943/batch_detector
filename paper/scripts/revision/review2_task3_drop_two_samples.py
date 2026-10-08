#!/usr/bin/env python3
"""
Review2 Task 3 — pooled PERMANOVA after dropping the two samples that are not
in the response-label curation (SRR6000943 in C3, ERR10290768 in C4), n = 281.
Same settings as the main n = 283 analysis: genus level, Aitchison distance,
CLR with pseudocount 1e-6, 999 permutations.

This reproduces results/revision/review2/permanova_n281_dropped2.tsv. The
original run of this step was done interactively; this script was written
afterwards from the same settings.
"""
import sys
sys.path.insert(0, "scripts/revision")
import pandas as pd
from lib import clr_transform, aitchison_dist, permanova, E0, cohort_of_n283

OUT = "results/revision/review2"
DROP = ["SRR6000943", "ERR10290768"]
N_PERMS = 999
SEED = 42

raw = pd.read_csv("results/ml/n283_4cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
labels = pd.read_csv("results/ml/n283_4cohort/response_labels_n283.tsv", sep="\t").set_index("run_accession")

raw = raw.drop(index=DROP)
resp = labels.reindex(raw.index)["response"]
cohort = pd.Series([cohort_of_n283(s) for s in raw.index], index=raw.index)
n = len(raw)

D = aitchison_dist(clr_transform(raw.values, pseudocount=1e-6))
res_r = permanova(D, resp.values, n_perms=N_PERMS, seed=SEED)
res_c = permanova(D, cohort.values, n_perms=N_PERMS, seed=SEED + 1)

rows = []
for factor, res, g in [("response", res_r, 2), ("cohort", res_c, cohort.nunique())]:
    e0 = E0(g, n)
    rows.append(dict(dataset="genus_n281_2sample_dropped", factor=factor, groups=g,
                     R2=round(res["R2"], 6), E0=round(e0, 6),
                     delta_R2=round(res["R2"] - e0, 6), p_value=res["p_value"]))

out = pd.DataFrame(rows)
out.to_csv(f"{OUT}/permanova_n281_dropped2_rerun.tsv", sep="\t", index=False)
print(out.to_string(index=False))
