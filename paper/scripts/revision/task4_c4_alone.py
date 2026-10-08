#!/usr/bin/env python3
"""Task 4 — C4 (Lee 2022, n=165) alone: PERMANOVA of response and of
recruiting site (5 sites), with chance-calibrated E0 and delta_R2."""
import sys, time
sys.path.insert(0, "scripts/revision")
import pandas as pd
from lib import clr_transform, aitchison_dist, permanova, E0

OUT = "results/revision"
N_PERMS = 999
SEED = 42

raw = pd.read_csv("results/ml/lee2022/X_genus_raw.tsv", sep="\t", index_col="run_accession")
resp = pd.read_csv("metadata/lee2022_labels.tsv", sep="\t").set_index("run_accession").reindex(raw.index)["response"]
site = pd.read_csv("metadata/lee2022_sites.tsv", sep="\t").set_index("run_accession").reindex(raw.index)["site"]

keep = resp.notna() & site.notna()
raw = raw.loc[keep]; resp = resp.loc[keep]; site = site.loc[keep]
n = len(raw)
print(f"n={n}, sites={site.value_counts().to_dict()}, response={resp.value_counts().to_dict()}")

clr = clr_transform(raw.values)
D = aitchison_dist(clr)

res_r = permanova(D, resp.values, n_perms=N_PERMS, seed=SEED)
e0_r = E0(2, n)
res_s = permanova(D, site.values, n_perms=N_PERMS, seed=SEED + 1)
e0_s = E0(site.nunique(), n)

rows = [
    dict(dataset="C4_alone", n=n, factor="response", groups=2,
         R2=round(res_r["R2"], 6), E0=round(e0_r, 6), delta_R2=round(res_r["R2"] - e0_r, 6),
         p_value=res_r["p_value"]),
    dict(dataset="C4_alone", n=n, factor="site", groups=site.nunique(),
         R2=round(res_s["R2"], 6), E0=round(e0_s, 6), delta_R2=round(res_s["R2"] - e0_s, 6),
         p_value=res_s["p_value"]),
]
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/task4_c4_alone_permanova.tsv", sep="\t", index=False)
print(df.to_string(index=False))
print(f"\n[{time.strftime('%H:%M:%S')}] Task 4 complete.")
