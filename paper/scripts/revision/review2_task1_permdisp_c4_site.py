#!/usr/bin/env python3
"""
Review2 Task 1 — PERMDISP for recruiting site in C4 (n=165, 5 sites).
Same settings as Task 4: genus level, Aitchison distance, CLR pseudocount=1e-6,
999 permutations. PERMDISP implementation copied verbatim (Anderson 2006,
PCoA + distance-to-group-centroid + permutation) from
scripts/steps_1_4_c4_integration.py, for methodological consistency with the
rest of the project.
"""
import sys
sys.path.insert(0, "scripts/revision")
import numpy as np
import pandas as pd
from lib import clr_transform, aitchison_dist

OUT = "results/revision/review2"
N_PERMS = 999
SEED = 42


def permdisp(D, grouping, n_perms=999, seed=42):
    grp = np.asarray(grouping)
    n = D.shape[0]
    groups = np.unique(grp)

    D2 = D ** 2
    A = -0.5 * D2
    ones = np.ones((n, 1))
    H = np.eye(n) - ones @ ones.T / n
    B = H @ A @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos_mask = eigvals > 1e-10
    coords = eigvecs[:, pos_mask] * np.sqrt(eigvals[pos_mask])

    def dist_to_centroid(coords, grp):
        d = np.zeros(n)
        for g in groups:
            mask = grp == g
            centroid = coords[mask].mean(axis=0)
            diffs = coords[mask] - centroid
            d[mask] = np.sqrt((diffs ** 2).sum(axis=1))
        return d

    d_obs = dist_to_centroid(coords, grp)
    group_disp = {str(g): round(float(d_obs[grp == g].mean()), 4) for g in groups}

    grand_mean = d_obs.mean()
    q = len(groups)
    SS_A = sum((grp == g).sum() * (d_obs[grp == g].mean() - grand_mean) ** 2 for g in groups)
    SS_W = sum(((d_obs[grp == g] - d_obs[grp == g].mean()) ** 2).sum() for g in groups)
    F_obs = (SS_A / (q - 1)) / (SS_W / (n - q))

    rng = np.random.default_rng(seed)
    perm_F = np.empty(n_perms)
    for i in range(n_perms):
        g_perm = rng.permutation(grp)
        d_p = dist_to_centroid(coords, g_perm)
        gm_p = d_p.mean()
        ssa_p = sum((g_perm == g).sum() * (d_p[g_perm == g].mean() - gm_p) ** 2 for g in groups)
        ssw_p = sum(((d_p[g_perm == g] - d_p[g_perm == g].mean()) ** 2).sum() for g in groups)
        perm_F[i] = (ssa_p / (q - 1)) / (ssw_p / (n - q))

    p_val = float((perm_F >= F_obs).sum() + 1) / (n_perms + 1)
    return dict(n=n, n_groups=q, F_stat=round(float(F_obs), 4), p_value=round(p_val, 4),
                n_perms=n_perms, group_dispersions=group_disp,
                dispersion_homogeneous=p_val > 0.05)


raw = pd.read_csv("results/ml/lee2022/X_genus_raw.tsv", sep="\t", index_col="run_accession")
resp = pd.read_csv("metadata/lee2022_labels.tsv", sep="\t").set_index("run_accession").reindex(raw.index)["response"]
site = pd.read_csv("metadata/lee2022_sites.tsv", sep="\t").set_index("run_accession").reindex(raw.index)["site"]
keep = resp.notna() & site.notna()
raw, site = raw.loc[keep], site.loc[keep]
n = len(raw)

clr = clr_transform(raw.values)
D = aitchison_dist(clr)

res = permdisp(D, site.values, n_perms=N_PERMS, seed=SEED)
print(f"n={n}, sites={site.value_counts().to_dict()}")
print(res)

out = pd.DataFrame([dict(dataset="C4_alone", n=n, factor="site", n_groups=res["n_groups"],
                          F_stat=res["F_stat"], p_value=res["p_value"], n_perms=N_PERMS,
                          dispersion_homogeneous=res["dispersion_homogeneous"],
                          group_dispersions=str(res["group_dispersions"]))])
out.to_csv(f"{OUT}/permdisp_c4_site.tsv", sep="\t", index=False)
print(f"\nSaved: {OUT}/permdisp_c4_site.tsv")
