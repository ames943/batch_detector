#!/usr/bin/env python3
"""
Sensitivity analysis — V4 and V5 only.
V1-V3 completed in c4_sensitivity_analysis.py (results hardcoded below).
V4 rarefaction bug fixed: probs clipped and renormalized before multinomial draw.

Run from cancer_project/:
    python3 -u scripts/c4_sensitivity_v4v5.py 2>&1 | tee -a logs/c4_sensitivity.log
"""

import os, glob, time, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings("ignore", category=ConvergenceWarning)

CLR_PATH    = "results/ml/lee2022/X_genus_clr.tsv"
RAW_PATH    = "results/ml/lee2022/X_genus_raw.tsv"
LABELS_PATH = "metadata/lee2022_labels.tsv"
SITES_PATH  = "metadata/lee2022_sites.tsv"
REPORTS_DIR = "results/kraken_reports/lee2022"
OUT_DIR     = "results/ml/lee2022/sensitivity"
os.makedirs(OUT_DIR, exist_ok=True)

EN_C, EN_L1 = 1.0, 0.7
EN_KWARGS = dict(penalty="elasticnet", solver="saga", C=EN_C, l1_ratio=EN_L1,
                 class_weight="balanced", max_iter=5000, tol=1e-3, random_state=42)
DEFAULT_PREV_FRAC = 0.10
DEFAULT_VAR_FRAC  = 0.50
DEFAULT_TOP_N     = 100
SEED, N_PERMS     = 42, 200
RAREFY_TARGET     = 1_000_000
EXCLUDED_GENERA   = {"Homo"}

# ── Hardcoded V1–V3 results (from completed run) ──────────────────────────────
PRIOR_RESULTS = [
    {"variant": "V1 Baseline (prev=10%, top-100)",          "n": 165, "auc": 0.5702,
     "perm_mean": 0.4944, "perm_std": 0.0658, "perm_p": 0.1194, "n_perms": 200, "significant": False},
    {"variant": "V2 Stricter prevalence (prev=20%, top-100)","n": 165, "auc": 0.5933,
     "perm_mean": 0.4927, "perm_std": 0.0657, "perm_p": 0.0796, "n_perms": 200, "significant": False},
    {"variant": "V3 Top-50 features (prev=10%, top-50)",    "n": 165, "auc": 0.5619,
     "perm_mean": 0.4957, "perm_std": 0.0834, "perm_p": 0.1841, "n_perms": 200, "significant": False},
]

print(f"[{time.strftime('%H:%M:%S')}] Loading data …")
clr      = pd.read_csv(CLR_PATH,    sep="\t", index_col="run_accession")
raw      = pd.read_csv(RAW_PATH,    sep="\t", index_col="run_accession")
labels   = pd.read_csv(LABELS_PATH, sep="\t").set_index("run_accession")
sites    = pd.read_csv(SITES_PATH,  sep="\t").set_index("run_accession")
response = labels["response"]
print(f"  n={len(clr)}, R={(response=='R').sum()}, NR={(response=='NR').sum()}")


def feature_select(clr_tr, raw_tr, resp_tr, prev_frac, top_n):
    n = len(clr_tr)
    min_prev  = max(2, int(np.ceil(prev_frac * n)))
    pres      = (raw_tr > 0).sum(axis=0)
    prevalent = pres[pres >= min_prev].index.tolist()
    vv        = clr_tr[prevalent].var(axis=0)
    high_var  = vv[vv >= vv.quantile(1.0 - DEFAULT_VAR_FRAC)].index.tolist()
    y_bin     = (resp_tr == "R").astype(int).values
    pbs       = {col: abs(stats.pointbiserialr(y_bin, clr_tr[col].values)[0])
                 for col in high_var}
    return pd.Series(pbs).sort_values(ascending=False).head(top_n).index.tolist()


def run_loocv(clr_df, raw_df, resp_ser, prev_frac=DEFAULT_PREV_FRAC, top_n=DEFAULT_TOP_N):
    samples = clr_df.index.tolist()
    y_true, y_prob = [], []
    for held_out in samples:
        train_idx = [s for s in samples if s != held_out]
        tr_clr = clr_df.loc[train_idx]; tr_raw = raw_df.loc[train_idx]
        tr_resp = resp_ser.loc[train_idx]; te_clr = clr_df.loc[[held_out]]
        sel = feature_select(tr_clr, tr_raw, tr_resp, prev_frac, top_n)
        m = LogisticRegression(**EN_KWARGS)
        m.fit(tr_clr[sel].values, tr_resp.values)
        classes = list(m.classes_)
        y_true.append(1 if resp_ser[held_out] == "R" else 0)
        y_prob.append(float(m.predict_proba(te_clr[sel].values)[0, classes.index("R")]))
    return np.array(y_true), np.array(y_prob)


def permutation_test(clr_df, raw_df, resp_ser, obs_auc, n_perms,
                     prev_frac=DEFAULT_PREV_FRAC, top_n=DEFAULT_TOP_N):
    samples = clr_df.index.tolist()
    rng = np.random.default_rng(SEED + 200)
    perm_aucs = []
    for perm_i in range(n_perms):
        perm_resp = resp_ser.copy()
        perm_resp.iloc[:] = rng.permutation(resp_ser.values)
        y_t, y_p = run_loocv(clr_df, raw_df, perm_resp, prev_frac, top_n)
        perm_aucs.append(roc_auc_score(y_t, y_p) if len(np.unique(y_t)) == 2 else 0.5)
        if (perm_i + 1) % 20 == 0:
            print(f"    perm {perm_i+1}/{n_perms}  mean_perm_auc={np.mean(perm_aucs):.3f}", flush=True)
    return perm_aucs


def report_results(variant_name, y_true, y_prob, perm_aucs, obs_auc, n):
    p_val = float((np.array(perm_aucs) >= obs_auc).sum() + 1) / (len(perm_aucs) + 1)
    sig   = "(*)" if p_val < 0.05 else "(ns)"
    print(f"\n  ── {variant_name} ──")
    print(f"     n={n}  AUC={obs_auc:.4f}  perm_mean={np.mean(perm_aucs):.4f}"
          f"±{np.std(perm_aucs):.4f}  p={p_val:.4f} {sig}", flush=True)
    return {"variant": variant_name, "n": n, "auc": round(obs_auc, 4),
            "perm_mean": round(float(np.mean(perm_aucs)), 4),
            "perm_std":  round(float(np.std(perm_aucs)), 4),
            "perm_p": round(p_val, 4), "n_perms": len(perm_aucs),
            "significant": p_val < 0.05}


def _parse_genus_name(stripped):
    if stripped.startswith("Candidatus "):
        return "Candidatus_" + stripped.split()[1]
    elif " (" in stripped:
        return stripped.split(" (")[0]
    return stripped.split()[0]


def build_rarefied_clr(target=RAREFY_TARGET, seed=SEED):
    print(f"  Building rarefied CLR (target={target:,} reads/sample) …", flush=True)
    known_genera = set(clr.columns)
    sample_counts = {}
    for fp in sorted(glob.glob(f"{REPORTS_DIR}/*_report.txt")):
        sample = os.path.basename(fp).replace("_report.txt", "")
        genus_reads = {}
        with open(fp) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6 or parts[3] != "G":
                    continue
                try:
                    clade_reads = int(parts[1])
                except ValueError:
                    continue
                name = _parse_genus_name(parts[5].strip())
                if name in EXCLUDED_GENERA or not name:
                    continue
                genus_reads[name] = genus_reads.get(name, 0) + clade_reads
        sample_counts[sample] = genus_reads

    count_df = (pd.DataFrame(sample_counts).T
                  .fillna(0)                                      # fill within-df NaN (absent genera)
                  .reindex(columns=sorted(known_genera), fill_value=0)
                  .reindex(clr.index, fill_value=0))

    row_totals = count_df.sum(axis=1)
    print(f"  Genus-level read depth: min={row_totals.min():,.0f}  "
          f"mean={row_totals.mean():,.0f}  max={row_totals.max():,.0f}", flush=True)
    rarefy_depth = target if (row_totals >= target).all() else int(row_totals.min())
    if rarefy_depth < target:
        print(f"  WARNING: using min depth {rarefy_depth:,} (some samples below target)", flush=True)

    rng = np.random.default_rng(seed)
    rare_counts = np.zeros_like(count_df.values, dtype=float)
    for i in range(len(count_df)):
        row = count_df.values[i].astype(float)
        row = np.maximum(row, 0.0)          # clip any floating-point negatives
        total = row.sum()
        if total <= 0:
            continue
        probs = row / total
        probs = np.maximum(probs, 0.0)      # clip again after division
        probs /= probs.sum()                # renormalize to exact 1.0
        rare_counts[i] = rng.multinomial(rarefy_depth, probs)

    rare_df = pd.DataFrame(rare_counts, index=count_df.index, columns=count_df.columns)
    X = rare_df.values + 1.0
    lx = np.log(X)
    clr_vals = lx - lx.mean(axis=1, keepdims=True)
    return (pd.DataFrame(clr_vals, index=rare_df.index, columns=rare_df.columns), rare_df)


# =============================================================================
# Run V4 and V5
# =============================================================================

new_results = []

# ── Variant 4: Rarefaction ────────────────────────────────────────────────────
print(f"\n[{time.strftime('%H:%M:%S')}] Variant 4: Rarefaction to {RAREFY_TARGET:,} reads …")
t0 = time.time()
rare_clr, rare_raw = build_rarefied_clr(target=RAREFY_TARGET, seed=SEED)
y_t4, y_p4 = run_loocv(rare_clr, rare_raw, response)
auc4 = roc_auc_score(y_t4, y_p4)
print(f"  LOOCV done in {time.time()-t0:.0f}s  AUC={auc4:.4f}  Running permutation test …", flush=True)
t0 = time.time()
perm4 = permutation_test(rare_clr, rare_raw, response, auc4, N_PERMS)
print(f"  Permutation test done in {time.time()-t0:.0f}s")
r4 = report_results("V4 Rarefaction (1M reads, prev=10%, top-100)", y_t4, y_p4, perm4, auc4, len(clr))
new_results.append(r4)

# ── Variant 5: Subset (drop Barcelona + Leeds) ────────────────────────────────
print(f"\n[{time.strftime('%H:%M:%S')}] Variant 5: Subset (drop Barcelona + Leeds, n=135) …")
keep_idx = sites[sites["site"].isin({"PRIMM-UK", "PRIMM-NL", "Manchester"})].index.intersection(clr.index)
clr5, raw5, resp5 = clr.loc[keep_idx], raw.loc[keep_idx], response.loc[keep_idx]
print(f"  Subset: n={len(clr5)}, R={(resp5=='R').sum()}, NR={(resp5=='NR').sum()}")
print(f"  Sites: {sites.loc[keep_idx,'site'].value_counts().to_dict()}")
t0 = time.time()
y_t5, y_p5 = run_loocv(clr5, raw5, resp5)
auc5 = roc_auc_score(y_t5, y_p5)
print(f"  LOOCV done in {time.time()-t0:.0f}s  AUC={auc5:.4f}  Running permutation test …", flush=True)
t0 = time.time()
perm5 = permutation_test(clr5, raw5, resp5, auc5, N_PERMS)
print(f"  Permutation test done in {time.time()-t0:.0f}s")
r5 = report_results("V5 Subset PRIMM-UK+NL+MCR (n=135)", y_t5, y_p5, perm5, auc5, len(clr5))
new_results.append(r5)

# =============================================================================
# Combined summary (V1–V5)
# =============================================================================

all_results = PRIOR_RESULTS + new_results

print(f"\n{'='*65}")
print("SENSITIVITY ANALYSIS SUMMARY — C4 Lee 2022, ElasticNet")
print(f"Fixed EN: C={EN_C}, l1_ratio={EN_L1} | Perms: N={N_PERMS}")
print(f"{'='*65}")
print(f"  {'Variant':<46} {'n':>4} {'AUC':>6} {'p':>7} {'Sig':>5}")
print(f"  {'-'*46} {'-'*4} {'-'*6} {'-'*7} {'-'*5}")
for r in all_results:
    sig = "(*)" if r["significant"] else "(ns)"
    print(f"  {r['variant']:<46} {r['n']:>4} {r['auc']:>6.4f} {r['perm_p']:>7.4f} {sig:>5}")

pd.DataFrame(all_results).to_csv(f"{OUT_DIR}/sensitivity_summary.tsv", sep="\t", index=False)
pd.DataFrame({"V4_rarefy": perm4, "V5_subset": perm5}).to_csv(
    f"{OUT_DIR}/sensitivity_perm_aucs_v4v5.tsv", sep="\t", index=False)

prob_df = pd.read_csv("results/ml/lee2022/X_genus_clr.tsv", sep="\t",
                      index_col="run_accession", usecols=["run_accession"])
pd.DataFrame({"run_accession": clr.index, "actual": response.values,
              "V4_prob_R": y_p4}).to_csv(
    f"{OUT_DIR}/sensitivity_probs_v4_n165.tsv", sep="\t", index=False)
pd.DataFrame({"run_accession": clr5.index, "actual": resp5.values,
              "V5_prob_R": y_p5}).to_csv(
    f"{OUT_DIR}/sensitivity_probs_v5_n135.tsv", sep="\t", index=False)

print(f"\nSaved to {OUT_DIR}/")
print(f"[{time.strftime('%H:%M:%S')}] DONE")
