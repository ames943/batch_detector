#!/usr/bin/env python3
"""
Sensitivity analysis on Lee 2022 (C4, n=165): 5 preprocessing variants.
ElasticNet LOOCV + permutation test (N=200) for each variant.
Fixed EN hyperparameters (C=1.0, l1_ratio=0.7 — consensus from baseline nested CV).

Run from cancer_project/:
    python3 scripts/c4_sensitivity_analysis.py 2>&1 | tee logs/c4_sensitivity.log
"""

import os, glob, time, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings("ignore", category=ConvergenceWarning)

# ── Paths ─────────────────────────────────────────────────────────────────────
CLR_PATH      = "results/ml/lee2022/X_genus_clr.tsv"
RAW_PATH      = "results/ml/lee2022/X_genus_raw.tsv"
LABELS_PATH   = "metadata/lee2022_labels.tsv"
SITES_PATH    = "metadata/lee2022_sites.tsv"
REPORTS_DIR   = "results/kraken_reports/lee2022"
OUT_DIR       = "results/ml/lee2022/sensitivity"
os.makedirs(OUT_DIR, exist_ok=True)

# ── Fixed hyperparameters (consensus from baseline nested CV) ─────────────────
EN_C       = 1.0
EN_L1      = 0.7
EN_KWARGS  = dict(
    penalty="elasticnet", solver="saga",
    C=EN_C, l1_ratio=EN_L1,
    class_weight="balanced",
    max_iter=5000, tol=1e-3, random_state=42,
)

# ── Default feature-selection parameters ──────────────────────────────────────
DEFAULT_PREV_FRAC  = 0.10
DEFAULT_VAR_FRAC   = 0.50
DEFAULT_TOP_N      = 100
SEED               = 42
N_PERMS            = 200
RAREFY_TARGET      = 1_000_000
EXCLUDED_GENERA    = {"Homo"}

print(f"[{time.strftime('%H:%M:%S')}] Loading base data …")
clr    = pd.read_csv(CLR_PATH,    sep="\t", index_col="run_accession")
raw    = pd.read_csv(RAW_PATH,    sep="\t", index_col="run_accession")
labels = pd.read_csv(LABELS_PATH, sep="\t").set_index("run_accession")
sites  = pd.read_csv(SITES_PATH,  sep="\t").set_index("run_accession")

response = labels["response"]
print(f"  n={len(clr)}, R={(response=='R').sum()}, NR={(response=='NR').sum()}")
print(f"  CLR features: {clr.shape[1]}")


# =============================================================================
# Core LOOCV helpers
# =============================================================================

def feature_select(clr_tr, raw_tr, resp_tr, prev_frac, top_n):
    """Leak-free 3-step feature selection on training fold."""
    n = len(clr_tr)
    min_prev  = max(2, int(np.ceil(prev_frac * n)))
    pres      = (raw_tr > 0).sum(axis=0)
    prevalent = pres[pres >= min_prev].index.tolist()
    vv        = clr_tr[prevalent].var(axis=0)
    var_frac  = DEFAULT_VAR_FRAC
    high_var  = vv[vv >= vv.quantile(1.0 - var_frac)].index.tolist()
    y_bin     = (resp_tr == "R").astype(int).values
    pbs       = {col: abs(stats.pointbiserialr(y_bin, clr_tr[col].values)[0])
                 for col in high_var}
    return pd.Series(pbs).sort_values(ascending=False).head(top_n).index.tolist()


def run_loocv(clr_df, raw_df, resp_ser, prev_frac=DEFAULT_PREV_FRAC,
              top_n=DEFAULT_TOP_N):
    """LOOCV with fixed EN + per-fold feature selection. Returns (y_true, y_prob)."""
    samples  = clr_df.index.tolist()
    y_true   = []
    y_prob   = []

    for held_out in samples:
        train_idx = [s for s in samples if s != held_out]
        tr_clr    = clr_df.loc[train_idx]
        tr_raw    = raw_df.loc[train_idx]
        tr_resp   = resp_ser.loc[train_idx]
        te_clr    = clr_df.loc[[held_out]]

        sel    = feature_select(tr_clr, tr_raw, tr_resp, prev_frac, top_n)
        X_tr   = tr_clr[sel].values
        y_tr   = tr_resp.values
        X_te   = te_clr[sel].values

        model = LogisticRegression(**EN_KWARGS)
        model.fit(X_tr, y_tr)
        classes  = list(model.classes_)
        r_prob   = float(model.predict_proba(X_te)[0, classes.index("R")])

        y_true.append(1 if resp_ser[held_out] == "R" else 0)
        y_prob.append(r_prob)

    return np.array(y_true), np.array(y_prob)


def permutation_test(clr_df, raw_df, resp_ser, obs_auc, n_perms,
                     prev_frac=DEFAULT_PREV_FRAC, top_n=DEFAULT_TOP_N):
    """Permutation test: shuffle labels, re-run LOOCV. Returns perm AUCs list."""
    samples = clr_df.index.tolist()
    rng = np.random.default_rng(SEED + 200)
    perm_aucs = []

    for perm_i in range(n_perms):
        perm_resp = resp_ser.copy()
        perm_resp.iloc[:] = rng.permutation(resp_ser.values)

        y_t, y_p = run_loocv(clr_df, raw_df, perm_resp, prev_frac, top_n)
        if len(np.unique(y_t)) == 2:
            perm_aucs.append(roc_auc_score(y_t, y_p))
        else:
            perm_aucs.append(0.5)

        if (perm_i + 1) % 20 == 0:
            print(f"    perm {perm_i+1}/{n_perms}  mean_perm_auc={np.mean(perm_aucs):.3f}",
                  flush=True)

    return perm_aucs


def report_results(variant_name, y_true, y_prob, perm_aucs, obs_auc, n):
    p_val = float((np.array(perm_aucs) >= obs_auc).sum() + 1) / (len(perm_aucs) + 1)
    sig   = "(*)" if p_val < 0.05 else "(ns)"
    print(f"\n  ── {variant_name} ──")
    print(f"     n={n}  AUC={obs_auc:.4f}  perm_mean={np.mean(perm_aucs):.4f}"
          f"±{np.std(perm_aucs):.4f}  p={p_val:.4f} {sig}")
    return {"variant": variant_name, "n": n,
            "auc": round(obs_auc, 4),
            "perm_mean": round(float(np.mean(perm_aucs)), 4),
            "perm_std":  round(float(np.std(perm_aucs)), 4),
            "perm_p": round(p_val, 4),
            "n_perms": len(perm_aucs),
            "significant": p_val < 0.05}


# =============================================================================
# Variant 4 helper: rarefaction from Kraken2 reports
# =============================================================================

def _parse_genus_name(stripped):
    """Apply same genus-name normalisation as build_matrix.py."""
    if stripped.startswith("Candidatus "):
        return "Candidatus_" + stripped.split()[1]
    elif " (" in stripped:
        return stripped.split(" (")[0]
    else:
        return stripped.split()[0]


def build_rarefied_clr(target=RAREFY_TARGET, seed=SEED):
    """
    Re-read Kraken2 reports, extract genus clade reads, rarefy to 'target' reads,
    return CLR DataFrame aligned to the existing CLR matrix column vocabulary.
    """
    print(f"  Building rarefied CLR (target={target:,} reads/sample) …", flush=True)
    report_files = sorted(glob.glob(f"{REPORTS_DIR}/*_report.txt"))
    known_genera = set(clr.columns)

    sample_counts = {}
    for fp in report_files:
        sample = os.path.basename(fp).replace("_report.txt", "")
        genus_reads = {}
        with open(fp) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                if parts[3] != "G":
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

    # Align to known genera (same vocabulary as CLR matrix)
    count_df = pd.DataFrame(sample_counts).T.reindex(
        columns=sorted(known_genera), fill_value=0
    )
    count_df = count_df.reindex(clr.index)   # same sample order

    # Check minimum depth
    row_totals = count_df.sum(axis=1)
    print(f"  Genus-level read depth: min={row_totals.min():,.0f}  "
          f"mean={row_totals.mean():,.0f}  max={row_totals.max():,.0f}", flush=True)
    n_below = (row_totals < target).sum()
    if n_below > 0:
        print(f"  WARNING: {n_below} samples have < {target:,} genus reads. "
              f"Using min depth = {int(row_totals.min()):,} instead.", flush=True)
        rarefy_depth = int(row_totals.min())
    else:
        rarefy_depth = target

    # Rarefy using multinomial sampling
    rng = np.random.default_rng(seed)
    rare_counts = np.zeros_like(count_df.values, dtype=float)
    for i in range(len(count_df)):
        row = count_df.values[i]
        total = row.sum()
        if total == 0:
            continue
        probs = row / total
        rare_counts[i] = rng.multinomial(rarefy_depth, probs)

    rare_df = pd.DataFrame(rare_counts, index=count_df.index, columns=count_df.columns)

    # CLR with pseudocount 1 (integer counts after rarefaction)
    X = rare_df.values + 1.0
    lx = np.log(X)
    clr_vals = lx - lx.mean(axis=1, keepdims=True)
    rare_clr = pd.DataFrame(clr_vals, index=rare_df.index, columns=rare_df.columns)

    # Raw rarefied for feature selection (prevalence filter needs 0/nonzero)
    return rare_clr, rare_df


# =============================================================================
# Run all variants
# =============================================================================

all_results = []

print(f"\n{'='*65}")
print("C4 Lee 2022 — ElasticNet Sensitivity Analysis")
print(f"EN hyperparameters: C={EN_C}, l1_ratio={EN_L1}  (consensus from baseline nested CV)")
print(f"Permutation test: N={N_PERMS}")
print(f"{'='*65}")

# ── Variant 1: Baseline (re-run with fixed EN for comparability) ───────────────
print(f"\n[{time.strftime('%H:%M:%S')}] Variant 1: Baseline (prev=10%, top-100) …")
t0 = time.time()
y_t1, y_p1 = run_loocv(clr, raw, response,
                        prev_frac=0.10, top_n=100)
auc1 = roc_auc_score(y_t1, y_p1)
print(f"  LOOCV done in {time.time()-t0:.0f}s  AUC={auc1:.4f}  Running permutation test …")
t0 = time.time()
perm1 = permutation_test(clr, raw, response, auc1, N_PERMS, prev_frac=0.10, top_n=100)
print(f"  Permutation test done in {time.time()-t0:.0f}s")
r1 = report_results("V1 Baseline (prev=10%, top-100)", y_t1, y_p1, perm1, auc1, len(clr))
all_results.append(r1)

# ── Variant 2: Stricter prevalence filter (≥20%) ──────────────────────────────
print(f"\n[{time.strftime('%H:%M:%S')}] Variant 2: Stricter prevalence (prev=20%, top-100) …")
t0 = time.time()
y_t2, y_p2 = run_loocv(clr, raw, response,
                        prev_frac=0.20, top_n=100)
auc2 = roc_auc_score(y_t2, y_p2)
print(f"  LOOCV done in {time.time()-t0:.0f}s  AUC={auc2:.4f}  Running permutation test …")
t0 = time.time()
perm2 = permutation_test(clr, raw, response, auc2, N_PERMS, prev_frac=0.20, top_n=100)
print(f"  Permutation test done in {time.time()-t0:.0f}s")
r2 = report_results("V2 Stricter prevalence (prev=20%, top-100)", y_t2, y_p2, perm2, auc2, len(clr))
all_results.append(r2)

# ── Variant 3: Top 50 features ────────────────────────────────────────────────
print(f"\n[{time.strftime('%H:%M:%S')}] Variant 3: Top-50 features (prev=10%, top-50) …")
t0 = time.time()
y_t3, y_p3 = run_loocv(clr, raw, response,
                        prev_frac=0.10, top_n=50)
auc3 = roc_auc_score(y_t3, y_p3)
print(f"  LOOCV done in {time.time()-t0:.0f}s  AUC={auc3:.4f}  Running permutation test …")
t0 = time.time()
perm3 = permutation_test(clr, raw, response, auc3, N_PERMS, prev_frac=0.10, top_n=50)
print(f"  Permutation test done in {time.time()-t0:.0f}s")
r3 = report_results("V3 Top-50 features (prev=10%, top-50)", y_t3, y_p3, perm3, auc3, len(clr))
all_results.append(r3)

# ── Variant 4: Rarefaction to 1M reads ────────────────────────────────────────
print(f"\n[{time.strftime('%H:%M:%S')}] Variant 4: Rarefaction to {RAREFY_TARGET:,} reads …")
t0 = time.time()
rare_clr, rare_raw = build_rarefied_clr(target=RAREFY_TARGET, seed=SEED)
# Use rarefied CLR for both feature values and for feature selection
# (for prevalence filter on rarefied raw counts)
y_t4, y_p4 = run_loocv(rare_clr, rare_raw, response,
                        prev_frac=0.10, top_n=100)
auc4 = roc_auc_score(y_t4, y_p4)
print(f"  LOOCV done in {time.time()-t0:.0f}s  AUC={auc4:.4f}  Running permutation test …")
t0 = time.time()
perm4 = permutation_test(rare_clr, rare_raw, response, auc4, N_PERMS, prev_frac=0.10, top_n=100)
print(f"  Permutation test done in {time.time()-t0:.0f}s")
r4 = report_results("V4 Rarefaction (1M reads, prev=10%, top-100)", y_t4, y_p4, perm4, auc4, len(clr))
all_results.append(r4)

# ── Variant 5: Subset — remove Barcelona (n=12) + Leeds (n=18) ───────────────
print(f"\n[{time.strftime('%H:%M:%S')}] Variant 5: Subset (drop Barcelona + Leeds, n=135) …")
keep_sites = {"PRIMM-UK", "PRIMM-NL", "Manchester"}
keep_mask  = sites["site"].isin(keep_sites)
keep_idx   = sites[keep_mask].index.intersection(clr.index)

clr5  = clr.loc[keep_idx]
raw5  = raw.loc[keep_idx]
resp5 = response.loc[keep_idx]
print(f"  Subset: n={len(clr5)}, R={(resp5=='R').sum()}, NR={(resp5=='NR').sum()}")
print(f"  Site distribution: {sites.loc[keep_idx, 'site'].value_counts().to_dict()}")

t0 = time.time()
y_t5, y_p5 = run_loocv(clr5, raw5, resp5,
                        prev_frac=0.10, top_n=100)
auc5 = roc_auc_score(y_t5, y_p5)
print(f"  LOOCV done in {time.time()-t0:.0f}s  AUC={auc5:.4f}  Running permutation test …")
t0 = time.time()
perm5 = permutation_test(clr5, raw5, resp5, auc5, N_PERMS, prev_frac=0.10, top_n=100)
print(f"  Permutation test done in {time.time()-t0:.0f}s")
r5 = report_results("V5 Subset PRIMM-UK+NL+MCR (n=135)", y_t5, y_p5, perm5, auc5, len(clr5))
all_results.append(r5)

# =============================================================================
# Summary table
# =============================================================================

print(f"\n{'='*65}")
print("SENSITIVITY ANALYSIS SUMMARY — C4 Lee 2022, ElasticNet")
print(f"Fixed EN: C={EN_C}, l1_ratio={EN_L1} | Perms: N={N_PERMS}")
print(f"{'='*65}")
print(f"  {'Variant':<46} {'n':>4} {'AUC':>6} {'p':>7} {'Sig':>5}")
print(f"  {'-'*46} {'-'*4} {'-'*6} {'-'*7} {'-'*5}")
for r in all_results:
    sig = "(*)" if r["significant"] else "(ns)"
    print(f"  {r['variant']:<46} {r['n']:>4} {r['auc']:>6.4f} {r['perm_p']:>7.4f} {sig:>5}")

# ── Save results ──────────────────────────────────────────────────────────────
summary_df = pd.DataFrame(all_results)
summary_path = f"{OUT_DIR}/sensitivity_summary.tsv"
summary_df.to_csv(summary_path, sep="\t", index=False)

# Save per-variant perm AUCs for inspection
perm_df = pd.DataFrame({
    "V1_baseline": perm1,
    "V2_prev20":   perm2,
    "V3_top50":    perm3,
    "V4_rarefy":   perm4,
    "V5_subset":   perm5,
})
perm_df.to_csv(f"{OUT_DIR}/sensitivity_perm_aucs.tsv", sep="\t", index=False)

# Save per-sample predicted probabilities for each variant
prob_df = pd.DataFrame({
    "run_accession": clr.index.tolist(),
    "actual": response.values,
    "V1_prob_R": y_p1,
    "V2_prob_R": y_p2,
    "V3_prob_R": y_p3,
    "V4_prob_R": y_p4,
})
# V5 has different sample set — save separately
prob_df.to_csv(f"{OUT_DIR}/sensitivity_probs_n165.tsv", sep="\t", index=False)

prob_v5 = pd.DataFrame({
    "run_accession": clr5.index.tolist(),
    "actual": resp5.values,
    "V5_prob_R": y_p5,
})
prob_v5.to_csv(f"{OUT_DIR}/sensitivity_probs_v5_n135.tsv", sep="\t", index=False)

print(f"\nSaved:\n  {summary_path}")
print(f"  {OUT_DIR}/sensitivity_perm_aucs.tsv")
print(f"  {OUT_DIR}/sensitivity_probs_n165.tsv")
print(f"  {OUT_DIR}/sensitivity_probs_v5_n135.tsv")
print(f"\n[{time.strftime('%H:%M:%S')}] DONE")
