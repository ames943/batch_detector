#!/usr/bin/env python3
"""
Task 5 — Repeated stratified 5-fold CV (new primary protocol), replacing LOOCV.

10 repeats x stratified 5-fold, inner 5-fold hyperparameter tuning, same
leak-free pipeline as scripts/nested_cv_n283.py (per-fold prevalence->variance
->top-100-|point-biserial-r| feature selection; per-fold mean-centering batch
correction for multi-cohort datasets). Models: ElasticNet, RandomForest, XGBoost.

Permutation test (N>=100): re-tuning the full inner grid under every permuted
label draw would multiply cost ~100x for negligible benefit, so — following
the precedent already set in this codebase (scripts/nested_cv_n283.py's own
permutation test) — the null uses FIXED hyperparameters: the consensus
(most-frequent) best params from the observed run's outer folds. Feature
selection and batch correction are still re-done per fold per permutation
(no leakage). This is documented explicitly, not silently substituted.

Usage:
    python3 scripts/revision/task5_repeated_cv.py --dataset {c1,n118,n283,c4}
"""
import argparse, os, sys, time, json, collections, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning
import xgboost as xgb

warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

SEED = 42
N_REPEATS = 10
OUTER_K = 5
INNER_K = 5
N_PERMS = 100
PREVALENCE_FRAC = 0.10
VARIANCE_KEEP_FRACTION = 0.50
TOP_N = 100
RF_N_JOBS = 2

EN_GRID = [{"C": round(1.0 / a, 6), "l1_ratio": l1} for a in [0.01, 0.1, 1.0] for l1 in [0.3, 0.5, 0.7]]
RF_GRID = [{"n_estimators": ne, "max_depth": md} for ne in [200, 500] for md in [5, 10, None]]
XGB_GRID = [{"n_estimators": ne, "max_depth": md, "learning_rate": lr}
            for ne in [100, 200] for md in [2, 3] for lr in [0.05, 0.1]]


def cohort_of_n283(sid):
    if sid.startswith("SRR5930"): return "cohort1"
    if sid.startswith("SRR11413"): return "cohort2"
    if sid.startswith("SRR6000"): return "cohort3"
    return "cohort4"


def load_dataset(name):
    if name == "c1":
        clr = pd.read_csv("results/ml/n118_3cohort/X_genus_clr.tsv", sep="\t", index_col="run_accession")
        raw = pd.read_csv("results/ml/n118_3cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
        labels = pd.read_csv("metadata/response_labels_3cohort.tsv", sep="\t").set_index("run_accession")
        mask = clr.index.str.startswith("SRR5930")
        clr, raw, labels = clr.loc[mask], raw.loc[mask], labels.reindex(clr.loc[mask].index)
        batch = None
    elif name == "n118":
        clr = pd.read_csv("results/ml/n118_3cohort/X_genus_clr.tsv", sep="\t", index_col="run_accession")
        raw = pd.read_csv("results/ml/n118_3cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
        labels = pd.read_csv("metadata/response_labels_3cohort.tsv", sep="\t").set_index("run_accession").reindex(clr.index)
        batch = labels["cohort"]
    elif name == "n283":
        clr = pd.read_csv("results/ml/n283_4cohort/X_genus_clr.tsv", sep="\t", index_col="run_accession")
        raw = pd.read_csv("results/ml/n283_4cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
        labels = pd.read_csv("results/ml/n283_4cohort/response_labels_n283.tsv", sep="\t").set_index("run_accession").reindex(clr.index)
        batch = pd.Series([cohort_of_n283(s) for s in clr.index], index=clr.index)
    elif name == "c4":
        clr = pd.read_csv("results/ml/lee2022/X_genus_clr.tsv", sep="\t", index_col="run_accession")
        raw = pd.read_csv("results/ml/lee2022/X_genus_raw.tsv", sep="\t", index_col="run_accession")
        labels = pd.read_csv("metadata/lee2022_labels.tsv", sep="\t").set_index("run_accession").reindex(clr.index)
        batch = None
    else:
        raise ValueError(name)
    keep = labels["response"].notna()
    clr, raw, labels = clr.loc[keep], raw.loc[keep], labels.loc[keep]
    if batch is not None:
        batch = batch.loc[keep]
    return clr, raw, labels["response"], batch


def batch_correct_fold(tr_clr, te_clr, tr_batch, te_batch):
    if tr_batch is None:
        return tr_clr.values, te_clr.values
    gm = tr_clr.mean(axis=0)
    offsets = {b: tr_clr.loc[tr_batch == b].mean(axis=0) - gm for b in tr_batch.unique()}
    corr_tr = tr_clr.copy()
    for b, off in offsets.items():
        corr_tr.loc[tr_batch == b] = tr_clr.loc[tr_batch == b].values - off.values
    corr_te = te_clr.copy()
    for b in te_batch.unique():
        m = te_batch == b
        if b in offsets:
            corr_te.loc[m] = te_clr.loc[m].values - offsets[b].values
    return corr_tr.values, corr_te.values


def feature_select(clr_tr_df, raw_tr, resp_tr):
    n = len(clr_tr_df)
    min_prev = max(2, int(np.ceil(PREVALENCE_FRAC * n)))
    pres = (raw_tr > 0).sum(axis=0)
    prevalent = pres[pres >= min_prev].index.tolist()
    vv = clr_tr_df[prevalent].var(axis=0)
    high_var = vv[vv >= vv.quantile(1.0 - VARIANCE_KEEP_FRACTION)].index.tolist()
    y_bin = (resp_tr == "R").astype(int).values
    pbs = {c: abs(stats.pointbiserialr(y_bin, clr_tr_df[c].values)[0]) for c in high_var}
    return pd.Series(pbs).sort_values(ascending=False).head(TOP_N).index.tolist()


def fit_predict_en(X_tr, y_tr, X_te, params):
    m = LogisticRegression(penalty="elasticnet", solver="saga", C=params["C"],
                            l1_ratio=params["l1_ratio"], class_weight="balanced",
                            max_iter=3000, tol=1e-3, random_state=SEED)
    m.fit(X_tr, y_tr)
    idx = list(m.classes_).index("R")
    return m.predict_proba(X_te)[:, idx]


def fit_predict_rf(X_tr, y_tr, X_te, params):
    m = RandomForestClassifier(n_estimators=params["n_estimators"], max_depth=params["max_depth"],
                                class_weight="balanced", random_state=SEED, n_jobs=RF_N_JOBS)
    m.fit(X_tr, y_tr)
    idx = list(m.classes_).index("R")
    return m.predict_proba(X_te)[:, idx]


def fit_predict_xgb(X_tr, y_tr_bin, X_te, params):
    m = xgb.XGBClassifier(n_estimators=params["n_estimators"], max_depth=params["max_depth"],
                           learning_rate=params["learning_rate"], eval_metric="logloss",
                           n_jobs=1, verbosity=0, random_state=SEED)
    m.fit(X_tr, y_tr_bin)
    return m.predict_proba(X_te)[:, 1]


def inner_tune(clr_tr_df, raw_tr, resp_tr, batch_tr, fold_seed):
    y_bin = (resp_tr == "R").astype(int).values
    skf = StratifiedKFold(n_splits=INNER_K, shuffle=True, random_state=fold_seed)
    en_aucs = [[] for _ in EN_GRID]; rf_aucs = [[] for _ in RF_GRID]; xgb_aucs = [[] for _ in XGB_GRID]

    for itr, iva in skf.split(clr_tr_df.values, y_bin):
        i_clr_tr, i_clr_va = clr_tr_df.iloc[itr], clr_tr_df.iloc[iva]
        i_raw_tr = raw_tr.iloc[itr]
        i_resp_tr, i_resp_va = resp_tr.iloc[itr], resp_tr.iloc[iva]
        i_bat_tr = batch_tr.iloc[itr] if batch_tr is not None else None
        i_bat_va = batch_tr.iloc[iva] if batch_tr is not None else None

        c_tr, c_va = batch_correct_fold(i_clr_tr, i_clr_va, i_bat_tr, i_bat_va)
        c_tr_df = pd.DataFrame(c_tr, index=i_clr_tr.index, columns=i_clr_tr.columns)
        sel = feature_select(c_tr_df, i_raw_tr, i_resp_tr)
        feat_idx = [list(clr_tr_df.columns).index(f) for f in sel]
        X_tr, X_va = c_tr_df[sel].values, c_va[:, feat_idx]
        y_va_bin = (i_resp_va == "R").astype(int).values
        if len(np.unique(y_va_bin)) < 2:
            continue

        for gi, p in enumerate(EN_GRID):
            try:
                probs = fit_predict_en(X_tr, i_resp_tr.values, X_va, p)
                en_aucs[gi].append(roc_auc_score(y_va_bin, probs))
            except Exception:
                pass
        for gi, p in enumerate(RF_GRID):
            try:
                probs = fit_predict_rf(X_tr, i_resp_tr.values, X_va, p)
                rf_aucs[gi].append(roc_auc_score(y_va_bin, probs))
            except Exception:
                pass
        y_tr_bin = (i_resp_tr == "R").astype(int).values
        for gi, p in enumerate(XGB_GRID):
            try:
                probs = fit_predict_xgb(X_tr, y_tr_bin, X_va, p)
                xgb_aucs[gi].append(roc_auc_score(y_va_bin, probs))
            except Exception:
                pass

    def best(grid, aucs):
        means = [np.mean(a) if a else -1.0 for a in aucs]
        return grid[int(np.argmax(means))]

    return best(EN_GRID, en_aucs), best(RF_GRID, rf_aucs), best(XGB_GRID, xgb_aucs)


def run_outer_fold_tuned(clr, raw, resp, batch, tr_idx, te_idx, fold_seed):
    clr_tr, clr_te = clr.iloc[tr_idx], clr.iloc[te_idx]
    raw_tr = raw.iloc[tr_idx]
    resp_tr, resp_te = resp.iloc[tr_idx], resp.iloc[te_idx]
    bat_tr = batch.iloc[tr_idx] if batch is not None else None
    bat_te = batch.iloc[te_idx] if batch is not None else None

    corr_tr, corr_te = batch_correct_fold(clr_tr, clr_te, bat_tr, bat_te)
    corr_tr_df = pd.DataFrame(corr_tr, index=clr_tr.index, columns=clr_tr.columns)
    sel = feature_select(corr_tr_df, raw_tr, resp_tr)
    feat_idx = [list(clr.columns).index(f) for f in sel]
    X_tr, X_te = corr_tr_df[sel].values, corr_te[:, feat_idx]

    best_en, best_rf, best_xgb = inner_tune(corr_tr_df, raw_tr, resp_tr, bat_tr, fold_seed)

    y_te_bin = (resp_te == "R").astype(int).values
    p_en = fit_predict_en(X_tr, resp_tr.values, X_te, best_en)
    p_rf = fit_predict_rf(X_tr, resp_tr.values, X_te, best_rf)
    p_xgb = fit_predict_xgb(X_tr, (resp_tr == "R").astype(int).values, X_te, best_xgb)

    return dict(y_true=y_te_bin, p_en=p_en, p_rf=p_rf, p_xgb=p_xgb,
                best_en=best_en, best_rf=best_rf, best_xgb=best_xgb)


def run_outer_fold_fixed(clr, raw, resp, batch, tr_idx, te_idx, params_en, params_rf, params_xgb):
    clr_tr, clr_te = clr.iloc[tr_idx], clr.iloc[te_idx]
    raw_tr = raw.iloc[tr_idx]
    resp_tr, resp_te = resp.iloc[tr_idx], resp.iloc[te_idx]
    bat_tr = batch.iloc[tr_idx] if batch is not None else None
    bat_te = batch.iloc[te_idx] if batch is not None else None

    corr_tr, corr_te = batch_correct_fold(clr_tr, clr_te, bat_tr, bat_te)
    corr_tr_df = pd.DataFrame(corr_tr, index=clr_tr.index, columns=clr_tr.columns)
    sel = feature_select(corr_tr_df, raw_tr, resp_tr)
    feat_idx = [list(clr.columns).index(f) for f in sel]
    X_tr, X_te = corr_tr_df[sel].values, corr_te[:, feat_idx]

    y_te_bin = (resp_te == "R").astype(int).values
    p_en = fit_predict_en(X_tr, resp_tr.values, X_te, params_en)
    p_rf = fit_predict_rf(X_tr, resp_tr.values, X_te, params_rf)
    p_xgb = fit_predict_xgb(X_tr, (resp_tr == "R").astype(int).values, X_te, params_xgb)
    return y_te_bin, p_en, p_rf, p_xgb


def consensus(param_list):
    counts = collections.Counter(json.dumps(p, sort_keys=True) for p in param_list)
    return json.loads(counts.most_common(1)[0][0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=["c1", "n118", "n283", "c4"])
    args = ap.parse_args()
    OUT = "results/revision"
    os.makedirs(OUT, exist_ok=True)

    t0 = time.time()
    clr, raw, resp, batch = load_dataset(args.dataset)
    n = len(clr)
    print(f"[{time.strftime('%H:%M:%S')}] dataset={args.dataset} n={n} "
          f"R={int((resp=='R').sum())} NR={int((resp=='NR').sum())} "
          f"batch={'yes' if batch is not None else 'no'}", flush=True)

    y_bin_full = (resp == "R").astype(int).values
    idx_all = np.arange(n)

    fold_aucs = []          # per (repeat, fold, model) AUC
    best_params_log = {"en": [], "rf": [], "xgb": []}
    all_probs = {"en": [], "rf": [], "xgb": []}   # concatenated across all outer test folds (all repeats) for pooled AUC/sanity
    all_true = []

    for rep in range(N_REPEATS):
        skf = StratifiedKFold(n_splits=OUTER_K, shuffle=True, random_state=SEED + rep)
        for fold_i, (tr_idx, te_idx) in enumerate(skf.split(idx_all, y_bin_full)):
            fold_seed = SEED + rep * 100 + fold_i
            res = run_outer_fold_tuned(clr, raw, resp, batch, tr_idx, te_idx, fold_seed)
            y_te = res["y_true"]
            if len(np.unique(y_te)) < 2:
                continue
            auc_en = roc_auc_score(y_te, res["p_en"])
            auc_rf = roc_auc_score(y_te, res["p_rf"])
            auc_xgb = roc_auc_score(y_te, res["p_xgb"])
            fold_aucs.append(dict(repeat=rep, fold=fold_i, n_test=len(te_idx),
                                   auc_en=auc_en, auc_rf=auc_rf, auc_xgb=auc_xgb))
            best_params_log["en"].append(res["best_en"])
            best_params_log["rf"].append(res["best_rf"])
            best_params_log["xgb"].append(res["best_xgb"])

        elapsed = time.time() - t0
        print(f"[{time.strftime('%H:%M:%S')}] repeat {rep+1}/{N_REPEATS} done "
              f"(elapsed={elapsed:.0f}s)", flush=True)

    fa = pd.DataFrame(fold_aucs)
    fa.to_csv(f"{OUT}/task5_{args.dataset}_fold_aucs.tsv", sep="\t", index=False)

    summary = dict(
        dataset=args.dataset, n=n, n_folds=len(fa),
        mean_auc_en=fa["auc_en"].mean(), sd_auc_en=fa["auc_en"].std(),
        mean_auc_rf=fa["auc_rf"].mean(), sd_auc_rf=fa["auc_rf"].std(),
        mean_auc_xgb=fa["auc_xgb"].mean(), sd_auc_xgb=fa["auc_xgb"].std(),
    )
    print(f"\n[{time.strftime('%H:%M:%S')}] OBSERVED (n_folds={len(fa)}):")
    for m in ["en", "rf", "xgb"]:
        print(f"  {m.upper():4s}: mean AUC = {summary[f'mean_auc_{m}']:.4f} "
              f"+/- {summary[f'sd_auc_{m}']:.4f}  ({time.time()-t0:.0f}s elapsed)", flush=True)

    cons_en = consensus(best_params_log["en"])
    cons_rf = consensus(best_params_log["rf"])
    cons_xgb = consensus(best_params_log["xgb"])
    print(f"  Consensus params: EN={cons_en}  RF={cons_rf}  XGB={cons_xgb}", flush=True)

    # ── Permutation test (fixed consensus params; see module docstring) ─────
    print(f"\n[{time.strftime('%H:%M:%S')}] Permutation test (N={N_PERMS}, fixed consensus params)...", flush=True)
    rng = np.random.default_rng(SEED + 777)
    perm_means = {"en": [], "rf": [], "xgb": []}
    t_perm0 = time.time()
    for pi in range(N_PERMS):
        perm_resp = resp.copy()
        perm_resp.iloc[:] = rng.permutation(resp.values)
        y_bin_perm = (perm_resp == "R").astype(int).values
        skf = StratifiedKFold(n_splits=OUTER_K, shuffle=True, random_state=SEED + 1000 + pi)
        aucs_en, aucs_rf, aucs_xgb = [], [], []
        for tr_idx, te_idx in skf.split(idx_all, y_bin_perm):
            y_te, p_en, p_rf, p_xgb = run_outer_fold_fixed(
                clr, raw, perm_resp, batch, tr_idx, te_idx, cons_en, cons_rf, cons_xgb)
            if len(np.unique(y_te)) < 2:
                continue
            aucs_en.append(roc_auc_score(y_te, p_en))
            aucs_rf.append(roc_auc_score(y_te, p_rf))
            aucs_xgb.append(roc_auc_score(y_te, p_xgb))
        perm_means["en"].append(np.mean(aucs_en))
        perm_means["rf"].append(np.mean(aucs_rf))
        perm_means["xgb"].append(np.mean(aucs_xgb))
        if (pi + 1) % 20 == 0:
            print(f"  perm {pi+1}/{N_PERMS}  ({time.time()-t_perm0:.0f}s elapsed)", flush=True)

    perm_summary = []
    for m in ["en", "rf", "xgb"]:
        pa = np.array(perm_means[m])
        obs = summary[f"mean_auc_{m}"]
        p = float((pa >= obs).sum() + 1) / (N_PERMS + 1)
        perm_summary.append(dict(model=m, observed_mean_auc=obs, observed_sd_auc=summary[f"sd_auc_{m}"],
                                  perm_mean=float(pa.mean()), perm_std=float(pa.std()),
                                  p_value=p, n_perms=N_PERMS))
        print(f"  {m.upper():4s}: obs={obs:.4f}  perm_mean={pa.mean():.4f}+/-{pa.std():.4f}  p={p:.4f}", flush=True)

    pd.DataFrame(perm_summary).to_csv(f"{OUT}/task5_{args.dataset}_permutation_summary.tsv", sep="\t", index=False)
    pd.DataFrame(perm_means).to_csv(f"{OUT}/task5_{args.dataset}_permutation_aucs.tsv", sep="\t", index=False)

    with open(f"{OUT}/task5_{args.dataset}_consensus_params.json", "w") as f:
        json.dump(dict(en=cons_en, rf=cons_rf, xgb=cons_xgb), f, indent=2)

    total_t = time.time() - t0
    print(f"\n[{time.strftime('%H:%M:%S')}] DONE dataset={args.dataset}  total_time={total_t/60:.1f}min", flush=True)


if __name__ == "__main__":
    main()
