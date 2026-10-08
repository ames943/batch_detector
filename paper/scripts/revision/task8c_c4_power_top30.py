#!/usr/bin/env python3
"""
Task 3 (final checks) — does dimensionality explain the power gap between
C4's full genus space (task8_c4_power_v2.tsv) and Cohort 1's top-30 genera
(task8_power_curve.tsv)?

Same method as task8b_c4_power_v2.py, but restricted to C4's OWN top-30 most
prevalent genera (matching the original task8 TOP_P=30 convention), so the
only thing that changes vs. task8_c4_power_v2.py is dimensionality (p=30 vs
p=3005) -- the underlying samples (C4, n=165) are identical.

Settings: n=165 only, delta_R2 in {0.005, 0.0078}, 200 reps.
Output: results/revision/final/task8c_c4_power_top30.tsv
"""
import sys, time
sys.path.insert(0, "scripts/revision")
import numpy as np
import pandas as pd
from lib import clr_transform, E0
from scipy.spatial.distance import pdist, squareform

OUT = "results/revision/final"
SEED = 42
TOP_P = 30
N_SIG_FEAT = 10
TARGETS = [0.005, 0.0078]
N_VALUES = [165]
N_CALIB_REPS = 15
BISECT_ITER = 18
N_POWER_REPS = 200
N_PERMS_INNER = 99


def r2_two_group(D2, y):
    n = D2.shape[0]
    ST = float(np.sum(np.triu(D2, k=1))) / n
    sw = 0.0
    for g in (0, 1):
        m = y == g
        ng = int(m.sum())
        if ng < 2:
            return 0.0
        sw += float(np.sum(np.triu(D2[np.ix_(m, m)], k=1))) / ng
    return (ST - sw) / ST if ST > 0 else 0.0


def perm_p_two_group(D2, y, n_perms, rng):
    n = D2.shape[0]
    ST = float(np.sum(np.triu(D2, k=1))) / n

    def sw_of(yy):
        sw = 0.0
        for g in (0, 1):
            m = yy == g
            ng = int(m.sum())
            if ng < 2:
                return None
            sw += float(np.sum(np.triu(D2[np.ix_(m, m)], k=1))) / ng
        return sw

    sw0 = sw_of(y)
    if not sw0:
        return 1.0
    F0 = (ST - sw0) / sw0 * (n - 2)
    n_ge = 0
    for _ in range(n_perms):
        swp = sw_of(rng.permutation(y))
        Fp = (ST - swp) / swp * (n - 2) if swp else 0.0
        if Fp >= F0:
            n_ge += 1
    return n_ge / n_perms


def main():
    t0 = time.time()
    raw = pd.read_csv("results/ml/lee2022/X_genus_raw.tsv", sep="\t", index_col="run_accession")
    resp = pd.read_csv("metadata/lee2022_labels.tsv", sep="\t").set_index("run_accession").reindex(raw.index)["response"]
    site = pd.read_csv("metadata/lee2022_sites.tsv", sep="\t").set_index("run_accession").reindex(raw.index)["site"]
    keep = resp.notna() & site.notna()
    raw = raw.loc[keep]

    prevalence_full = (raw.values > 0).mean(axis=0)
    top_idx = np.argsort(prevalence_full)[::-1][:TOP_P]
    base_raw = raw.values[:, top_idx].astype(float)
    n_obs, p_top = base_raw.shape
    print(f"[{time.strftime('%H:%M:%S')}] C4 TOP-{TOP_P} base: n_obs={n_obs}, p={p_top} "
          f"(restricted to C4's own most-prevalent genera, same samples as task8_c4_power_v2.py)",
          flush=True)

    prevalence = (base_raw > 0).mean(axis=0)
    sig_cols = np.argsort(prevalence)[::-1][:N_SIG_FEAT]

    def draw(n_target, delta_sig, rng):
        idx = rng.integers(0, n_obs, size=n_target)
        alpha = base_raw[idx] + 1e-6
        G = rng.standard_gamma(alpha)
        comp = G / G.sum(axis=1, keepdims=True)
        X_clr = clr_transform(comp, pseudocount=1e-6)
        y = np.zeros(n_target, dtype=int)
        n_r = n_target // 2
        y[:n_r] = 1
        if delta_sig > 0:
            X_clr[:n_r][:, sig_cols] += delta_sig
        D2 = squareform(pdist(X_clr, metric="euclidean")) ** 2
        return D2, y

    def mean_r2(n_target, delta_sig, n_reps, rng):
        return float(np.mean([r2_two_group(*draw(n_target, delta_sig, rng)) for _ in range(n_reps)]))

    def calibrate(n_target, target_r2, seed):
        if target_r2 <= 0:
            return 0.0
        rng = np.random.default_rng(seed)
        lo, hi = 0.0, 5.0
        while mean_r2(n_target, hi, N_CALIB_REPS, rng) < target_r2:
            hi *= 2.0
            if hi > 200:
                break
        for _ in range(BISECT_ITER):
            mid = (lo + hi) / 2.0
            if mean_r2(n_target, mid, N_CALIB_REPS, rng) < target_r2:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0

    rows = []
    for n_target in N_VALUES:
        e0 = E0(2, n_target)
        for target_delta in TARGETS:
            target_r2 = e0 + target_delta
            s = SEED + int(n_target) + int(target_delta * 1e6)
            delta_sig = calibrate(n_target, target_r2, s)
            rng = np.random.default_rng(s + 1)
            n_sig = 0
            for _ in range(N_POWER_REPS):
                D2, y = draw(n_target, delta_sig, rng)
                if perm_p_two_group(D2, y, N_PERMS_INNER, rng) < 0.05:
                    n_sig += 1
            power = n_sig / N_POWER_REPS
            rows.append(dict(n=n_target, target_delta_R2=target_delta, E0=round(e0, 6),
                              target_R2=round(target_r2, 6), calibrated_delta_sig=round(delta_sig, 5),
                              power=round(power, 4), p_features=p_top))
            print(f"[{time.strftime('%H:%M:%S')}] n={n_target:4d} target_dR2={target_delta:.4f} "
                  f"-> delta_sig={delta_sig:.4f}  power={power:.3f}  ({time.time()-t0:.0f}s elapsed)",
                  flush=True)

    df = pd.DataFrame(rows)
    out_path = f"{OUT}/task8c_c4_power_top30.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nSaved: {out_path}")
    print(f"\n[{time.strftime('%H:%M:%S')}] Total time: {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
