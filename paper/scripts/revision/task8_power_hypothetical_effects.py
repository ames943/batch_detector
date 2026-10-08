#!/usr/bin/env python3
"""
Task 8 — Power over hypothetical effects, using the Dirichlet generator on a
single real cohort (C1, Frankel n=39, top-30 most-prevalent genera for
tractability at n up to 1000).

Part A: power(PERMANOVA, alpha=0.05) for true delta_R2 in {0.0025,0.005,0.01,0.02}
at n in {40,80,165,300,500,1000}. A signal delta is calibrated per (n, target
delta_R2) via binary search (same binary-search-calibration idea as
scripts/phase2_simulation.py, but applied to Dirichlet-resampled real
compositions rather than a fitted Gaussian) so that the *expected observed*
R2 = E0(2,n) + target_delta_R2.

Part B: 95% upper bound on C4's (Lee 2022, n=165) response delta_R2 via
Dirichlet bootstrap of the OBSERVED data (no injected effect — this quantifies
sampling uncertainty in the effect we actually measured), then the n needed
for 80% power at that bound, read off the Part-A curve (nearest calibrated
delta_R2, interpolated over n).
"""
import sys, time
sys.path.insert(0, "scripts/revision")
import numpy as np
import pandas as pd
from lib import clr_transform, E0
from scipy.spatial.distance import pdist, squareform

OUT = "results/revision"
SEED = 42
TOP_P = 30
N_SIG_FEAT = 10
TARGET_DELTAS = [0.0025, 0.005, 0.01, 0.02]
N_VALUES = [40, 80, 165, 300, 500, 1000]
N_CALIB_REPS = 15
BISECT_ITER = 18
N_POWER_REPS = 100
N_PERMS_INNER = 99


def load_base(path_raw):
    raw = pd.read_csv(path_raw, sep="\t", index_col="run_accession")
    prevalence = (raw > 0).mean(axis=0)
    top = prevalence.nlargest(TOP_P).index.tolist()
    return raw[top].values.astype(float)  # (n_obs, TOP_P) raw counts


def r2_two_group(D2, y):
    n = D2.shape[0]
    ST = float(np.sum(np.triu(D2, k=1))) / n
    SW = 0.0
    for g in (0, 1):
        m = y == g
        ng = int(m.sum())
        if ng < 2:
            return 0.0
        SW += float(np.sum(np.triu(D2[np.ix_(m, m)], k=1))) / ng
    return (ST - SW) / ST if ST > 0 else 0.0


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
    if sw0 is None or sw0 == 0:
        return 1.0
    F0 = (ST - sw0) / sw0 * (n - 2)
    n_ge = 0
    for _ in range(n_perms):
        yp = rng.permutation(y)
        swp = sw_of(yp)
        if swp is None or swp == 0:
            Fp = 0.0
        else:
            Fp = (ST - swp) / swp * (n - 2)
        if Fp >= F0:
            n_ge += 1
    return n_ge / n_perms


def draw_dataset(base_raw, n_target, delta_sig, rng):
    """Dirichlet-resample n_target samples from base_raw, balanced R/NR (first
    half R), with an injected CLR-space signal shift of delta_sig on the first
    N_SIG_FEAT features for R-labeled samples."""
    n_obs = base_raw.shape[0]
    idx = rng.integers(0, n_obs, size=n_target)
    alpha = base_raw[idx] + 1e-6
    G = rng.standard_gamma(alpha)
    comp = G / G.sum(axis=1, keepdims=True)
    X_clr = clr_transform(comp, pseudocount=1e-6)
    y = np.zeros(n_target, dtype=int)
    n_r = n_target // 2
    y[:n_r] = 1
    if delta_sig > 0:
        X_clr[:n_r, :N_SIG_FEAT] += delta_sig
    D2 = squareform(pdist(X_clr, metric="euclidean")) ** 2
    return D2, y


def mean_r2_at_delta(base_raw, n_target, delta_sig, n_reps, rng):
    r2s = []
    for _ in range(n_reps):
        D2, y = draw_dataset(base_raw, n_target, delta_sig, rng)
        r2s.append(r2_two_group(D2, y))
    return float(np.mean(r2s))


def calibrate(base_raw, n_target, target_r2, rng_seed):
    if target_r2 <= 0:
        return 0.0
    rng = np.random.default_rng(rng_seed)
    lo, hi = 0.0, 5.0
    while mean_r2_at_delta(base_raw, n_target, hi, N_CALIB_REPS, rng) < target_r2:
        hi *= 2.0
        if hi > 200:
            break
    for _ in range(BISECT_ITER):
        mid = (lo + hi) / 2.0
        r2mid = mean_r2_at_delta(base_raw, n_target, mid, N_CALIB_REPS, rng)
        if r2mid < target_r2:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def power_at(base_raw, n_target, delta_sig, n_reps, n_perms_inner, rng_seed):
    rng = np.random.default_rng(rng_seed)
    n_sig = 0
    for _ in range(n_reps):
        D2, y = draw_dataset(base_raw, n_target, delta_sig, rng)
        p = perm_p_two_group(D2, y, n_perms_inner, rng)
        if p < 0.05:
            n_sig += 1
    return n_sig / n_reps


def main():
    t0 = time.time()
    base_c1 = load_base("results/ml/n118_3cohort/X_genus_raw.tsv")
    # restrict to cohort1 rows for the "single cohort" base
    raw118 = pd.read_csv("results/ml/n118_3cohort/X_genus_raw.tsv", sep="\t", index_col="run_accession")
    mask = raw118.index.str.startswith("SRR5930")
    prevalence = (raw118.loc[mask] > 0).mean(axis=0)
    top = prevalence.nlargest(TOP_P).index.tolist()
    base_c1 = raw118.loc[mask, top].values.astype(float)
    print(f"[{time.strftime('%H:%M:%S')}] base C1: n_obs={base_c1.shape[0]}, p={base_c1.shape[1]}", flush=True)

    rows = []
    for n_target in N_VALUES:
        e0 = E0(2, n_target)
        for target_delta in TARGET_DELTAS:
            target_r2 = e0 + target_delta
            seed_cal = SEED + int(n_target) + int(target_delta * 100000)
            delta_sig = calibrate(base_c1, n_target, target_r2, seed_cal)
            power = power_at(base_c1, n_target, delta_sig, N_POWER_REPS, N_PERMS_INNER, seed_cal + 1)
            rows.append(dict(n=n_target, target_delta_R2=target_delta, E0=round(e0, 6),
                              target_R2=round(target_r2, 6), calibrated_delta_sig=round(delta_sig, 5),
                              power=round(power, 4)))
            print(f"[{time.strftime('%H:%M:%S')}] n={n_target:4d} target_dR2={target_delta:.4f} "
                  f"-> delta_sig={delta_sig:.4f}  power={power:.3f}  "
                  f"({time.time()-t0:.0f}s elapsed)", flush=True)

    dfA = pd.DataFrame(rows)
    dfA.to_csv(f"{OUT}/task8_power_curve.tsv", sep="\t", index=False)
    print("\n" + dfA.to_string(index=False))

    # ── Part B: 95% upper bound on C4's response delta_R2 (Dirichlet bootstrap, no injected effect) ──
    print(f"\n[{time.strftime('%H:%M:%S')}] Part B: C4 Dirichlet bootstrap upper bound ...", flush=True)
    raw_lee = pd.read_csv("results/ml/lee2022/X_genus_raw.tsv", sep="\t", index_col="run_accession")
    resp_lee = pd.read_csv("metadata/lee2022_labels.tsv", sep="\t").set_index("run_accession").reindex(raw_lee.index)["response"]
    keep = resp_lee.notna()
    raw_lee = raw_lee.loc[keep]; resp_lee = resp_lee.loc[keep]
    prevalence_lee = (raw_lee > 0).mean(axis=0)
    top_lee = prevalence_lee.nlargest(TOP_P).index.tolist()
    X_lee_raw = raw_lee[top_lee].values.astype(float)
    y_lee = (resp_lee == "R").astype(int).values
    n_lee = len(X_lee_raw)

    rng = np.random.default_rng(SEED + 555)
    N_BOOT = 500
    boot_r2s = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n_lee, size=n_lee)
        alpha = X_lee_raw[idx] + 1e-6
        G = rng.standard_gamma(alpha)
        comp = G / G.sum(axis=1, keepdims=True)
        X_clr = clr_transform(comp, pseudocount=1e-6)
        y_boot = y_lee[idx]
        if y_boot.sum() in (0, n_lee):
            continue
        D2 = squareform(pdist(X_clr, metric="euclidean")) ** 2
        boot_r2s.append(r2_two_group(D2, y_boot))
    boot_r2s = np.array(boot_r2s)
    e0_lee = E0(2, n_lee)
    boot_delta = boot_r2s - e0_lee
    upper_95 = float(np.percentile(boot_delta, 95))
    print(f"  n={n_lee}, E0={e0_lee:.5f}, observed delta_R2={0.000533:.5f} (from Task 4), "
          f"bootstrap delta_R2 median={np.median(boot_delta):.5f}, 95th pct={upper_95:.5f}", flush=True)

    # n needed for 80% power at this bound: interpolate from the calibrated curve (nearest n, delta)
    # Build an approximate power(n, delta) surface using the 4 calibrated deltas per n, linearly
    # interpolate power vs n at delta closest to the bound; report first n crossing 0.80.
    closest_delta = min(TARGET_DELTAS, key=lambda d: abs(d - upper_95))
    sub = dfA[dfA["target_delta_R2"] == closest_delta].sort_values("n")
    n80 = None
    for i in range(len(sub) - 1):
        p0, p1 = sub.iloc[i]["power"], sub.iloc[i + 1]["power"]
        n0, n1 = sub.iloc[i]["n"], sub.iloc[i + 1]["n"]
        if p0 < 0.80 <= p1:
            frac = (0.80 - p0) / (p1 - p0)
            n80 = n0 + frac * (n1 - n0)
            break
    if n80 is None and (sub["power"] >= 0.80).any():
        n80 = int(sub[sub["power"] >= 0.80]["n"].min())

    partB = dict(n_c4=n_lee, E0_c4=e0_lee, observed_delta_R2_c4=0.000533,
                 dirichlet_boot_median_delta_R2=float(np.median(boot_delta)),
                 dirichlet_boot_95pct_upper_bound_delta_R2=upper_95,
                 closest_calibrated_delta_used_for_lookup=closest_delta,
                 n_for_80pct_power_at_bound=n80, n_bootstrap_reps=len(boot_r2s))
    pd.DataFrame([partB]).to_csv(f"{OUT}/task8_c4_upper_bound.tsv", sep="\t", index=False)
    print("\n" + str(partB))
    print(f"\n[{time.strftime('%H:%M:%S')}] Task 8 complete. Total time: {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
