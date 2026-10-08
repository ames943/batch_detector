#!/usr/bin/env python3
"""
Task A (power-check fix) — recompute power in C4's OWN feature space.

Bug being fixed (see results/revision/SUMMARY.md Section "Power check v2" for
the full writeup): scripts/revision/task8_power_hypothetical_effects.py built
its power curve (Part A) from COHORT 1's top-30 genera, then looked up the
n-for-80%-power for C4's 95% upper-bound effect (delta_R2=0.0078, computed
from C4's own data in Part B) against that C1-space curve. Because 0.0078 is
not one of the 4 calibrated target deltas [0.0025, 0.005, 0.01, 0.02], the
lookup silently substituted the NEAREST one -- 0.01 (|0.01-0.0078|=0.0022 <
|0.005-00.0078|=0.0028) -- interpolated n=80..165 on THAT curve, and reported
the result as if it were "at delta_R2=0.0078". See
scripts/revision/task8_power_hypothetical_effects.py lines 213-224.

This script:
  - uses C4's own genus CLR matrix, with NO prevalence filter -- confirmed by
    reading scripts/revision/task4_c4_alone.py: it calls
    clr_transform(raw.values) directly on all columns of
    results/ml/lee2022/X_genus_raw.tsv, with no prevalence/variance selection
    step. "The same prevalence filtering as Task 4" is therefore no filtering
    at all; the Aitchison distance here is computed on the FULL genus set to
    match Task 4 exactly. (The 10 signal-injection columns are chosen as the
    most-prevalent genera purely so the injected shift lands on a
    non-degenerate column -- this does not filter the distance computation,
    which still uses every genus column, exactly as Task 4 did.)
  - runs the SAME injected-shift binary-search calibration as Task 8's
    power_over_effect_grid / power_curve code (Dirichlet resample + CLR-space
    shift on signal columns for R-labeled samples), so results are
    methodologically comparable to the original Task 8 approach.
  - targets ACTUALLY computed: delta_R2 in {0.0025, 0.005, 0.0078, 0.01}
    (0.0078 now included directly, not approximated).
  - n in {80, 131, 165, 250, 400}, 200 reps/cell, alpha=0.05.

Output: results/revision/task8_c4_power_v2.tsv (does not touch the original
results/revision/task8_power_curve.tsv or task8_c4_upper_bound.tsv).
"""
import sys, time
sys.path.insert(0, "scripts/revision")
import numpy as np
import pandas as pd
from lib import clr_transform, E0
from scipy.spatial.distance import pdist, squareform

OUT = "results/revision"
SEED = 42
N_SIG_FEAT = 10
TARGETS = [0.0025, 0.005, 0.0078, 0.01]
N_VALUES = [80, 131, 165, 250, 400]
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
    n_obs, p_full = raw.shape
    print(f"[{time.strftime('%H:%M:%S')}] C4 base: n_obs={n_obs}, p={p_full} "
          f"(FULL genus set, no prevalence filter -- matches task4_c4_alone.py)", flush=True)

    base_raw = raw.values.astype(float)
    prevalence = (base_raw > 0).mean(axis=0)
    sig_cols = np.argsort(prevalence)[::-1][:N_SIG_FEAT]  # most-prevalent, for injection only

    def draw(n_target, delta_sig, rng):
        idx = rng.integers(0, n_obs, size=n_target)
        alpha = base_raw[idx] + 1e-6
        G = rng.standard_gamma(alpha)
        comp = G / G.sum(axis=1, keepdims=True)
        X_clr = clr_transform(comp, pseudocount=1e-6)   # FULL p_full-dim CLR, matches Task 4
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
                              power=round(power, 4)))
            print(f"[{time.strftime('%H:%M:%S')}] n={n_target:4d} target_dR2={target_delta:.4f} "
                  f"-> delta_sig={delta_sig:.4f}  power={power:.3f}  ({time.time()-t0:.0f}s elapsed)",
                  flush=True)

    df = pd.DataFrame(rows)
    out_path = f"{OUT}/task8_c4_power_v2.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print("\n" + df.to_string(index=False))
    print(f"\nSaved: {out_path}")

    # ── Answers to the specific questions asked ──────────────────────────────
    pw_165_0078 = df[(df["n"] == 165) & (df["target_delta_R2"] == 0.0078)]["power"]
    print("\n--- ANSWERS ---")
    print(f"Power at n=165 for delta_R2=0.0078: "
          f"{pw_165_0078.iloc[0] if len(pw_165_0078) else 'N/A'}")

    sub78 = df[df["target_delta_R2"] == 0.0078].sort_values("n")
    n80 = None
    for i in range(len(sub78) - 1):
        p0, p1 = sub78.iloc[i]["power"], sub78.iloc[i + 1]["power"]
        n0, n1 = sub78.iloc[i]["n"], sub78.iloc[i + 1]["n"]
        if p0 < 0.80 <= p1:
            frac = (0.80 - p0) / (p1 - p0)
            n80 = n0 + frac * (n1 - n0)
            break
    if n80 is None and (sub78["power"] >= 0.80).any():
        n80 = int(sub78[sub78["power"] >= 0.80]["n"].min())
    print(f"n reaching 80% power at delta_R2=0.0078 (from computed grid points; "
          f"linear interpolation only between ACTUALLY COMPUTED n's): {n80}")

    sub165 = df[df["n"] == 165].sort_values("target_delta_R2")
    min_det = None
    for i in range(len(sub165) - 1):
        p0, p1 = sub165.iloc[i]["power"], sub165.iloc[i + 1]["power"]
        d0, d1 = sub165.iloc[i]["target_delta_R2"], sub165.iloc[i + 1]["target_delta_R2"]
        if p0 < 0.80 <= p1:
            frac = (0.80 - p0) / (p1 - p0)
            min_det = d0 + frac * (d1 - d0)
            break
    if min_det is None and (sub165["power"] >= 0.80).any():
        min_det = float(sub165[sub165["power"] >= 0.80]["target_delta_R2"].min())
    print(f"Minimum detectable delta_R2 (80% power) at n=165 "
          f"(interpolated between computed target deltas at n=165): {min_det}")

    print(f"\n[{time.strftime('%H:%M:%S')}] Total time: {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
