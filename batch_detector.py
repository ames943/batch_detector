#!/usr/bin/env python3
"""batch_detector: compare cohort (batch) effects with the response effect
before pooling several cohorts of the same kind of omics data.

For both the cohort labels and the response labels it runs PERMANOVA on
Aitchison distances, subtracts the chance level of R2, and writes a JSON
report, a PDF figure and a power table. The report ends with a GO / CAUTION /
NO-GO pooling recommendation. The cutoffs behind that recommendation are rules
of thumb, not calibrated thresholds (see _pooling_recommendation_v2).

Example:
    python batch_detector.py --input features.tsv --labels meta.tsv \\
        --batch cohort --output out/ --clr

Run with --help for all options. README.md describes the input files and the
output fields.

References
----------
Anderson MJ (2001). Austral Ecology 26(1):32-46.
Aitchison J (1986). The Statistical Analysis of Compositional Data.
McArdle BH, Anderson MJ (2001). Ecology 82(1):290-297.
Phipson B, Smyth GK (2010). Stat Appl Genet Mol Biol 9:Article 39.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import norm as sp_norm

__version__ = "2.0.0"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

_PSEUDOCOUNT = 1e-6


# PERMANOVA

def _ss_total(d2):
    # d2 is the matrix of squared distances
    return float(np.sum(np.triu(d2, k=1))) / d2.shape[0]


def _ss_within(d2, grp):
    sw = 0.0
    for g in np.unique(grp):
        m = grp == g
        ng = int(m.sum())
        if ng < 2:
            continue
        sub = d2[np.ix_(m, m)]
        sw += float(np.sum(np.triu(sub, k=1))) / ng
    return sw


def permanova(D, grouping, n_perms=999, seed=42):
    """One-way PERMANOVA (Anderson 2001) on a distance matrix D.

    Returns a dict with R2, pseudo-F, the permutation p-value and the sums of
    squares.
    """
    grp = np.asarray(grouping)
    n = D.shape[0]
    q = len(np.unique(grp))
    if q < 2:
        raise ValueError("PERMANOVA needs at least 2 distinct groups.")

    d2 = D ** 2
    ST = _ss_total(d2)
    SW = _ss_within(d2, grp)
    SA = ST - SW
    F = (SA / (q - 1)) / (SW / (n - q)) if SW > 0 else 0.0
    R2 = SA / ST if ST > 0 else 0.0

    rng = np.random.default_rng(seed)
    perm_F = np.empty(n_perms)
    for i in range(n_perms):
        g_p = rng.permutation(grp)
        sw_p = _ss_within(d2, g_p)
        sa_p = ST - sw_p
        perm_F[i] = (sa_p / (q - 1)) / (sw_p / (n - q)) if sw_p > 0 else 0.0

    # add-one p-value (Phipson & Smyth 2010)
    p_val = float((perm_F >= F).sum() + 1) / (n_perms + 1)

    return {
        "R2": round(float(R2), 6),
        "F_stat": round(float(F), 4),
        "p_value": round(p_val, 4),
        "SS_total": round(float(ST), 4),
        "SS_between": round(float(SA), 4),
        "SS_within": round(float(SW), 4),
        "n_samples": n,
        "n_groups": q,
        "n_perms": n_perms,
        "perm_F_mean": round(float(perm_F.mean()), 4),
        "perm_F_std": round(float(perm_F.std()), 4),
    }


def chance_expectation_r2(n_groups, n_samples):
    """Expected PERMANOVA R2 when the labels are random: (g - 1) / (n - 1).

    Raw R2 is above zero even with no real effect, and the amount depends on
    the number of groups and on n, so R2 values from different factors or
    sample sizes should be compared after subtracting this.
    """
    return (n_groups - 1) / (n_samples - 1)


# Two-factor (cohort + response) PERMANOVA with marginal effects. This is a
# plain numpy version of vegan's
#   adonis2(dist ~ cohort + response, by = "margin",
#           permutations = how(blocks = cohort))
# built on Gower centering and a hat matrix (McArdle & Anderson 2001).
# Response is permuted within cohorts only. Cohort is permuted freely, because
# it is the blocking factor.

def _dummies(labels):
    levels = sorted(pd.unique(labels))
    Z = np.zeros((len(labels), len(levels) - 1))
    for j, lev in enumerate(levels[1:]):
        Z[:, j] = (labels == lev).astype(float)
    return Z


def _gower_center(D):
    n = D.shape[0]
    A = -0.5 * (D ** 2)
    J = np.eye(n) - np.ones((n, n)) / n
    return J @ A @ J


def _hat_trace(G, X):
    # trace(H G) with H = X (X'X)^+ X', computed as trace((X'X)^+ X'GX) so it
    # stays cheap inside the permutation loops
    M_inv = np.linalg.pinv(X.T @ X)
    XtGX = X.T @ G @ X
    return float(np.trace(M_inv @ XtGX))


def margin_permanova_two_factor(D, cohort, response, n_perms=999, seed=42):
    n = D.shape[0]
    G = _gower_center(D)
    SS_total = float(np.trace(G))
    ones = np.ones((n, 1))
    Z_coh, Z_resp = _dummies(cohort), _dummies(response)

    SS_full = _hat_trace(G, np.hstack([ones, Z_coh, Z_resp]))
    SS_coh_only = _hat_trace(G, np.hstack([ones, Z_coh]))
    SS_resp_only = _hat_trace(G, np.hstack([ones, Z_resp]))

    R2_resp = (SS_full - SS_coh_only) / SS_total
    R2_coh = (SS_full - SS_resp_only) / SS_total

    rng = np.random.default_rng(seed)
    strata = {c: np.where(cohort == c)[0] for c in np.unique(cohort)}

    perm_resp = np.empty(n_perms)
    for i in range(n_perms):
        resp_p = response.copy()
        for c, idxs in strata.items():
            resp_p[idxs] = rng.permutation(response[idxs])
        SS_full_p = _hat_trace(G, np.hstack([ones, Z_coh, _dummies(resp_p)]))
        perm_resp[i] = (SS_full_p - SS_coh_only) / SS_total
    p_resp = float((perm_resp >= R2_resp).sum() + 1) / (n_perms + 1)

    perm_coh = np.empty(n_perms)
    for i in range(n_perms):
        coh_p = rng.permutation(cohort)
        SS_full_p = _hat_trace(G, np.hstack([ones, _dummies(coh_p), Z_resp]))
        perm_coh[i] = (SS_full_p - SS_resp_only) / SS_total
    p_coh = float((perm_coh >= R2_coh).sum() + 1) / (n_perms + 1)

    return dict(R2_response_margin=float(R2_resp), p_response_margin=p_resp,
                R2_cohort_margin=float(R2_coh), p_cohort_margin=p_coh,
                response_perm_restricted_within_cohort=True,
                cohort_perm_restricted=False, n_perms=n_perms)


def rho_permutation_null(D, cohort, response, n_perms=999, seed=42):
    """Null distribution of rho = R2_cohort / R2_response.

    Both labelings are shuffled independently on every draw. The percentile of
    the observed rho within that null shows whether the ratio is bigger than
    what unrelated random labels give.
    """
    d2 = D ** 2
    ST = _ss_total(d2)
    rng = np.random.default_rng(seed)

    def r2(grp):
        sw = _ss_within(d2, grp)
        return (ST - sw) / ST if ST > 0 else 0.0

    obs_r2_resp, obs_r2_coh = r2(response), r2(cohort)
    obs_rho = obs_r2_coh / obs_r2_resp if obs_r2_resp > 0 else np.inf

    null_rho = np.empty(n_perms)
    for i in range(n_perms):
        r2r = r2(rng.permutation(response))
        r2c = r2(rng.permutation(cohort))
        null_rho[i] = r2c / r2r if r2r > 0 else np.nan
    valid = null_rho[np.isfinite(null_rho)]
    pctile = float((valid <= obs_rho).mean() * 100) if len(valid) else float("nan")

    return dict(observed_rho=float(obs_rho) if np.isfinite(obs_rho) else None,
                null_rho_mean=float(np.nanmean(null_rho)),
                null_rho_median=float(np.nanmedian(null_rho)),
                percentile_of_observed=pctile, n_perms=n_perms)


# CLR and Aitchison distance

def clr_transform(X, pseudocount=_PSEUDOCOUNT):
    """Centered log-ratio transform of each row, after adding a pseudocount."""
    lx = np.log(X + pseudocount)
    return lx - lx.mean(axis=1, keepdims=True)


def aitchison_dist(X_clr):
    return squareform(pdist(X_clr, metric="euclidean"))


def clr_multiplicative_replacement(X_raw, delta=1e-5):
    """Multiplicative replacement of zeros (Martin-Fernandez et al. 2003),
    then CLR."""
    comp = X_raw / X_raw.sum(axis=1, keepdims=True)
    n = comp.shape[0]
    zero_mask = comp == 0
    n_zeros = zero_mask.sum(axis=1)
    out = comp.copy()
    for i in range(n):
        if n_zeros[i] == 0:
            continue
        out[i, zero_mask[i]] = delta
        nonzero = ~zero_mask[i]
        out[i, nonzero] = comp[i, nonzero] * (1 - n_zeros[i] * delta)
    lx = np.log(out)
    return lx - lx.mean(axis=1, keepdims=True)


def clr_transform_with_method(X_raw, zero_handling, pseudocount):
    if zero_handling == "multiplicative_replacement":
        return clr_multiplicative_replacement(X_raw, delta=pseudocount)
    return clr_transform(X_raw, pseudocount=pseudocount)


# Dirichlet bootstrap and power

def _dirichlet_draw(alpha, rng):
    # a Dirichlet draw per row, via normalized gamma variates
    G = rng.standard_gamma(alpha)
    return G / G.sum(axis=1, keepdims=True)


def _wilson_ci(k, n, alpha=0.05):
    if n == 0:
        return 0.0, 1.0
    z = sp_norm.ppf(1 - alpha / 2)
    p = k / n
    c = (p + z**2 / (2 * n)) / (1 + z**2 / n)
    m = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / (1 + z**2 / n)
    return max(0.0, c - m), min(1.0, c + m)


def _permanova_pval_fast(d2, y, n_perms_inner, rng):
    """Permutation p-value for a two-group split, used inside the power loops."""
    n = d2.shape[0]
    ST = float(np.sum(np.triu(d2, k=1))) / n
    SW = 0.0
    for g in (0, 1):
        m = y == g
        ng = int(m.sum())
        if ng < 2:
            return 1.0
        sub = d2[np.ix_(m, m)]
        SW += float(np.sum(np.triu(sub, k=1))) / ng

    if SW == 0:
        return 1.0
    F0 = (ST - SW) / SW * (n - 2)

    n_ge = 0
    for _ in range(n_perms_inner):
        y_p = rng.permutation(y)
        sw_p = 0.0
        for g in (0, 1):
            m = y_p == g
            ng = int(m.sum())
            if ng < 2:
                continue
            sub = d2[np.ix_(m, m)]
            sw_p += float(np.sum(np.triu(sub, k=1))) / ng
        F_p = (ST - sw_p) / sw_p * (n - 2) if sw_p > 0 else 0.0
        if F_p >= F0:
            n_ge += 1

    # no +1 correction in the inner loop
    return n_ge / n_perms_inner


def dirichlet_power_analysis(X_raw, y_binary, n_values, n_reps=200,
                             n_perms_inner=199, pseudocount=_PSEUDOCOUNT,
                             seed=42):
    """Legacy (v1) power curve for the response effect actually observed.

    For each target n, resample rows with replacement, draw a fresh Dirichlet
    composition for every resampled row (so there are no zero-distance
    duplicates), run PERMANOVA, and count how often p < 0.05. Returns one dict
    per n.
    """
    rng = np.random.default_rng(seed)
    alpha_mat = X_raw + pseudocount
    n_obs = len(y_binary)
    y = np.asarray(y_binary, dtype=int)

    rows = []
    for n_target in n_values:
        n_sig = 0
        for _ in range(n_reps):
            idx = rng.integers(0, n_obs, size=n_target)
            alpha_draw = alpha_mat[idx]
            X_sim = _dirichlet_draw(alpha_draw, rng)
            X_clr = clr_transform(X_sim, pseudocount)
            D_sim = squareform(pdist(X_clr, metric="euclidean"))
            d2 = D_sim ** 2
            y_sim = y[idx]

            if y_sim.sum() in (0, n_target):
                continue

            p = _permanova_pval_fast(d2, y_sim, n_perms_inner, rng)
            if p < 0.05:
                n_sig += 1

        power = n_sig / n_reps
        ci_lo, ci_hi = _wilson_ci(n_sig, n_reps)
        rows.append({
            "n": n_target,
            "power": round(power, 4),
            "ci_lower": round(ci_lo, 4),
            "ci_upper": round(ci_hi, 4),
            "n_reps": n_reps,
            "n_perms_inner": n_perms_inner,
        })
        log.info("  n=%4d  power=%.3f  95%% CI [%.3f, %.3f]",
                 n_target, power, ci_lo, ci_hi)

    return rows


def power_over_effect_grid(X_raw, n_values, delta_grid, n_reps=100,
                           n_perms_inner=99, n_calib_reps=15,
                           calib_bisect_iter=18, top_p=30, n_sig_feat=10,
                           seed=42):
    """Default (v2) power analysis over hypothetical effect sizes.

    The v1 curve resamples the effect that was observed, which is circular
    when that effect is close to zero. Here we ask instead how large a
    response effect (excess R2 over chance) has to be, and at what n, to be
    detected with 80% power.

    The top_p most prevalent features are used as the Dirichlet base. For
    every (n, target delta R2) a CLR shift on n_sig_feat features is found by
    bisection so the expected R2 is E0 + delta, then power is the share of
    simulated datasets with permutation p < 0.05.
    """
    prevalence = (X_raw > 0).mean(axis=0)
    top_idx = np.argsort(prevalence)[::-1][:top_p]
    base = X_raw[:, top_idx]
    n_obs = base.shape[0]

    def draw(n_target, delta_sig, rng):
        idx = rng.integers(0, n_obs, size=n_target)
        alpha = base[idx] + 1e-6
        G = rng.standard_gamma(alpha)
        comp = G / G.sum(axis=1, keepdims=True)
        X_clr = clr_transform(comp, 1e-6)
        y = np.zeros(n_target, dtype=int)
        n_r = n_target // 2
        y[:n_r] = 1
        if delta_sig > 0:
            X_clr[:n_r, :n_sig_feat] += delta_sig
        d2 = squareform(pdist(X_clr, metric="euclidean")) ** 2
        return d2, y

    def r2_of(d2, y):
        n = d2.shape[0]
        ST = float(np.sum(np.triu(d2, k=1))) / n
        sw = 0.0
        for g in (0, 1):
            m = y == g
            ng = int(m.sum())
            if ng < 2:
                return 0.0
            sw += float(np.sum(np.triu(d2[np.ix_(m, m)], k=1))) / ng
        return (ST - sw) / ST if ST > 0 else 0.0

    def perm_p(d2, y, n_perms, rng):
        n = d2.shape[0]
        ST = float(np.sum(np.triu(d2, k=1))) / n

        def sw_of(yy):
            sw = 0.0
            for g in (0, 1):
                m = yy == g
                ng = int(m.sum())
                if ng < 2:
                    return None
                sw += float(np.sum(np.triu(d2[np.ix_(m, m)], k=1))) / ng
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

    def calibrate(n_target, target_r2, rng_seed):
        if target_r2 <= 0:
            return 0.0
        rng = np.random.default_rng(rng_seed)
        lo, hi = 0.0, 5.0
        while np.mean([r2_of(*draw(n_target, hi, rng)) for _ in range(n_calib_reps)]) < target_r2:
            hi *= 2.0
            if hi > 200:
                break
        for _ in range(calib_bisect_iter):
            mid = (lo + hi) / 2.0
            r2mid = np.mean([r2_of(*draw(n_target, mid, rng)) for _ in range(n_calib_reps)])
            if r2mid < target_r2:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0

    rows = []
    for n_t in n_values:
        e0 = chance_expectation_r2(2, n_t)
        for target_delta in delta_grid:
            target_r2 = e0 + target_delta
            s = seed + int(n_t) + int(target_delta * 1e6)
            delta_sig = calibrate(n_t, target_r2, s)
            rng = np.random.default_rng(s + 1)
            n_sig = 0
            for _ in range(n_reps):
                d2, y = draw(n_t, delta_sig, rng)
                if perm_p(d2, y, n_perms_inner, rng) < 0.05:
                    n_sig += 1
            power = n_sig / n_reps
            rows.append(dict(n=n_t, target_delta_R2=target_delta, E0=round(e0, 6),
                             calibrated_delta_sig=round(delta_sig, 5), power=round(power, 4)))
            log.info("  n=%4d  target_dR2=%.4f  power=%.3f", n_t, target_delta, power)
    return rows


def _find_n_80pct(power_rows):
    """Smallest n whose power is at least 0.80, or None."""
    for r in power_rows:
        if r["power"] >= 0.80:
            return r["n"]
    return None


def simulate_at_current_n(X_raw, y_binary, n_target, n_sim_reps=50,
                          n_perms_inner=99, pseudocount=_PSEUDOCOUNT, seed=99):
    """Short independent rerun of the v1 power estimate at the current n."""
    rows = dirichlet_power_analysis(
        X_raw, y_binary,
        n_values=[n_target],
        n_reps=n_sim_reps,
        n_perms_inner=n_perms_inner,
        pseudocount=pseudocount,
        seed=seed,
    )
    return rows[0] if rows else {}


# Pooling recommendation

def _pooling_recommendation_v2(response_delta_r2, response_p, rho,
                               rho_percentile, single_cohort=False):
    """Return (GO | CAUTION | NO-GO, reason).

    The rule works on the chance-corrected response effect and on rho =
    R2_cohort / R2_response:

      NO-GO    response delta R2 <= 0, or response p >= 0.05, or rho > 4, or
               rho sits at or above the 95th percentile of its shuffle null
      CAUTION  significant response effect and 0.25 < rho <= 4
      GO       significant response effect and rho <= 0.25 (or one cohort only)

    The 0.25 and 4 cutoffs are conservative rules of thumb. They are not
    derived from the paper's simulation: the calibrated grid showed no clean
    rho-based threshold, and batch correction only helped at the highest
    simulated signal level, whatever the ratio was.
    """
    if response_delta_r2 is None or response_p is None:
        return "CAUTION", "Response PERMANOVA unavailable — cannot assess chance-calibrated signal."

    if response_delta_r2 <= 0:
        return "NO-GO", (
            f"Response R2 does not exceed its own chance expectation E0 "
            f"(delta_R2 = {response_delta_r2:+.5f} <= 0). There is no detectable "
            "excess biological signal at this n, independent of any batch effect — "
            "pooling cannot help what is not there to begin with."
        )
    if response_p >= 0.05:
        return "NO-GO", (
            f"Response delta_R2 is positive ({response_delta_r2:+.5f}) but not "
            f"statistically significant (p = {response_p:.3f}). Cannot distinguish "
            "the excess signal from sampling noise at this n."
        )
    if single_cohort:
        return "GO", (
            f"Single-cohort mode. Significant chance-calibrated excess signal "
            f"(delta_R2 = {response_delta_r2:+.5f}, p = {response_p:.3f}) — no batch "
            "factor to weigh it against."
        )
    if rho is None:
        return "CAUTION", (
            f"Significant chance-calibrated signal (delta_R2={response_delta_r2:+.5f}, "
            f"p={response_p:.3f}) but no batch factor was supplied to assess rho."
        )
    if rho_percentile is not None and rho_percentile >= 95:
        return "NO-GO", (
            f"rho = batch_R2/response_R2 = {rho:.2f}x sits at/above the "
            f"{rho_percentile:.1f}th percentile of its own independent-shuffle null — "
            "the batch effect is more extreme than pure chance from unrelated random "
            "labels would produce. This is a structural batch dominance, not noise."
        )
    if rho <= 0.25:
        return "GO", (
            f"Significant chance-calibrated signal (delta_R2={response_delta_r2:+.5f}, "
            f"p={response_p:.3f}) with rho={rho:.2f}x <= the 0.25x heuristic cutoff. "
            "Batch effect is small relative to the biological signal; "
            "pooling appears appropriate."
        )
    if rho <= 4.0:
        return "CAUTION", (
            f"Significant chance-calibrated signal (delta_R2={response_delta_r2:+.5f}, "
            f"p={response_p:.3f}) but rho={rho:.2f}x is between the 0.25x and 4x "
            "heuristic cutoffs. Apply and independently validate batch correction; "
            "re-run this tool on the corrected data to confirm delta_R2 has not "
            "collapsed."
        )
    return "NO-GO", (
        f"rho={rho:.2f}x exceeds the 4x heuristic cutoff. "
        "Batch variance dominates the (real, significant) biological signal."
    )


def _pooling_recommendation(batch_signal_ratio, response_p):
    """Original (v1) rule, kept for --legacy-recommendation.

    Uses the raw batch:signal R2 ratio with cutoffs 3 and 8, without
    subtracting chance levels.
    """
    if batch_signal_ratio is None:
        if response_p < 0.05:
            return "GO", (
                f"Single-cohort mode. Significant biological signal detected "
                f"(PERMANOVA p = {response_p:.3f})."
            )
        return "CAUTION", (
            f"Single-cohort mode. Biological signal not significant "
            f"(PERMANOVA p = {response_p:.3f}). Check the power curve to "
            "determine the same-protocol n needed for 80% power."
        )

    if batch_signal_ratio > 8:
        return "NO-GO", (
            f"Batch variance dominates biological signal by {batch_signal_ratio:.1f}×. "
            "Naive pooling across these cohorts is not recommended. "
            "Seek a same-protocol cohort or a substantially larger n "
            "(see power curve); validate any batch-correction method independently."
        )
    if batch_signal_ratio >= 3:
        return "CAUTION", (
            f"Moderate batch:signal ratio ({batch_signal_ratio:.1f}×). "
            "Apply and independently validate batch correction before downstream "
            "analysis. Re-run this tool on the corrected data to confirm batch R² "
            "has decreased without collapsing response R²."
        )
    if response_p < 0.05:
        return "GO", (
            f"Low batch:signal ratio ({batch_signal_ratio:.1f}×) with significant "
            f"biological signal (p = {response_p:.3f}). Pooling appears appropriate."
        )
    return "CAUTION", (
        f"Low batch:signal ratio ({batch_signal_ratio:.1f}×) but biological signal "
        f"is not significant (p = {response_p:.3f}). Dataset may be underpowered "
        "at the current n — check the power curve."
    )


# Figure

_PALETTE = {
    "batch": "#E53935",
    "response": "#1E88E5",
    "GO": "#43A047",
    "CAUTION": "#FB8C00",
    "NO-GO": "#E53935",
    "power_line": "#1565C0",
    "power_fill": "#90CAF9",
    "threshold_80": "#C62828",
    "current_n": "#E65100",
    "n_80pct": "#AD1457",
}


def _sig_stars(p):
    if p < 0.001:
        return "p < 0.001 (***)"
    if p < 0.01:
        return f"p = {p:.3f} (**)"
    if p < 0.05:
        return f"p = {p:.3f} (*)"
    if p < 0.10:
        return f"p = {p:.3f} (.)"
    return f"p = {p:.3f} (ns)"


def generate_figure(res_response, res_batch, power_rows, n_current, n_80pct,
                    recommendation, out_path):
    """Two panels: PERMANOVA R2 bars (left), v1 power curve (right).

    The power panel is only filled when power_rows is non-empty, which is the
    case with --legacy-power.
    """
    rec_color = _PALETTE.get(recommendation, "#555555")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5))
    fig.subplots_adjust(wspace=0.38)

    # left panel: R2 bars
    ax = axes[0]
    bar_labels, bar_heights, bar_colors, bar_pvals = [], [], [], []

    if res_response:
        bar_labels.append("Response\n(biological signal)")
        bar_heights.append(res_response["R2"])
        bar_colors.append(_PALETTE["response"])
        bar_pvals.append(res_response["p_value"])

    if res_batch:
        bar_labels.append("Batch\n(technical noise)")
        bar_heights.append(res_batch["R2"])
        bar_colors.append(_PALETTE["batch"])
        bar_pvals.append(res_batch["p_value"])

    x_pos = np.arange(len(bar_labels))
    bars = ax.bar(
        x_pos, bar_heights,
        color=bar_colors, width=0.45,
        edgecolor="white", linewidth=1.4, zorder=3,
    )

    y_max = max(bar_heights) if bar_heights else 0.10
    for bar, pval, h in zip(bars, bar_pvals, bar_heights):
        stars = "***" if pval < 0.001 else "**" if pval < 0.01 \
            else "*" if pval < 0.05 else "ns"
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h + y_max * 0.03,
            stars,
            ha="center", va="bottom",
            fontsize=13, fontweight="bold",
            color="black",
        )
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h / 2,
            f"R² = {h:.4f}",
            ha="center", va="center",
            fontsize=9, color="white", fontweight="bold",
        )

    if res_batch and res_response and res_response["R2"] > 0:
        ratio = res_batch["R2"] / res_response["R2"]
        ax.text(
            0.97, 0.97,
            f"Batch : Signal = {ratio:.1f}×",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10, fontweight="bold", color=rec_color,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor=rec_color, linewidth=1.5, alpha=0.90),
        )

    ax.set_xticks(x_pos)
    ax.set_xticklabels(bar_labels, fontsize=11)
    ax.set_ylabel("PERMANOVA R²  (Aitchison distance)", fontsize=11)
    ax.set_ylim(0, y_max * 1.30)
    ax.set_title("(A)  Batch vs. Biological Signal", fontsize=12, fontweight="bold", pad=10)
    ax.grid(True, axis="y", alpha=0.28, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    if res_response or res_batch:
        legend_lines = []
        if res_response:
            legend_lines.append(f"Response: {_sig_stars(res_response['p_value'])}")
        if res_batch:
            legend_lines.append(f"Batch: {_sig_stars(res_batch['p_value'])}")
        ax.text(
            0.03, 0.97, "\n".join(legend_lines),
            transform=ax.transAxes, ha="left", va="top",
            fontsize=8.5, color="#333333",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#F5F5F5",
                      edgecolor="#BDBDBD", alpha=0.90),
        )

    # right panel: power curve
    ax2 = axes[1]

    if power_rows:
        ns_pw = [r["n"] for r in power_rows]
        ps_pw = [r["power"] for r in power_rows]
        ci_lo = [r["ci_lower"] for r in power_rows]
        ci_hi = [r["ci_upper"] for r in power_rows]

        ax2.fill_between(ns_pw, ci_lo, ci_hi,
                         color=_PALETTE["power_fill"], alpha=0.40, zorder=2,
                         label="95% Wilson CI")
        ax2.plot(ns_pw, ps_pw, "o-",
                 color=_PALETTE["power_line"], linewidth=2.0, markersize=6,
                 zorder=4, label="Estimated power")

        ax2.axhline(0.80, color=_PALETTE["threshold_80"], linestyle="--",
                    linewidth=1.4, zorder=3, label="80% power target")

        pw_at_current = next((r["power"] for r in power_rows if r["n"] == n_current), None)
        ax2.axvline(n_current, color=_PALETTE["current_n"], linestyle=":",
                    linewidth=1.8, zorder=3, label=f"Current n = {n_current}")
        if pw_at_current is not None:
            ax2.scatter([n_current], [pw_at_current],
                        color=_PALETTE["current_n"], s=80, zorder=5)

        if n_80pct is not None:
            ax2.axvline(n_80pct, color=_PALETTE["n_80pct"], linestyle="-.",
                        linewidth=1.3, zorder=3, alpha=0.85,
                        label=f"n for 80% power ≈ {n_80pct}")

        ax2.set_xlabel("Same-protocol sample size  n", fontsize=11)
        ax2.set_ylabel("Empirical power  P(PERMANOVA p < 0.05)", fontsize=11)
        ax2.set_ylim(-0.04, 1.08)
        ax2.set_xlim(min(ns_pw) * 0.92, max(ns_pw) * 1.05)
        ax2.legend(fontsize=8.5, loc="lower right", framealpha=0.92)
        ax2.grid(True, alpha=0.25, zorder=0)
        ax2.spines[["top", "right"]].set_visible(False)
    else:
        ax2.text(0.5, 0.5,
                 "Power curve drawn only with --legacy-power.\n"
                 "Default power grid: batch_detector_power_grid.tsv",
                 ha="center", va="center", transform=ax2.transAxes,
                 fontsize=11, color="#777777")

    ax2.set_title("(B)  Dirichlet-Corrected Power Curve", fontsize=12,
                  fontweight="bold", pad=10)

    fig.suptitle(
        f"BatchDetector  |  Recommendation:  {recommendation}  "
        f"(n = {n_current})",
        fontsize=12, fontweight="bold", color=rec_color, y=1.02,
    )

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Figure saved: %s", out_path)


# Input handling

def _load_feature_matrix(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    log.info("Feature matrix: %d samples x %d features (%s)",
             df.shape[0], df.shape[1], path)
    return df


def _load_labels(path, response_col):
    df = pd.read_csv(path, sep="\t", index_col=0)
    if response_col not in df.columns:
        raise ValueError(
            f"Response column '{response_col}' not found in labels file. "
            f"Available columns: {list(df.columns)}"
        )
    return df[response_col]


def _load_batch(batch_arg, labels_path):
    """--batch is either a column of the labels file or a 2-column TSV."""
    labels_df = pd.read_csv(labels_path, sep="\t", index_col=0)
    if batch_arg in labels_df.columns:
        log.info("Batch taken from column '%s' of the labels file.", batch_arg)
        return labels_df[batch_arg]
    try:
        batch_df = pd.read_csv(batch_arg, sep="\t", index_col=0)
        col = batch_df.columns[0]
        log.info("Batch taken from file '%s' (column '%s').", batch_arg, col)
        return batch_df[col]
    except (FileNotFoundError, pd.errors.ParserError, IsADirectoryError) as exc:
        raise ValueError(
            f"--batch '{batch_arg}' is neither a column in the labels file "
            f"nor a readable TSV file. Error: {exc}"
        ) from exc


def _preprocess(X, response, batch, apply_clr, pseudocount=_PSEUDOCOUNT,
                zero_handling="pseudocount"):
    """Align samples across files, fill NaN with 0, drop constant features and
    optionally apply CLR.

    Returns (X_out, X_raw, response_values, batch_values, sample_ids).
    X_raw is kept because the Dirichlet power simulation needs untransformed
    values.
    """
    common = X.index
    if response is not None:
        common = common.intersection(response.dropna().index)
    if batch is not None:
        common = common.intersection(batch.dropna().index)

    n_dropped = len(X) - len(common)
    if n_dropped > 0:
        log.warning("Dropped %d samples with missing labels after alignment.",
                    n_dropped)
    if len(common) == 0:
        raise ValueError(
            "No samples remain after aligning feature matrix with labels. "
            "Check that sample IDs match across files."
        )

    X = X.loc[common]
    response = response.loc[common] if response is not None else None
    batch = batch.loc[common] if batch is not None else None

    n_nan = int(X.isna().sum().sum())
    if n_nan > 0:
        log.warning("Imputing %d NaN feature values with 0.", n_nan)
        X = X.fillna(0)

    feat_var = X.var(axis=0)
    n_zero_var = int((feat_var == 0).sum())
    if n_zero_var > 0:
        log.info("Dropping %d zero-variance features.", n_zero_var)
        X = X.loc[:, feat_var > 0]
    if X.shape[1] == 0:
        raise ValueError("All features have zero variance. Cannot proceed.")

    X_raw = X.values.astype(float)

    if np.any(X_raw < 0):
        log.warning(
            "Feature matrix contains negative values. The Dirichlet power "
            "analysis needs non-negative counts or relative abundances and "
            "will be skipped."
        )

    if apply_clr:
        X_out = clr_transform_with_method(X_raw, zero_handling, pseudocount)
        log.info("CLR applied (zero_handling=%s, pseudocount/delta=%.0e).",
                 zero_handling, pseudocount)
    else:
        X_out = X_raw.copy()
        log.info("No CLR; using the values as given.")

    log.info("Final dataset: n = %d samples, p = %d features.", *X_out.shape)

    return (
        X_out,
        X_raw,
        response.values if response is not None else None,
        batch.values if batch is not None else None,
        list(common),
    )


def _build_parser():
    p = argparse.ArgumentParser(
        prog="batch_detector.py",
        description=(
            "Compare cohort (batch) effects with the response effect before "
            "pooling cohorts.\n\n"
            "Runs PERMANOVA on Aitchison distances for response and for cohort, "
            "subtracts the chance level of R2, and writes a JSON report, a PDF "
            "figure and a power table with a GO / CAUTION / NO-GO pooling "
            "recommendation. The recommendation cutoffs are heuristic.\n\n"
            "Inputs are tab-separated files: a samples x features table and a "
            "labels table, both with sample IDs in the first column. See "
            "README.md for details."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  # three cohorts, raw counts, CLR applied by the tool
  python batch_detector.py --input features.tsv --labels meta.tsv \\
      --batch cohort --output out/ --clr

  # already-transformed data (no --clr), cohort labels in a separate file
  python batch_detector.py --input features_clr.tsv --labels meta.tsv \\
      --batch cohort_file.tsv --output out/

  # one cohort only (no --batch)
  python batch_detector.py --input features.tsv --labels meta.tsv \\
      --output out/ --clr

  # old (v1) power curve plus a 50-rep check at the current n
  python batch_detector.py --input features.tsv --labels meta.tsv \\
      --batch cohort --output out/ --clr --legacy-power --simulate
        """,
    )

    p.add_argument("--input", required=True,
                   help="TSV, samples x features, first column = sample ID.")
    p.add_argument("--labels", required=True,
                   help="TSV, first column = sample ID; must contain --response-col.")
    p.add_argument("--output", required=True,
                   help="Output directory (created if absent).")

    p.add_argument("--batch", default=None,
                   help="Cohort/batch labels: a column name in --labels, or the "
                        "path to a 2-column TSV (sample_id, batch).")
    p.add_argument("--response-col", default="response",
                   help="Column holding the binary outcome (default: response).")
    p.add_argument("--clr", action="store_true",
                   help="Apply the CLR transform (use for counts or relative "
                        "abundances).")
    p.add_argument("--pseudocount", type=float, default=1e-6,
                   help="Pseudocount, or delta for multiplicative replacement, "
                        "used for CLR and the Dirichlet draws (default: 1e-6).")
    p.add_argument("--zero-handling", choices=["pseudocount", "multiplicative_replacement"],
                   default="pseudocount",
                   help="How zeros are replaced before CLR (default: pseudocount).")

    p.add_argument("--n-perms", type=int, default=999,
                   help="PERMANOVA permutations (default: 999).")

    p.add_argument("--legacy-power", action="store_true",
                   help="Use the v1 power analysis (resamples the observed "
                        "effect at several n) instead of the v2 default "
                        "(hypothetical chance-corrected effect sizes).")
    p.add_argument("--delta-grid", type=str, default="0.0025,0.005,0.01,0.02",
                   help="v2 power: comma-separated target delta R2 values "
                        "(default: 0.0025,0.005,0.01,0.02).")
    p.add_argument("--n-grid", type=str, default=None,
                   help="v2 power: comma-separated n values "
                        "(default: current n times 0.5, 1, 2 and 4).")
    p.add_argument("--power-reps", type=int, default=100,
                   help="v2 power: simulated datasets per (n, delta R2) cell "
                        "(default: 100).")
    p.add_argument("--power-perms-inner", type=int, default=99,
                   help="v2 power: permutations per simulated dataset "
                        "(default: 99).")
    p.add_argument("--calib-reps", type=int, default=15,
                   help="v2 power: simulated datasets per calibration step "
                        "(default: 15).")
    p.add_argument("--calib-iters", type=int, default=18,
                   help="v2 power: bisection steps per calibration "
                        "(default: 18).")
    p.add_argument("--n-boot", type=int, default=200,
                   help="v1 power: bootstrap reps per n (default: 200).")
    p.add_argument("--n-perms-inner", type=int, default=199,
                   help="v1 power: permutations per rep (default: 199).")

    p.add_argument("--legacy-recommendation", action="store_true",
                   help="Use the v1 GO/CAUTION/NO-GO rule (raw batch:signal "
                        "ratio, cutoffs 3 and 8) instead of the v2 default "
                        "(chance-corrected, heuristic rho cutoffs 0.25 and 4).")

    p.add_argument("--simulate", action="store_true",
                   help="Check the v1 power estimate with a 50-rep Dirichlet "
                        "bootstrap at the current n. Only used together with "
                        "--legacy-power.")

    p.add_argument("--seed", type=int, default=42,
                   help="Random seed (default: 42).")
    p.add_argument("--version", action="version",
                   version=f"batch_detector {__version__}")

    return p


def main(argv=None):
    args = _build_parser().parse_args(argv)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("batch_detector %s", __version__)
    log.info("Input       : %s", args.input)
    log.info("Labels      : %s  (column: '%s')", args.labels, args.response_col)
    log.info("Batch       : %s", args.batch or "(none, single-cohort mode)")
    log.info("CLR         : %s  (pseudocount = %.0e)", args.clr, args.pseudocount)
    log.info("PERMANOVA   : n_perms = %d", args.n_perms)
    log.info("Output dir  : %s", out_dir)

    X_df = _load_feature_matrix(args.input)
    response_series = _load_labels(args.labels, args.response_col)
    batch_series = _load_batch(args.batch, args.labels) if args.batch else None

    X_out, X_raw, response_arr, batch_arr, sample_ids = _preprocess(
        X_df, response_series, batch_series,
        apply_clr=args.clr, pseudocount=args.pseudocount,
        zero_handling=args.zero_handling,
    )
    n = len(sample_ids)

    if n < 4:
        log.error("Need at least 4 samples for PERMANOVA. Got n = %d.", n)
        sys.exit(1)

    log.info("Computing the %d x %d Aitchison distance matrix", n, n)
    D = aitchison_dist(X_out)
    log.info("  max = %.2f   mean (off-diagonal) = %.2f",
             D.max(), D[D > 0].mean())

    # response PERMANOVA
    res_response = None
    unique_resp = np.unique(response_arr) if response_arr is not None else []
    if response_arr is not None and len(unique_resp) >= 2:
        log.info("PERMANOVA, response (n_perms = %d)", args.n_perms)
        res_response = permanova(D, response_arr, n_perms=args.n_perms, seed=args.seed)
        log.info("  R2 = %.4f   F = %.4f   %s",
                 res_response["R2"], res_response["F_stat"],
                 _sig_stars(res_response["p_value"]))
    elif response_arr is not None:
        log.warning("Response PERMANOVA skipped: fewer than 2 unique labels "
                    "(%s). Check --response-col.", unique_resp)

    # batch PERMANOVA
    res_batch = None
    if batch_arr is not None:
        n_batches = len(np.unique(batch_arr))
        if n_batches < 2:
            log.warning("Only 1 unique batch value, treating as single-cohort mode.")
            batch_arr = None
        else:
            log.info("PERMANOVA, batch (%d groups, n_perms = %d)",
                     n_batches, args.n_perms)
            res_batch = permanova(D, batch_arr, n_perms=args.n_perms,
                                  seed=args.seed + 1)
            log.info("  R2 = %.4f   F = %.4f   %s",
                     res_batch["R2"], res_batch["F_stat"],
                     _sig_stars(res_batch["p_value"]))

    batch_signal_ratio = None
    if res_batch and res_response and res_response["R2"] > 0:
        batch_signal_ratio = round(res_batch["R2"] / res_response["R2"], 2)
        log.info("Batch:signal ratio = %.2fx", batch_signal_ratio)

    # chance-corrected excess R2
    response_E0 = response_delta_R2 = None
    if res_response:
        response_E0 = chance_expectation_r2(res_response["n_groups"], n)
        response_delta_R2 = res_response["R2"] - response_E0
        log.info("Response chance level E0 = %.5f, delta_R2 = %+.5f",
                 response_E0, response_delta_R2)
    batch_E0 = batch_delta_R2 = None
    if res_batch:
        batch_E0 = chance_expectation_r2(res_batch["n_groups"], n)
        batch_delta_R2 = res_batch["R2"] - batch_E0
        log.info("Batch chance level    E0 = %.5f, delta_R2 = %+.5f",
                 batch_E0, batch_delta_R2)

    # joint model and rho null, only when there is a batch factor
    joint_model = None
    rho_result = None
    if res_batch and response_arr is not None:
        log.info("Joint (marginal) model, response permuted within cohort")
        joint_model = margin_permanova_two_factor(D, batch_arr, response_arr,
                                                  n_perms=args.n_perms, seed=args.seed)
        log.info("  R2_response|cohort = %.5f  p=%.4f   R2_cohort|response = %.5f  p=%.4f",
                 joint_model["R2_response_margin"], joint_model["p_response_margin"],
                 joint_model["R2_cohort_margin"], joint_model["p_cohort_margin"])

        log.info("Permutation null for rho = batch_R2 / response_R2")
        rho_result = rho_permutation_null(D, batch_arr, response_arr,
                                          n_perms=args.n_perms, seed=args.seed + 3)
        log.info("  observed rho = %s   null mean = %.2f   percentile of observed = %.1f",
                 f"{rho_result['observed_rho']:.2f}x" if rho_result["observed_rho"] is not None else "inf",
                 rho_result["null_rho_mean"], rho_result["percentile_of_observed"])

    resp_p = res_response["p_value"] if res_response else 1.0
    if args.legacy_recommendation:
        rec, rec_reason = _pooling_recommendation(batch_signal_ratio, resp_p)
    else:
        rec, rec_reason = _pooling_recommendation_v2(
            response_delta_r2=response_delta_R2,
            response_p=(res_response["p_value"] if res_response else None),
            rho=(rho_result["observed_rho"] if rho_result else None),
            rho_percentile=(rho_result["percentile_of_observed"] if rho_result else None),
            single_cohort=(res_batch is None),
        )
    log.info("Recommendation: %s", rec)
    log.info("Reason: %s", rec_reason)

    # power analysis needs a binary response
    y_binary = None
    if response_arr is not None and len(unique_resp) == 2:
        y_binary = (response_arr == unique_resp[0]).astype(int)
    elif response_arr is not None and len(unique_resp) != 2:
        log.warning("Power analysis requires a binary response "
                    "(%d unique values). Skipping.", len(unique_resp))

    power_rows = []
    power_grid_rows = []
    n_80pct = None
    n_80pct_by_delta = {}

    if y_binary is not None and not np.any(X_raw < 0):
        if args.legacy_power:
            log.info("Power analysis v1 (observed effect, several n): "
                     "%d reps x %d permutations per n", args.n_boot, args.n_perms_inner)
            n_min = max(10, n // 2)
            n_max = min(500, max(200, n * 3))
            step = max(10, (n_max - n_min) // 12)
            n_values = sorted(set(range(n_min, n_max + 1, step)) | {n})
            power_rows = dirichlet_power_analysis(
                X_raw, y_binary, n_values=n_values, n_reps=args.n_boot,
                n_perms_inner=args.n_perms_inner, pseudocount=args.pseudocount, seed=args.seed + 2,
            )
            n_80pct = _find_n_80pct(power_rows)
            if n_80pct is not None:
                log.info("  80%% power reached at n = %d", n_80pct)
            else:
                log.info("  80%% power not reached up to n = %d.", n_values[-1] if n_values else 0)
        else:
            log.info("Power analysis v2 (hypothetical chance-corrected effects)")
            delta_grid = [float(x) for x in args.delta_grid.split(",")]
            if args.n_grid:
                n_values = sorted(set(int(x) for x in args.n_grid.split(",")))
            else:
                n_values = sorted(set(max(10, int(round(n * f))) for f in [0.5, 1.0, 2.0, 4.0]))
            log.info("  n_grid=%s  delta_grid=%s  (%d reps, %d permutations per cell)",
                     n_values, delta_grid, args.power_reps, args.power_perms_inner)
            power_grid_rows = power_over_effect_grid(
                X_raw, n_values=n_values, delta_grid=delta_grid,
                n_reps=args.power_reps, n_perms_inner=args.power_perms_inner,
                n_calib_reps=args.calib_reps, calib_bisect_iter=args.calib_iters,
                seed=args.seed + 2,
            )
            grid_df = pd.DataFrame(power_grid_rows)
            for d in delta_grid:
                sub = sorted(grid_df[grid_df["target_delta_R2"] == d].to_dict("records"),
                             key=lambda r: r["n"])
                n80 = _find_n_80pct(sub)
                n_80pct_by_delta[d] = n80
                log.info("  target_delta_R2=%.4f -> n for 80%% power ~= %s", d, n80)
            # headline number: the smallest target delta that reaches 80% at
            # some n on the grid
            finite = {d: v for d, v in n_80pct_by_delta.items() if v is not None}
            if finite:
                n_80pct = finite[min(finite)]
    elif y_binary is not None:
        log.info("Power analysis skipped (negative feature values).")

    # optional check of the v1 power estimate
    sim_result = None
    if args.simulate and not args.legacy_power:
        log.warning("--simulate only checks the v1 power curve; add "
                    "--legacy-power or drop --simulate. Skipping.")
    elif args.simulate and y_binary is not None and not np.any(X_raw < 0):
        log.info("Simulation check: 50 reps at current n = %d", n)
        sim_result = simulate_at_current_n(
            X_raw, y_binary,
            n_target=n,
            n_sim_reps=50,
            n_perms_inner=99,
            pseudocount=args.pseudocount,
            seed=args.seed + 999,
        )
        log.info(
            "  Simulated power at n = %d: %.3f  95%% CI [%.3f, %.3f]  "
            "(50 reps, 99 permutations)",
            n,
            sim_result.get("power", float("nan")),
            sim_result.get("ci_lower", float("nan")),
            sim_result.get("ci_upper", float("nan")),
        )
        pw_full = next((r for r in power_rows if r["n"] == n), None)
        if pw_full:
            delta = abs(sim_result.get("power", 0) - pw_full["power"])
            log.info(
                "  Full curve at n = %d: %.3f, difference %.3f (%s)",
                n, pw_full["power"], delta,
                "within sampling noise" if delta < 0.10 else "larger than expected",
            )

    report = {
        "tool": "batch_detector.py",
        "version": __version__,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "input_file": str(args.input),
        "labels_file": str(args.labels),
        "response_col": args.response_col,
        "batch_arg": str(args.batch) if args.batch else None,
        "n_samples": n,
        "n_features": X_out.shape[1],
        "n_batch_groups": (int(len(np.unique(batch_arr)))
                           if batch_arr is not None else None),
        "clr_applied": args.clr,
        "pseudocount": args.pseudocount,
        "permanova_n_perms": args.n_perms,
        "zero_handling": args.zero_handling,
        "response_R2": res_response["R2"] if res_response else None,
        "response_F": res_response["F_stat"] if res_response else None,
        "permanova_p_response": res_response["p_value"] if res_response else None,
        "response_E0": response_E0,
        "response_delta_R2": response_delta_R2,
        "batch_R2": res_batch["R2"] if res_batch else None,
        "batch_F": res_batch["F_stat"] if res_batch else None,
        "permanova_p_batch": res_batch["p_value"] if res_batch else None,
        "batch_E0": batch_E0,
        "batch_delta_R2": batch_delta_R2,
        "batch_signal_ratio": batch_signal_ratio,
        "joint_margin_model": joint_model,
        "rho_permutation_null": rho_result,
        "recommendation_version": "legacy_v1" if args.legacy_recommendation else "v2_chance_calibrated",
        "pooling_recommendation": rec,
        "recommendation_reason": rec_reason,
        "power_analysis_version": "legacy_v1" if args.legacy_power else "v2_hypothetical_effect_grid",
        "min_n_for_80pct_power": n_80pct,
        "n_for_80pct_power_by_target_delta_R2": (n_80pct_by_delta if not args.legacy_power else None),
        "dirichlet_bootstrap_n_reps": args.n_boot if (power_rows and args.legacy_power) else None,
        "simulation_validation": sim_result if args.simulate else None,
        "permanova_response_detail": res_response,
        "permanova_batch_detail": res_batch,
    }

    json_path = out_dir / "batch_detector_report.json"
    with open(json_path, "w") as fh:
        json.dump(report, fh, indent=2)
    log.info("Report saved: %s", json_path)

    power_tsv_path = None
    if power_rows:
        power_tsv_path = out_dir / "batch_detector_power.tsv"
        pd.DataFrame(power_rows).to_csv(power_tsv_path, sep="\t", index=False)
        log.info("Power curve (v1): %s", power_tsv_path)
    if power_grid_rows:
        power_tsv_path = out_dir / "batch_detector_power_grid.tsv"
        pd.DataFrame(power_grid_rows).to_csv(power_tsv_path, sep="\t", index=False)
        log.info("Power grid (v2): %s", power_tsv_path)

    fig_path = out_dir / "batch_detector_figure.pdf"
    generate_figure(
        res_response=res_response,
        res_batch=res_batch,
        power_rows=power_rows,
        n_current=n,
        n_80pct=n_80pct,
        recommendation=rec,
        out_path=fig_path,
    )

    log.info("Summary")
    if res_response:
        log.info("  Response PERMANOVA : R2 = %.4f   %s",
                 res_response["R2"], _sig_stars(res_response["p_value"]))
    if res_batch:
        log.info("  Batch PERMANOVA    : R2 = %.4f   %s",
                 res_batch["R2"], _sig_stars(res_batch["p_value"]))
    if response_delta_R2 is not None:
        log.info("  Response delta_R2  : %+.5f  (E0=%.5f)", response_delta_R2, response_E0)
    if batch_delta_R2 is not None:
        log.info("  Batch delta_R2     : %+.5f  (E0=%.5f)", batch_delta_R2, batch_E0)
    if batch_signal_ratio is not None:
        log.info("  Batch:signal ratio : %.2fx", batch_signal_ratio)
    if rho_result is not None:
        log.info("  rho null percentile: %.1f", rho_result["percentile_of_observed"])
    log.info("  Recommendation     : %s  (%s)", rec,
             "legacy_v1" if args.legacy_recommendation else "v2_chance_calibrated")
    if n_80pct is not None:
        log.info("  Min n (80%% power)  : %d same-protocol samples", n_80pct)
    else:
        log.info("  Min n (80%% power)  : not reached on the grid")
    log.info("Outputs: %s", json_path)
    if power_tsv_path:
        log.info("         %s", power_tsv_path)
    log.info("         %s", fig_path)


if __name__ == "__main__":
    main()
