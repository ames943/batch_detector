"""
Shared statistics library for the BIBM camera-ready revision.

Implements: CLR transform + zero-handling variants, Aitchison distance,
single-factor PERMANOVA (Anderson 2001), chance-expectation E0 = (g-1)/(n-1),
a from-scratch Type-III ("by=margin") two-factor PERMANOVA with restricted
(within-cohort) permutation for the response term (Python substitute for
R vegan::adonis2(... by="margin", permutations=how(blocks=cohort))),
Hedges' g effect sizes, DerSimonian-Laird random-effects meta-analysis,
BH-FDR, and a Dirichlet parametric bootstrap (Gamma/sum-normalize identity)
used throughout for zero-distance-free resampling and power analysis.

No R/vegan is installed in this environment (checked: `Rscript` not found).
Section 2 of the revision therefore uses the Python implementation below,
as explicitly permitted by the task instructions.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy import stats

PSEUDOCOUNT_DEFAULT = 1e-6


# ─────────────────────────── CLR / zero-handling ────────────────────────────

def clr_transform(X: np.ndarray, pseudocount: float = PSEUDOCOUNT_DEFAULT) -> np.ndarray:
    lx = np.log(X + pseudocount)
    return lx - lx.mean(axis=1, keepdims=True)


def clr_pseudocount(X_raw: np.ndarray, eps: float) -> np.ndarray:
    """Additive pseudocount on raw counts, then CLR (no further pseudocount needed)."""
    Xp = X_raw + eps
    lx = np.log(Xp)
    return lx - lx.mean(axis=1, keepdims=True)


def clr_multiplicative_replacement(X_raw: np.ndarray, delta: float = 1e-5) -> np.ndarray:
    """
    Multiplicative replacement (Martin-Fernandez et al. 2003) then CLR.
    Falls back to a manual implementation if skbio is unavailable.
    """
    try:
        from skbio.stats.composition import multiplicative_replacement, clr as skbio_clr
        # skbio operates on closed compositions (rows summing to 1)
        comp = X_raw / X_raw.sum(axis=1, keepdims=True)
        comp_r = multiplicative_replacement(comp, delta=delta)
        return skbio_clr(comp_r)
    except Exception:
        comp = X_raw / X_raw.sum(axis=1, keepdims=True)
        n, p = comp.shape
        zero_mask = comp == 0
        n_zeros = zero_mask.sum(axis=1)
        out = comp.copy()
        for i in range(n):
            if n_zeros[i] == 0:
                continue
            out[i, zero_mask[i]] = delta
            nonzero = ~zero_mask[i]
            correction = 1 - n_zeros[i] * delta
            out[i, nonzero] = comp[i, nonzero] * correction
        lx = np.log(out)
        return lx - lx.mean(axis=1, keepdims=True)


def aitchison_dist(X_clr: np.ndarray) -> np.ndarray:
    return squareform(pdist(X_clr, metric="euclidean"))


# ─────────────────────────── PERMANOVA (1-factor) ───────────────────────────

def _ss_total(d2: np.ndarray) -> float:
    return float(np.sum(np.triu(d2, k=1))) / d2.shape[0]


def _ss_within(d2: np.ndarray, grp: np.ndarray) -> float:
    sw = 0.0
    for g in np.unique(grp):
        m = grp == g
        ng = int(m.sum())
        if ng < 2:
            continue
        sw += float(np.sum(np.triu(d2[np.ix_(m, m)], k=1))) / ng
    return sw


def permanova(D: np.ndarray, grouping, n_perms: int = 999, seed: int = 42) -> dict:
    """Distance-based PERMANOVA (Anderson 2001). D = (n,n) distance matrix (not squared)."""
    grp = np.asarray(grouping)
    n = D.shape[0]
    q = len(np.unique(grp))
    if q < 2:
        raise ValueError("PERMANOVA requires >=2 groups")
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
    p_val = float((perm_F >= F).sum() + 1) / (n_perms + 1)

    return dict(R2=float(R2), F_stat=float(F), p_value=float(p_val),
                n_samples=n, n_groups=q, n_perms=n_perms,
                perm_F=perm_F)


def E0(g: int, n: int) -> float:
    """Chance expectation of PERMANOVA R^2 under a random g-group label, per reviewer note."""
    return (g - 1) / (n - 1)


# ───────────────── Type-III ("margin") two-factor PERMANOVA ─────────────────

def _dummies(labels: np.ndarray) -> np.ndarray:
    """Treatment-coded dummy matrix (drop first level), for use alongside an intercept column."""
    levels = sorted(pd.unique(labels))
    Z = np.zeros((len(labels), len(levels) - 1))
    for j, lev in enumerate(levels[1:]):
        Z[:, j] = (labels == lev).astype(float)
    return Z


def _gower_center(D: np.ndarray) -> np.ndarray:
    n = D.shape[0]
    A = -0.5 * (D ** 2)
    J = np.eye(n) - np.ones((n, n)) / n
    return J @ A @ J


def _hat_trace(G: np.ndarray, X: np.ndarray) -> float:
    """trace(H G) with H = X (X'X)^+ X', computed as trace((X'X)^+ (X'GX))
    via the cyclic property of trace -- O(n k^2) instead of O(n^3), which
    matters because this is called inside permutation loops."""
    M_inv = np.linalg.pinv(X.T @ X)
    XtGX = X.T @ G @ X
    return float(np.trace(M_inv @ XtGX))


def margin_permanova_two_factor(D: np.ndarray, cohort: np.ndarray, response: np.ndarray,
                                  n_perms: int = 999, seed: int = 42) -> dict:
    """
    Type-III ("by=margin") two-factor PERMANOVA — Python substitute for
    R vegan::adonis2(dist ~ cohort + response, by="margin",
                      permutations=how(nperm=999, blocks=cohort)).

    - Marginal R^2 for each term = unique variance explained given the other term
      (SS_full - SS_reduced-without-that-term) / SS_total.
    - response p-value: permutations RESTRICTED within cohort strata (matches
      vegan's blocks=cohort — appropriate because response is nested for testing
      within the batch structure).
    - cohort p-value: permutations unrestricted across all samples (cohort is the
      blocking factor itself, so it cannot be tested under within-block restriction;
      this mirrors the plain single-factor batch PERMANOVA from Section 1).
    """
    n = D.shape[0]
    G = _gower_center(D)
    SS_total = float(np.trace(G))
    ones = np.ones((n, 1))

    Z_coh = _dummies(cohort)
    Z_resp = _dummies(response)

    X_full = np.hstack([ones, Z_coh, Z_resp])
    X_coh_only = np.hstack([ones, Z_coh])
    X_resp_only = np.hstack([ones, Z_resp])

    SS_full = _hat_trace(G, X_full)
    SS_coh_only = _hat_trace(G, X_coh_only)
    SS_resp_only = _hat_trace(G, X_resp_only)

    SS_margin_resp = SS_full - SS_coh_only
    SS_margin_coh = SS_full - SS_resp_only

    R2_resp = SS_margin_resp / SS_total
    R2_coh = SS_margin_coh / SS_total

    rng = np.random.default_rng(seed)

    # response: restricted permutation within cohort strata
    strata = {c: np.where(cohort == c)[0] for c in np.unique(cohort)}
    perm_resp = np.empty(n_perms)
    for i in range(n_perms):
        resp_p = response.copy()
        for c, idxs in strata.items():
            resp_p[idxs] = rng.permutation(response[idxs])
        Z_resp_p = _dummies(resp_p)
        X_full_p = np.hstack([ones, Z_coh, Z_resp_p])
        SS_full_p = _hat_trace(G, X_full_p)
        perm_resp[i] = (SS_full_p - SS_coh_only) / SS_total
    p_resp = float((perm_resp >= R2_resp).sum() + 1) / (n_perms + 1)

    # cohort: unrestricted permutation
    perm_coh = np.empty(n_perms)
    for i in range(n_perms):
        coh_p = rng.permutation(cohort)
        Z_coh_p = _dummies(coh_p)
        X_full_p = np.hstack([ones, Z_coh_p, Z_resp])
        SS_full_p = _hat_trace(G, X_full_p)
        perm_coh[i] = (SS_full_p - SS_resp_only) / SS_total
    p_coh = float((perm_coh >= R2_coh).sum() + 1) / (n_perms + 1)

    return dict(
        n=n, SS_total=SS_total,
        R2_response_margin=float(R2_resp), p_response_margin=p_resp,
        R2_cohort_margin=float(R2_coh), p_cohort_margin=p_coh,
        response_perm_restricted_within_cohort=True,
        cohort_perm_restricted=False,
        n_perms=n_perms,
    )


# ───────────────────────────── rho permutation null ─────────────────────────

def rho_permutation_null(D: np.ndarray, cohort: np.ndarray, response: np.ndarray,
                           n_perms: int = 999, seed: int = 42) -> dict:
    """
    Null distribution of rho = R2_cohort / R2_response under INDEPENDENT
    shuffles of cohort and response (each shuffled independently of the other
    and of the real data), 999 reps. Reports percentile of the observed rho.
    """
    d2 = D ** 2
    ST = _ss_total(d2)
    rng = np.random.default_rng(seed)

    def r2(grp):
        sw = _ss_within(d2, grp)
        return (ST - sw) / ST if ST > 0 else 0.0

    obs_r2_resp = r2(response)
    obs_r2_coh = r2(cohort)
    obs_rho = obs_r2_coh / obs_r2_resp if obs_r2_resp > 0 else np.inf

    null_rho = np.empty(n_perms)
    for i in range(n_perms):
        resp_p = rng.permutation(response)
        coh_p = rng.permutation(cohort)
        r2r = r2(resp_p)
        r2c = r2(coh_p)
        null_rho[i] = r2c / r2r if r2r > 0 else np.nan
    valid = null_rho[np.isfinite(null_rho)]
    percentile = float((valid <= obs_rho).mean() * 100) if len(valid) else float("nan")

    return dict(observed_r2_response=obs_r2_resp, observed_r2_cohort=obs_r2_coh,
                observed_rho=obs_rho, null_rho_mean=float(np.nanmean(null_rho)),
                null_rho_median=float(np.nanmedian(null_rho)),
                percentile_of_observed=percentile, n_perms=n_perms)


# ───────────────────────────── Dirichlet bootstrap ──────────────────────────

def dirichlet_draw(alpha: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    G = rng.standard_gamma(alpha)
    return G / G.sum(axis=1, keepdims=True)


# ───────────────────────────── Effect sizes / meta ──────────────────────────

def hedges_g(x1: np.ndarray, x2: np.ndarray) -> tuple:
    """Hedges' g (bias-corrected standardized mean difference) for x1 (group1) vs x2 (group2), with SE."""
    n1, n2 = len(x1), len(x2)
    if n1 < 2 or n2 < 2:
        return np.nan, np.nan
    m1, m2 = np.mean(x1), np.mean(x2)
    s1, s2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    if sp == 0:
        return 0.0, np.nan
    d = (m1 - m2) / sp
    J = 1 - 3 / (4 * (n1 + n2 - 2) - 1)
    g = J * d
    var_g = J ** 2 * ((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2 - 2)))
    return float(g), float(np.sqrt(var_g))


def dsl_meta(effects: np.ndarray, variances: np.ndarray) -> dict:
    w_fe = 1.0 / variances
    theta_fe = np.sum(w_fe * effects) / np.sum(w_fe)
    Q = float(np.sum(w_fe * (effects - theta_fe) ** 2))
    k = len(effects)
    df_Q = k - 1
    p_Q = float(1 - stats.chi2.cdf(Q, df=df_Q)) if df_Q > 0 else 1.0
    c_factor = np.sum(w_fe) - np.sum(w_fe ** 2) / np.sum(w_fe)
    tau2 = float(max(0.0, (Q - df_Q) / c_factor)) if c_factor > 0 else 0.0
    I2 = float(max(0.0, (Q - df_Q) / Q * 100)) if Q > df_Q else 0.0
    w_re = 1.0 / (variances + tau2)
    theta_re = float(np.sum(w_re * effects) / np.sum(w_re))
    se_re = float(np.sqrt(1.0 / np.sum(w_re)))
    z = theta_re / se_re if se_re > 0 else 0.0
    p_z = float(2 * (1 - stats.norm.cdf(abs(z))))
    return dict(theta_re=theta_re, se_re=se_re, ci_lo=theta_re - 1.96 * se_re,
                ci_hi=theta_re + 1.96 * se_re, tau2=tau2, I2=I2, Q=Q, Q_p=p_Q,
                z=z, p_z=p_z)


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    m = len(pvals)
    order = np.argsort(pvals)
    ranks = np.empty(m, dtype=int)
    ranks[order] = np.arange(1, m + 1)
    q = pvals * m / ranks
    q_mono = np.minimum.accumulate(q[order][::-1])[::-1]
    q_out = np.empty(m)
    q_out[order] = np.minimum(q_mono, 1.0)
    return q_out


def cohort_of_n283(sid: str) -> str:
    if sid.startswith("SRR5930"):
        return "cohort1"
    if sid.startswith("SRR11413"):
        return "cohort2"
    if sid.startswith("SRR6000"):
        return "cohort3"
    return "cohort4"
