"""Make the small synthetic dataset used in the examples and the tests.

120 samples from 3 cohorts, 50 features. Every cohort has its own mean
composition, so there is a strong cohort effect. The response label is a coin
flip that has nothing to do with the features or the cohort.

    python make_example_data.py            # writes into this folder
    python make_example_data.py some/dir   # writes somewhere else
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 2026
N_PER_COHORT = 40
N_FEATURES = 50
COHORTS = ["cohortA", "cohortB", "cohortC"]


def make_data(seed=SEED):
    rng = np.random.default_rng(seed)
    features = [f"taxon_{i:02d}" for i in range(1, N_FEATURES + 1)]

    # shared baseline composition, then a different shift of the log
    # abundances for each cohort
    baseline = rng.normal(0, 1.5, N_FEATURES)
    cohort_shift = {c: rng.normal(0, 0.8, N_FEATURES) for c in COHORTS}

    rows, ids, cohort = [], [], []
    for c in COHORTS:
        for i in range(N_PER_COHORT):
            # sample-to-sample variation inside a cohort
            logit = baseline + cohort_shift[c] + rng.normal(0, 0.5, N_FEATURES)
            p = np.exp(logit)
            p /= p.sum()
            depth = int(rng.integers(15000, 30000))
            rows.append(rng.multinomial(depth, p))
            ids.append(f"{c}_{i + 1:02d}")
            cohort.append(c)

    X = pd.DataFrame(rows, index=pd.Index(ids, name="sample_id"), columns=features)
    response = rng.choice(["R", "NR"], size=len(ids), p=[0.45, 0.55])
    meta = pd.DataFrame({"cohort": cohort, "response": response},
                        index=pd.Index(ids, name="sample_id"))
    return X, meta


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    out.mkdir(parents=True, exist_ok=True)
    X, meta = make_data()
    X.to_csv(out / "example_features.tsv", sep="\t")
    meta.to_csv(out / "example_metadata.tsv", sep="\t")
    print(f"wrote {X.shape[0]} samples x {X.shape[1]} features to {out}")
