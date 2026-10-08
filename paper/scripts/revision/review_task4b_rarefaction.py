#!/usr/bin/env python3
"""
Review Task 4b — rarefy all available samples to a common genus-level read
depth, rerun pooled genus PERMANOVA (response + cohort).

X_genus_raw.tsv (used everywhere else in this project) stores Kraken2's
PERCENTAGE column, not integer read counts -- true rarefaction (random
subsampling without replacement to a common depth) requires integer counts,
so this script re-parses the raw Kraken2 reports using column 2
(n_reads_clade, the standard integer count used for genus-level abundance in
this pipeline) with the SAME genus name-cleaning and Homo-exclusion logic as
scripts/build_matrix_3cohort.py, for direct comparability.

Source: results/kraken_reports/**/*_report.txt (278/283 available -- see
review_task4_depth_classification.py for the 5 missing).
"""
import glob, os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, "scripts/revision")
from lib import clr_transform, aitchison_dist, permanova, E0

OUT = "results/revision/review"
SEED = 42
N_PERMS = 999
EXCLUDED = {"Homo"}

labels283 = pd.read_csv("results/ml/n283_4cohort/response_labels_n283.tsv", sep="\t").set_index("run_accession")


def cohort_of(sid):
    if sid.startswith("SRR5930"): return "cohort1"
    if sid.startswith("SRR11413"): return "cohort2"
    if sid.startswith("SRR6000"): return "cohort3"
    return "cohort4"


def clean_name(stripped):
    if stripped.startswith("Candidatus "):
        return "Candidatus_" + stripped.split()[1]
    elif " (" in stripped:
        return stripped.split(" (")[0]
    else:
        return stripped.split()[0]


def parse_counts(fp):
    genus_counts = {}
    with open(fp) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6 or parts[3] != "G":
                continue
            n_reads_clade = int(parts[1].strip())
            name = clean_name(parts[5].strip())
            if not name or name in EXCLUDED:
                continue
            genus_counts[name] = genus_counts.get(name, 0) + n_reads_clade
    return genus_counts


report_files = sorted(glob.glob("results/kraken_reports/*_report.txt") +
                       glob.glob("results/kraken_reports/lee2022/*_report.txt"))
sample_counts = {}
for fp in report_files:
    sid = os.path.basename(fp).replace("_report.txt", "")
    if sid not in labels283.index:
        continue
    sample_counts[sid] = parse_counts(fp)

print(f"[{time.strftime('%H:%M:%S')}] Parsed integer genus counts for {len(sample_counts)}/{len(labels283)} samples", flush=True)

all_genera = sorted(set(g for d in sample_counts.values() for g in d))
totals = pd.Series({s: sum(d.values()) for s, d in sample_counts.items()}).sort_values()
totals.to_csv(f"{OUT}/rarefaction_per_sample_totals.tsv", sep="\t", header=["total_genus_reads"])
print(totals.describe())
print("\nLowest 10 totals:")
print(totals.head(10))

# ── Choose a sensible rarefaction depth ─────────────────────────────────────
# The true minimum is used unless it is a clear low-depth outlier (>3x below
# the next lowest value), in which case that sample is reported as dropped
# and the next value is used instead. Stated explicitly either way.
sorted_vals = totals.values
if len(sorted_vals) > 1 and sorted_vals[0] * 3 < sorted_vals[1]:
    DEPTH = int(sorted_vals[1])
    outlier = totals.index[0]
    print(f"\nLowest total ({sorted_vals[0]}, sample {outlier}) is a >3x outlier below "
          f"the next value ({sorted_vals[1]}) -- using the SECOND-lowest total as the "
          f"rarefaction depth and dropping the outlier sample.")
else:
    DEPTH = int(sorted_vals[0])
    outlier = None
    print(f"\nUsing the true minimum total genus-read count as the rarefaction depth: {DEPTH}")

rng = np.random.default_rng(SEED)
rarefied_rows = {}
dropped = []
for sid, counts in sample_counts.items():
    total = sum(counts.values())
    if total < DEPTH:
        dropped.append(sid)
        continue
    genera = list(counts.keys())
    colors = np.array([counts[g] for g in genera], dtype=np.int64)
    # Efficient multivariate-hypergeometric subsampling without replacement
    # to exactly DEPTH reads (does not materialize a per-read array).
    sub_counts = rng.multivariate_hypergeometric(colors, DEPTH)
    rarefied_rows[sid] = dict(zip(genera, sub_counts))

print(f"\nRarefied {len(rarefied_rows)} samples to depth={DEPTH}; dropped {len(dropped)} "
      f"below-depth samples: {dropped}")

raw_rare = pd.DataFrame.from_dict(rarefied_rows, orient="index", columns=all_genera).fillna(0.0)
raw_rare.index.name = "run_accession"
raw_rare.to_csv(f"{OUT}/X_genus_raw_rarefied.tsv", sep="\t")

resp = labels283.reindex(raw_rare.index)["response"]
coh = pd.Series([cohort_of(s) for s in raw_rare.index], index=raw_rare.index)

clr_rare = clr_transform(raw_rare.values, pseudocount=1e-6)
D = aitchison_dist(clr_rare)

res_r = permanova(D, resp.values, n_perms=N_PERMS, seed=SEED)
res_c = permanova(D, coh.values, n_perms=N_PERMS, seed=SEED + 1)
n = len(raw_rare)
e0_r = E0(2, n)
e0_c = E0(coh.nunique(), n)

rows = [
    dict(dataset=f"genus_rarefied_n{n}", factor="response", groups=2,
         R2=round(res_r["R2"], 6), E0=round(e0_r, 6), delta_R2=round(res_r["R2"] - e0_r, 6),
         p_value=res_r["p_value"]),
    dict(dataset=f"genus_rarefied_n{n}", factor="cohort", groups=coh.nunique(),
         R2=round(res_c["R2"], 6), E0=round(e0_c, 6), delta_R2=round(res_c["R2"] - e0_c, 6),
         p_value=res_c["p_value"]),
]
out_df = pd.DataFrame(rows)
out_df.to_csv(f"{OUT}/rarefied_permanova.tsv", sep="\t", index=False)
print("\n" + out_df.to_string(index=False))
print(f"\nn_rarefied={n}, n_dropped={len(dropped)} (of {len(sample_counts)} with reports; "
      f"{283-len(sample_counts)} samples had no report at all -- see task4a)")
print(f"[{time.strftime('%H:%M:%S')}] Done.")
