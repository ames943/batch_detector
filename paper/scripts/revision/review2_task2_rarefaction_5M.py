#!/usr/bin/env python3
"""
Review2 Task 2 — milder rarefaction to 5,000,000 reads (vs. the first
rarefaction's 2,098,257, which discarded ~90% of most samples' reads).
Re-parses raw Kraken2 reports for integer genus counts (same logic as
review_task4b_rarefaction.py -- X_genus_raw.tsv stores percentages, not
counts, so cannot be rarefied directly). Samples below 5,000,000 total
genus-assigned reads are dropped (rarefaction cannot upsample).
"""
import glob, os, sys, time
import numpy as np
import pandas as pd

sys.path.insert(0, "scripts/revision")
from lib import clr_transform, aitchison_dist, permanova, E0

OUT = "results/revision/review2"
SEED = 42
N_PERMS = 999
DEPTH = 5_000_000
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
totals = pd.Series({s: sum(d.values()) for s, d in sample_counts.items()})
totals_named = totals.to_frame("total_genus_reads")
totals_named["cohort"] = [cohort_of(s) for s in totals_named.index]

below = totals_named[totals_named["total_genus_reads"] < DEPTH].sort_values("total_genus_reads")
below.to_csv(f"{OUT}/dropped_below_5M.tsv", sep="\t")
print(f"\nSamples below {DEPTH:,} reads (to be dropped): {len(below)}")
print("By cohort:")
print(below["cohort"].value_counts())

rng = np.random.default_rng(SEED)
rarefied_rows = {}
for sid, counts in sample_counts.items():
    total = sum(counts.values())
    if total < DEPTH:
        continue
    genera = list(counts.keys())
    colors = np.array([counts[g] for g in genera], dtype=np.int64)
    sub_counts = rng.multivariate_hypergeometric(colors, DEPTH)
    rarefied_rows[sid] = dict(zip(genera, sub_counts))

print(f"\nRarefied {len(rarefied_rows)} samples to depth={DEPTH:,}; dropped {len(below)}.")

raw_rare = pd.DataFrame.from_dict(rarefied_rows, orient="index", columns=all_genera).fillna(0.0)
raw_rare.index.name = "run_accession"
raw_rare.to_csv(f"{OUT}/X_genus_raw_rarefied_5M.tsv", sep="\t")

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
    dict(dataset=f"genus_rarefied5M_n{n}", factor="response", groups=2,
         R2=round(res_r["R2"], 6), E0=round(e0_r, 6), delta_R2=round(res_r["R2"] - e0_r, 6),
         p_value=res_r["p_value"]),
    dict(dataset=f"genus_rarefied5M_n{n}", factor="cohort", groups=coh.nunique(),
         R2=round(res_c["R2"], 6), E0=round(e0_c, 6), delta_R2=round(res_c["R2"] - e0_c, 6),
         p_value=res_c["p_value"]),
]
out_df = pd.DataFrame(rows)
out_df.to_csv(f"{OUT}/rarefied_5M_permanova.tsv", sep="\t", index=False)
print("\n" + out_df.to_string(index=False))
print(f"\nFinal n={n} (dropped {len(below)} of {len(sample_counts)} parsed samples)")
print(f"[{time.strftime('%H:%M:%S')}] Done.")
