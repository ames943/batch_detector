#!/usr/bin/env python3
"""
Review Task 4a — per-cohort median depth (post-fastp) and median Kraken2
classification rate (excluding human), for all 283 samples.
Source: results/fastp/**/*_fastp.json, results/kraken_reports/**/*_report.txt.
"""
import json, re
from pathlib import Path
import pandas as pd
import numpy as np

OUT = "results/revision/review"

labels283 = pd.read_csv("results/ml/n283_4cohort/response_labels_n283.tsv", sep="\t")


def cohort_of(sid):
    if sid.startswith("SRR5930"): return "C1"
    if sid.startswith("SRR11413"): return "C2"
    if sid.startswith("SRR6000"): return "C3"
    return "C4"


labels283["cohort"] = labels283["run_accession"].apply(cohort_of)

FASTP_DIRS = ["results/fastp", "results/fastp/lee2022"]
KRAKEN_DIRS = ["results/kraken_reports", "results/kraken_reports/lee2022"]


def find_file(sid, dirs, suffix):
    for d in dirs:
        p = Path(d) / f"{sid}{suffix}"
        if p.exists():
            return p
    return None


def parse_kraken(path):
    unclassified_pct = None
    homo_pct = 0.0
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            pct, n_reads_clade, n_reads_direct, rank, taxid, name = parts[:6]
            if rank == "U":
                unclassified_pct = float(pct)
            if taxid.strip() == "9605":  # Homo (genus)
                homo_pct = float(pct)
    return unclassified_pct, homo_pct


rows = []
missing = []
for _, row in labels283.iterrows():
    sid = row["run_accession"]
    fastp_p = find_file(sid, FASTP_DIRS, "_fastp.json")
    kraken_p = find_file(sid, KRAKEN_DIRS, "_report.txt")
    if fastp_p is None or kraken_p is None:
        missing.append(sid)
        continue
    fj = json.load(open(fastp_p))
    total_reads = fj["summary"]["after_filtering"]["total_reads"]
    unclassified_pct, homo_pct = parse_kraken(kraken_p)
    classified_pct_excl_human = 100.0 - unclassified_pct - homo_pct
    rows.append(dict(run_accession=sid, cohort=row["cohort"],
                      total_reads_after_fastp=total_reads,
                      read_pairs_after_fastp=total_reads / 2,
                      pct_unclassified=unclassified_pct,
                      pct_homo=homo_pct,
                      pct_classified_excl_human=classified_pct_excl_human))

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/depth_classification_per_sample.tsv", sep="\t", index=False)

summary = df.groupby("cohort").agg(
    n=("run_accession", "size"),
    median_read_pairs=("read_pairs_after_fastp", "median"),
    median_total_reads=("total_reads_after_fastp", "median"),
    median_pct_classified_excl_human=("pct_classified_excl_human", "median"),
).reset_index()
summary.to_csv(f"{OUT}/depth_classification_summary.tsv", sep="\t", index=False)

print(summary.to_string(index=False))
print(f"\nTotal samples processed: {len(df)} / {len(labels283)}")
if missing:
    print(f"Missing fastp/kraken files for {len(missing)} samples: {missing[:20]}")
else:
    print("No missing files.")
