#!/usr/bin/env python3
"""
Review2 Task 4 — ENA lookup for the six unprocessed studies + Gopalakrishnan 2018.

1. Per-run ENA metadata for every accession listed in metadata/response_sheet.xlsx
   for the six studies (+ the 22 'Gopalakrishnan_2018' runs) -> ena_per_run_response_sheet_studies.tsv
2. Whole-BioProject totals (runs, library strategy, FASTQ.gz bytes, first_public)
   for each study's BioProject -> ena_project_totals.tsv

All values come from the ENA Portal API (https://www.ebi.ac.uk/ena/portal/api/search).
fastq_bytes are compressed FASTQ.gz sizes as reported by ENA.
"""
import io
import pandas as pd
import requests

OUT = "results/revision/review2"
API = "https://www.ebi.ac.uk/ena/portal/api/search"
STUDIES = ["Gopalakrishnan_2018", "Gunjur_2024", "Heshiki_2020", "Liu_2022",
           "McCulloch_2022", "Spencer_2021", "Tsakmaklis_2023"]
PROJECTS = ["PRJEB49516", "PRJNA494824", "PRJNA866654", "PRJNA762360", "PRJNA770295",
            "PRJNA1011235", "PRJEB22893", "PRJEB22894", "PRJEB22874", "PRJEB22895"]
RUN_FIELDS = ("run_accession,study_accession,sample_accession,experiment_title,library_strategy,"
              "library_source,library_layout,instrument_model,read_count,base_count,fastq_bytes,"
              "study_title,center_name,first_public")


def gb(series):
    return series.astype(str).apply(
        lambda s: sum(float(x) for x in s.split(";") if x not in ("", "nan"))) / 1e9


def ena(query, fields):
    r = requests.post(API, data={"result": "read_run", "query": query, "fields": fields,
                                 "format": "tsv", "limit": 0}, timeout=120)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), sep="\t")


sheet = pd.read_excel("metadata/response_sheet.xlsx")
per_run = []
for d in STUDIES:
    ids = sheet.loc[sheet.dataset == d, "sampleid"].astype(str).tolist()
    for i in range(0, len(ids), 100):
        q = " OR ".join(f"run_accession={x}" for x in ids[i:i + 100])
        t = ena(q, RUN_FIELDS)
        t["dataset"] = d
        per_run.append(t)
per_run = pd.concat(per_run)
per_run["fastq_gb"] = gb(per_run.fastq_bytes)
per_run.to_csv(f"{OUT}/ena_per_run_response_sheet_studies.tsv", sep="\t", index=False)
print(per_run.groupby("dataset").agg(runs=("run_accession", "size"),
                                     projects=("study_accession", lambda x: sorted(set(x))),
                                     fastq_gb=("fastq_gb", "sum")).round(1))

rows = []
for p in PROJECTS:
    t = ena(f"study_accession={p}", "run_accession,library_strategy,library_layout,"
                                    "instrument_model,fastq_bytes,first_public")
    t["gb"] = gb(t.fastq_bytes)
    g = (t.groupby(["library_strategy", "library_layout", "instrument_model"])
          .agg(runs=("run_accession", "size"), fastq_gb=("gb", "sum")).reset_index())
    g.insert(0, "project", p)
    g["first_public"] = t.first_public.min()
    rows.append(g)
tot = pd.concat(rows)
tot["fastq_gb"] = tot.fastq_gb.round(2)
tot.to_csv(f"{OUT}/ena_project_totals.tsv", sep="\t", index=False)
print(tot.to_string(index=False))
