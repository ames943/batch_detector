#!/usr/bin/env python3
"""
Prepare batch_detector.py inputs from microbiomeHD CRC cohorts.

Downloads 4 colorectal cancer cohorts from Duvallet et al. 2017 (microbiomeHD):
  - Zeller 2014  (n=156, CRC=53, H=88)
  - Baxter 2016  (n=490, CRC=120, H=172)
  - Zackular 2014 (n=90, CRC=30, H=30)
  - Zhao 2012    (n=102, CRC=46, H=56)

Produces:
  input_feature_matrix.tsv — samples × genus counts (raw, for CLR in batch_detector)
  input_labels.tsv         — sample_id, response (1=CRC, 0=H), study (batch)
"""

import re
import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path("microbiomeHD_crc")
OUT = Path("results/batch_detector_external_validation")
OUT.mkdir(parents=True, exist_ok=True)

# ── Dataset configuration ─────────────────────────────────────────────────────
DATASETS = {
    "zeller": {
        "dir": "crc_zeller_results",
        "prefix": "crc_zeller",
        "meta_id_col": "Subject ID",
        "rdp_file": "RDP/RDP_classifications.denovo.txt",
    },
    "baxter": {
        "dir": "crc_baxter_results",
        "prefix": "crc_baxter",
        "meta_id_col": "Sample_Name_s",
        "rdp_file": "RDP/RDP_classifications.denovo.txt",
    },
    "zackular": {
        "dir": "crc_zackular_results",
        "prefix": "crc_zackular",
        "meta_id_col": "sample_id",
        "rdp_file": "RDP/RDP_classifications.denovo.txt",
    },
    "zhao": {
        "dir": "crc_zhao_results",
        "prefix": "crc_zhao",
        "meta_id_col": "#SampleID",
        "rdp_file": "RDP/RDP_classifications.denovo.txt",
    },
}

RDP_CONF_THRESHOLD = 0.8  # minimum confidence to accept genus assignment


def parse_rdp(rdp_path: Path) -> dict:
    """
    Parse RDP fixed-column classifier output → OTU_ID: genus_name.

    RDP output format (tab-separated):
      OTU_ID  <blank>  Root  rootrank  conf  Bacteria  domain  conf  ... genus  genus  conf

    We look for the column where rank == 'genus' and grab the preceding name.
    Returns {} entry if confidence < threshold or genus not found.
    """
    otu2genus = {}
    with open(rdp_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            otu_id = parts[0]
            genus_name = None
            conf = 0.0
            # Walk through triplets after the first two columns (otu_id, blank)
            i = 2
            while i + 2 < len(parts):
                name = parts[i]
                rank = parts[i + 1]
                try:
                    c = float(parts[i + 2])
                except (ValueError, IndexError):
                    c = 0.0
                if rank == "genus":
                    genus_name = re.sub(r'"', "", name).strip()
                    conf = c
                    break
                i += 3
            if genus_name and conf >= RDP_CONF_THRESHOLD:
                otu2genus[otu_id] = genus_name
    print(f"  RDP: {len(otu2genus)} OTUs assigned to genus at conf ≥ {RDP_CONF_THRESHOLD}")
    return otu2genus


def load_otu_genus(ds_dir: Path, prefix: str, otu2genus: dict) -> pd.DataFrame:
    """
    Load OTU table, collapse to genus level.
    OTU table is OTU_ID × samples; we transpose and sum by genus.
    Returns DataFrame: samples × genus (raw counts, float).
    """
    otu_path = ds_dir / f"{prefix}.otu_table.100.denovo"
    otu = pd.read_csv(otu_path, sep="\t", index_col=0)
    otu = otu.apply(pd.to_numeric, errors="coerce").fillna(0)
    # otu is OTU_ID × samples — transpose → samples × OTU_ID
    otu = otu.T
    otu.index = otu.index.astype(str)

    # Map OTU IDs to genera (drop OTUs without a genus assignment)
    assigned = {otu_id: g for otu_id, g in otu2genus.items() if otu_id in otu.columns}
    otu_sub = otu[[c for c in otu.columns if c in assigned]]
    otu_sub = otu_sub.rename(columns=assigned)

    # Sum OTU counts sharing the same genus
    genus_df = otu_sub.groupby(level=0, axis=1).sum()
    print(f"  OTU table: {len(otu)} samples, {genus_df.shape[1]} genera after collapse")
    return genus_df


def load_meta(ds_dir: Path, prefix: str, meta_id_col: str) -> pd.DataFrame:
    """
    Load metadata, return df with index = sample_id (str) and 'response' column.
    response: 1 = CRC, 0 = H, NaN = excluded (nonCRC / other).
    """
    meta_path = ds_dir / f"{prefix}.metadata.txt"
    meta = pd.read_csv(meta_path, sep="\t", encoding="latin-1")
    meta[meta_id_col] = meta[meta_id_col].astype(str)
    meta = meta.set_index(meta_id_col)

    # Map DiseaseState → binary response
    def _map(state):
        if state == "CRC":
            return 1
        if state == "H":
            return 0
        return np.nan

    meta["response"] = meta["DiseaseState"].map(_map)
    return meta[["response"]]


# ── Main processing loop ──────────────────────────────────────────────────────

all_genus_dfs = []
all_label_dfs = []

for study, cfg in DATASETS.items():
    print(f"\n── {study} ────────────────────────────────")
    ds_dir = BASE / cfg["dir"]

    rdp_path = ds_dir / cfg["rdp_file"]
    otu2genus = parse_rdp(rdp_path)

    genus_df = load_otu_genus(ds_dir, cfg["prefix"], otu2genus)
    meta_df = load_meta(ds_dir, cfg["prefix"], cfg["meta_id_col"])

    # Align samples (intersection of OTU table and metadata)
    common = genus_df.index.intersection(meta_df.index)
    genus_df = genus_df.loc[common]
    meta_df = meta_df.loc[common]

    # Drop samples with no response label (nonCRC)
    mask = meta_df["response"].notna()
    genus_df = genus_df.loc[mask]
    meta_df = meta_df.loc[mask]

    # Add study label for batch
    meta_df = meta_df.copy()
    meta_df["study"] = study

    # Prefix sample IDs with study name to avoid collisions
    genus_df.index = [f"{study}__{sid}" for sid in genus_df.index]
    meta_df.index = [f"{study}__{sid}" for sid in meta_df.index]

    n_crc = int(meta_df["response"].sum())
    n_h = int((meta_df["response"] == 0).sum())
    print(f"  After filtering: {len(genus_df)} samples  (CRC={n_crc}, H={n_h})")
    print(f"  Genera: {genus_df.shape[1]}")

    all_genus_dfs.append(genus_df)
    all_label_dfs.append(meta_df)

# ── Combine and align genera ───────────────────────────────────────────────────

print("\n── Combining datasets ─────────────────────────────────")
feature_matrix = pd.concat(all_genus_dfs, axis=0, join="outer").fillna(0)
labels = pd.concat(all_label_dfs, axis=0)

print(f"Combined: {feature_matrix.shape[0]} samples × {feature_matrix.shape[1]} genera")
print(f"Labels: CRC={int(labels['response'].sum())}, H={int((labels['response']==0).sum())}")
print(f"Study distribution:\n{labels['study'].value_counts().to_string()}")

# Drop genera with zero variance across all samples
feat_var = feature_matrix.var(axis=0)
feature_matrix = feature_matrix.loc[:, feat_var > 0]
print(f"After removing zero-variance genera: {feature_matrix.shape[1]} genera")

# ── Save ───────────────────────────────────────────────────────────────────────

fm_path = OUT / "input_feature_matrix.tsv"
lb_path = OUT / "input_labels.tsv"

feature_matrix.index.name = "sample_id"
labels.index.name = "sample_id"

feature_matrix.to_csv(fm_path, sep="\t")
labels.to_csv(lb_path, sep="\t")

print(f"\nSaved:\n  {fm_path}\n  {lb_path}")
