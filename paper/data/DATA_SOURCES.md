# Data sources and what is redistributed here

Everything in `paper/` is derived from public data. This file says where each
input came from, what the terms are, and what is (and is not) copied into this
repo. Last checked October 2026; the license wording was read from the sources
linked below. Check the terms yourself before reusing any of it.

## Included in this repo

| files | size | original source | terms | note |
|---|---|---|---|---|
| `results/ml/n283_4cohort/X_genus_{raw,clr}.tsv`, `results/ml/lee2022/X_genus_{raw,clr}.tsv`, `results/ml/n118_3cohort/X_genus_{raw,clr}.tsv`, `results/ml/phase3c/X_{phylum,species}_clr.tsv` | 0.15 to 19 MB each (45 MB total) | Genus/phylum/species abundance tables **we computed** (fastp + Kraken2) from public reads: C1 Frankel 2017 (PRJNA397906), C2 Peng 2020 (PRJNA615114), C3 Matson 2018 (PRJNA399742), C4 Lee 2022 (PRJEB43119) | Raw reads are INSDC data: "uniform policy of free and unrestricted access" ([INSDC policy](https://www.insdc.org/policy/)); NCBI "places no restrictions on the use or distribution of the data" but notes submitters may claim rights ([NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)). The Lee 2022 and Matson 2018 data statements only give the repository accession; Frankel 2017 likewise (SRP115355). | No reads, only taxon tables. **The Peng 2020 data statement could not be read** (paywalled), so for that cohort only the INSDC terms were checked. |
| `metadata/response_labels*.tsv`, `results/ml/n283_4cohort/response_labels_n283.tsv`, `metadata/lee2022_labels.tsv` | 3 to 6 KB | R / NR labels from Orletskaia & Olekhnovich, *Comput Struct Biotechnol J* 2026, doi:10.34133/csbj.0065 (preprint doi:10.1101/2025.05.07.652660); two samples (SRR6000943, ERR10290768) from the original studies' label files | Article and preprint: CC BY 4.0. The authors' label tables are also in their GitHub repo [JeniaOle13/cancer-biomarkers](https://github.com/JeniaOle13/cancer-biomarkers) (MIT) | Only run accession, response and cohort are included. **Please cite Orletskaia & Olekhnovich when reusing the labels.** |
| `metadata/lee2022_sites.tsv` | 3 KB | Recruiting site of each Lee 2022 run. Identical, for all 165 runs, to the public ENA sample attribute `Cohort` of PRJEB43119 (checked October 2026) | ENA/INSDC terms as above; Lee et al. 2022 is CC BY 4.0 | Collection dates and coordinates that ENA also exposes were not copied. |
| `metadata/download_manifest_SRR11413xxx.tsv` | 7 KB | ENA file report for PRJNA615114 | public ENA metadata | |
| `results/**` TSV, JSON, MD, PDF, PNG (small result tables, figures, notes) | < 3.5 MB each | our own analysis outputs | MIT, like the rest of the repo | `results/revision/review*/X_genus_raw_rarefied*.tsv` are rarefied versions of the abundance tables above. |
| `results/ml/tumor/` result tables (AUCs, permutation nulls, PERMANOVA summary) | < 10 KB | aggregate statistics from the tumor analysis | no patient-level data inside | the per-patient feature tables are **not** included, see below |

## The four microbiome cohorts and their original studies

| cohort | INSDC accession | original study |
|---|---|---|
| C1 | PRJNA397906 (SRA study SRP115355) | Frankel AE et al. Metagenomic shotgun sequencing and unbiased metabolomic profiling identify specific human gut microbiota and metabolites associated with immune checkpoint therapy efficacy in melanoma patients. *Neoplasia* 19(10):848-855, 2017 |
| C2 | PRJNA615114 | Peng Z et al. The gut microbiome is associated with clinical response to anti-PD-1/PD-L1 immunotherapy in gastrointestinal cancer. *Cancer Immunol Res* 8(10):1251-1261, 2020 |
| C3 | PRJNA399742 (SRA study SRP116709) | Matson V et al. The commensal microbiome is associated with anti-PD-1 efficacy in metastatic melanoma patients. *Science* 359(6371):104-108, 2018 |
| C4 | PRJEB43119 | Lee KA et al. Cross-cohort gut microbiome associations with immune checkpoint inhibitor response in advanced melanoma. *Nat Med* 28(3):535-544, 2022 |

Genus tables were derived by the author from public INSDC reads; please cite the original studies if you use them.

## Not included (and how to get them)

| files | size | original source | terms | what to do |
|---|---|---|---|---|
| Raw FASTQ files | about 1 TB | the four INSDC projects above | public | `scripts/process_batch.py` and `scripts/download_lee2022.sh` download them from ENA one sample at a time and delete them after classification |
| Kraken2 reports (`results/kraken_reports/`) and fastp JSON | 173 MB and 168 MB | our output from the reads | would be fine to share, left out for size | regenerate with the commands in `paper/README.md` |
| Tumor feature tables (`results/ml/tumor/combined_tumor_features.tsv`, `riaz2017_`, `liu2019_`, `hugo2016_features.tsv`; raw JSON in `data/cbioportal/`) | 0.1 MB (features), 95 MB (raw JSON) | cBioPortal studies `mel_iatlas_riaz_nivolumab_2017` (Riaz 2017), `mel_dfci_2019` (Liu 2019), `mel_ucla_2016` (Hugo 2016) | cBioPortal data are ODbL unless noted ([cBioPortal FAQ](https://docs.cbioportal.org/user-guide/faq/)). ODbL is share-alike, so a table built from them would have to be released under ODbL. For the Riaz study no license file was found in the cBioPortal datahub. | Left out to be safe. Rebuild with `python scripts/cbioportal_download.py` then `python scripts/process_cbioportal.py` (network needed). |
| Duvallet CRC tables (`results/batch_detector_external_validation/input_*.tsv`) | 0.55 MB | microbiomeHD, Duvallet et al. 2017, Nat Commun 8:1784; Zeller 2014, Baxter 2016, Zackular 2014, Zhao 2012 | The Zenodo record 10.5281/zenodo.840333 is **CC BY-NC 4.0**. The GitHub repo cduvallet/microbiomeHD has no license file. | Left out. `bash data/get_duvallet_crc.sh` downloads the four archives and rebuilds the tables. The rebuilt files are byte-identical to the ones used in the paper (checked with md5). |
| `metadata/response_sheet.xlsx` (624 rows, 11 studies, with a `log_ratio` column) | 28 KB | our local copy of the label curation | The authors' public repos contain `data/meta_data.tsv` (sampleid, response, dataset, cancer_type) and `article/Table_S2.xlsx`, but no file named `response_sheet.xlsx` and no `log_ratio` column was found. Provenance of our copy is not confirmed. | Left out. Only `scripts/revision/review2_task4_ena_lookup.py` reads it; to rerun it, use `data/meta_data.tsv` from the repo above (same first columns). |

## Citations for reused data

- Frankel AE et al. *Neoplasia* 19(10):848-855, 2017.
- Peng Z et al. *Cancer Immunol Res* 8(10):1251-1261, 2020.
- Matson V et al. *Science* 359(6371):104-108, 2018.
- Lee KA et al. *Nat Med* 28(3):535-544, 2022.
- Orletskaia VA, Olekhnovich EI. *Comput Struct Biotechnol J*, 2026, doi:10.34133/csbj.0065.
- Duvallet C et al. *Nat Commun* 8:1784, 2017.
- Riaz N et al. *Cell* 171(4):934-949, 2017; Liu D et al. *Nat Med* 25:1916-1927, 2019; Hugo W et al. *Cell* 165(1):35-44, 2016; Gao J et al. *Sci Signal* 6(269):pl1, 2013 (cBioPortal).
