# Review-2 tasks: three small reruns + data-availability lookup

Rules followed: new outputs in `results/revision/review2/`; every number
appended to `SUMMARY.md` Section 15 with its source file.

---

## 1. PERMDISP for recruiting site in C4

Script: `scripts/revision/review2_task1_permdisp_c4_site.py`. Same settings
as Task 4 (genus level, Aitchison distance, CLR pseudocount=1e-6, 999
permutations). PERMDISP implementation copied verbatim from
`scripts/steps_1_4_c4_integration.py` (Anderson 2006: PCoA + distance-to-
group-centroid + permutation) for methodological consistency with the rest
of the project. Output: `results/revision/review2/permdisp_c4_site.tsv`.

| n | groups | F | p | homogeneous |
|---|---|---|---|---|
| 165 | 5 | 2.2805 | 0.120 | Yes (p>0.05) |

Per-site mean dispersion (distance to centroid): Barcelona=41.99, Leeds=43.03,
Manchester=46.34, PRIMM-NL=47.00, PRIMM-UK=46.46 — Barcelona and Leeds sit
somewhat lower than the other three, but not significantly so (p=0.120).
Unlike the batch/cohort PERMDISP results elsewhere in this project (which are
consistently heterogeneous, p<0.05), site dispersion within C4 is homogeneous
— the site PERMANOVA R² (Task 4: 9.40%, p=0.001) reflects genuine
centroid separation between sites, not a dispersion artifact.

---

## 2. Second, milder rarefaction (5,000,000 reads)

Script: `scripts/revision/review2_task2_rarefaction_5M.py`. Same re-parsing
of raw Kraken2 integer counts as the first rarefaction (`X_genus_raw.tsv`
stores percentages, not counts — see Section 13/Task 4b). This time rarefied
to **5,000,000 reads** instead of 2,098,257.

**29 samples dropped** for falling below 5,000,000 total genus-assigned
reads, **all 29 from C4** (cohort4/Lee 2022) — matches the expectation of
"about 30, mostly C4" almost exactly. List: `results/revision/review2/dropped_below_5M.tsv`.
Rarefied matrix: `results/revision/review2/X_genus_raw_rarefied_5M.tsv`.

Pooled PERMANOVA on the remaining 254 samples
(`results/revision/review2/rarefied_5M_permanova.tsv`):

| dataset | factor | groups | R2 | E0 | delta_R2 | p |
|---|---|---|---|---|---|---|
| genus_rarefied5M_n254 | response | 2 | 0.004382 | 0.003953 | +0.000429 | 0.186 |
| genus_rarefied5M_n254 | cohort | 4 | 0.095608 | 0.011858 | **+0.083751** | 0.001 |

Final n=254 (283 parsed, 29 dropped). Same conclusion as the deeper (2.1M)
rarefaction: response ΔR² stays near zero and non-significant (p=0.186),
cohort ΔR² stays large and significant (p=0.001) — robust across two very
different rarefaction depths (2.1M discarding ~90% of reads vs. 5M discarding
~0-70% depending on the sample), not an artifact of how aggressively the data
were rarefied.

---

## 3. The one-sample discrepancies

Compared this project's C3 (n=39) and C4 (n=165) run-accession lists directly
against `metadata/response_sheet.xlsx`'s `Matson_2019` (38 rows) and
`Lee_2022` (164 rows) subsets (same sample-ID scheme, run accessions, in both
— confirmed before comparing).

**C3**: `SRR6000943` is in this project's C3 but not in the response_sheet's
38-sample Matson_2019 list. Its label (R) came from
`metadata/response_labels_PRJNA399742_extended.tsv` — this project's own
Matson label file, sourced independently of the Orletskaia compilation
(patient ID "P13", "Responder"). Also present, identically, in
`response_labels_PRJNA399742.tsv` and `response_labels_3cohort.tsv`. This
sample was evidently excluded from Orletskaia & Olekhnovich's later 38-sample
curation of Matson 2018 for a reason not stated in any file available to this
project; it was not excluded by this project's own original labeling.

**C4**: `ERR10290768` is in this project's C4 but not in the response_sheet's
164-sample Lee_2022 list. Its label (NR) came from
`metadata/lee2022_labels.tsv` — this project's own Lee 2022 label file, and
its site (PRIMM-UK) is present in `metadata/lee2022_sites.tsv`. Same
situation as the C3 case: this project's own, independently-sourced label
file includes it; the Orletskaia compilation's 164-sample Lee_2022 subset
does not, for a reason not documented anywhere available to this project.

Both extra samples' labels trace to **this project's own per-cohort label
files, not to the Orletskaia compilation** — i.e., they were correctly
labeled by an independent source before the Orletskaia harmonization was
ever consulted, and nothing suggests they are mislabeled; they simply fall
outside the 624-sample curation's own (undocumented, from this project's
vantage point) inclusion criteria.

**Rerun with both dropped (n=283 -> n=281)**:
`results/revision/review2/permanova_n281_dropped2.tsv`.

| dataset | factor | groups | R2 | E0 | delta_R2 | p |
|---|---|---|---|---|---|---|
| genus_n281_2sample_dropped | response | 2 | 0.003868 | 0.003571 | +0.000296 | 0.305 |
| genus_n281_2sample_dropped | cohort | 4 | 0.105186 | 0.010714 | **+0.094471** | 0.001 |

Essentially unchanged from the full n=283 result (response ΔR²=+0.00033,
p=0.284; cohort ΔR²=+0.09499, p=0.001 — Section 1 of `SUMMARY.md`). The
one-sample discrepancies are not a confound.

---

## 4. Data availability for the six unprocessed studies (+ Gopalakrishnan 2018)

Script: `scripts/revision/review2_task4_ena_lookup.py`. Outputs:
`results/revision/review2/ena_per_run_response_sheet_studies.tsv` (per-run ENA
metadata for every accession `metadata/response_sheet.xlsx` lists for these
studies) and `results/revision/review2/ena_project_totals.tsv` (whole-BioProject
totals). Sizes are ENA-reported compressed FASTQ.gz bytes. Every claim below
carries the URL it came from. Papers' full text was read through Europe PMC /
NCBI E-utilities (the PMC web pages sit behind a reCAPTCHA); the URLs given
are the human-readable article links.

**Bottom line: all six studies, and Gopalakrishnan's shotgun data, have public
raw reads.** None is "controlled access" or "not deposited". Gopalakrishnan's
EGA entry is exome data only, not microbiome.

| Study | Raw-read accession | Access | Seq. type | FMT / combination | Whole-BioProject shotgun size | Subset in `response_sheet.xlsx` |
|---|---|---|---|---|---|---|
| Gunjur 2024 | PRJEB49516 (= ERP134027) | Public (ENA) | Shotgun, NovaSeq 6000 2x150 | **Combination** (nivolumab + ipilimumab, phase 2 CA209-538, NCT02923934); **no FMT** | 107 runs, 348.8 GB | 106 runs, 348.8 GB |
| Heshiki 2020 | PRJNA494824 | Public (ENA) | Shotgun, HiSeq 1500 PE100 | Chemo (n=15) or **chemo + immunotherapy** (n=11); **no FMT**; 8 cancer types, not an anti-PD-1 trial | 71 runs, 304.0 GB | 11 runs, 42.4 GB |
| Liu 2022 | PRJNA866654 | Public (ENA) | Shotgun, HiSeq (ENA: HiSeq 4000) 2x150 | NSCLC "immunotherapy" (agents not stated in the text read); **no FMT**; not a trial | 14 runs, 69.8 GB | 14 runs, 69.8 GB |
| McCulloch 2022 | PRJNA762360 | Public (SRA/ENA) | Shotgun (94 runs, NovaSeq 6000) **plus** 16S (80 runs) | Anti-PD-1 alone, or **pembrolizumab + peg-IFN** (14 pts, separate trial HCC 13-105); **no FMT** in this study | 94 shotgun runs, 397.5 GB (+ 80 16S runs, 1.8 GB) | 27 runs, 123.7 GB |
| Spencer 2021 | PRJNA770295 | Public (SRA/ENA) | Shotgun (309 runs) **plus** 16S (498 runs) | Observational diet/probiotic survey in ICB patients; **no human FMT** (FMT only in mouse experiments) | 309 shotgun runs, 557.4 GB (+ 498 16S runs, 7.5 GB) | 134 runs, 271.2 GB |
| Tsakmaklis 2023 | PRJNA1011235 ("ANIMA") | Public (SRA/ENA) | Shotgun, NovaSeq 6000 | Prospective, non-interventional; anti-PD-1 (15), **anti-PD-1 + anti-CTLA-4 (13)**, anti-CTLA-4 (1); **no FMT** | 29 runs, 92.8 GB | 29 runs, 92.8 GB |

### Per-study sources

**Gunjur 2024** — Nat Med, "A gut microbial signature for combination immune
checkpoint blockade across cancer types", doi:10.1038/s41591-024-02823-z.
- Data availability (verbatim): "All CA209-538 fecal shotgun metagenomic
  sequencing data (after first-pass human decontamination) have been deposited
  to the European Nucleotide Archive (study accession no. ERP134027)."
  https://europepmc.org/article/PMC/PMC10957475
- ERP134027 = PRJEB49516, 107 runs, first public 2024-01-03:
  https://www.ebi.ac.uk/ena/browser/view/PRJEB49516
- Combination ICB (nivolumab 3 mg/kg + ipilimumab 1 mg/kg), trial CA209-538,
  NCT02923934; 106 discovery-cohort baseline stool samples; NovaSeq 6000
  2x150: https://europepmc.org/article/PMC/PMC10957475
- FMT: zero mentions of FMT/microbiota transplant in the full text (keyword
  count on the Europe PMC XML).
- Note: the same statement lists two *other* cohorts it reanalysed under EGA
  accessions (EGAS00001006982 Simpson 2022; EGAD00001006734 Andrews 2021) and
  thanks the data owners "for permission to access" them. Those are outside
  the six studies asked about, and EGA itself was not checked.

**Heshiki 2020** — *Microbiome* 8:28 (not *ISME J*), "Predictable modulation of
cancer treatment outcomes by the gut microbiota", doi:10.1186/s40168-020-00811-2.
- Data availability (verbatim): "The shotgun metagenomic sequences have been
  deposited in the European Nucleotide Archive under accession number
  PRJNA494824." https://europepmc.org/article/PMC/PMC7059390
- 71 fecal samples from 26 patients (31 baseline + 40 on-treatment), HiSeq
  1500 PE100; 26 patients: chemotherapy only (n=15) or chemo + immunotherapy
  (n=11); eight cancer types:
  https://europepmc.org/article/PMC/PMC7059390
- ENA: PRJNA494824 has 71 runs, 304.0 GB:
  https://www.ebi.ac.uk/ena/browser/view/PRJNA494824
- FMT: zero mentions in the full text.
- Only 11 of the 71 runs are in the compilation; which 11 and why is not
  stated in anything accessible.

**Liu 2022** — *Cancers* 14(21):5401, "Exploring Gut Microbiome in Predicting
the Efficacy of Immunotherapy in Non-Small Cell Lung Cancer",
doi:10.3390/cancers14215401.
- Data availability (verbatim): "Raw data used for analysis has been uploaded
  to the (NCBI) BioSample database under BioProject ID PRJNA866654, currently
  pending review and publication." The reads are now public:
  ENA lists PRJNA866654 with first_public 2022-10-01, 14 paired WGS runs
  (69.8 GB) plus one 0-GB single-end "OTHER" run:
  https://www.ebi.ac.uk/ena/browser/view/PRJNA866654
  Paper: https://europepmc.org/article/PMC/PMC9656313
- The 14 patients ("DS1") are NSCLC patients "who received immunotherapy";
  the specific agents and whether chemotherapy was combined are **not stated
  in the extracted text, so neither can be confirmed.** The paper's second
  dataset (DS2, 65 samples) comes from other public SRA projects that were
  not identified; PRJNA866654 holds DS1 only.
- FMT: zero mentions. Not a trial: samples selected from "our previous study".
- Note: the compilation tags this study cancer_type "Other"; it is NSCLC.

**McCulloch 2022** — Nat Med 28:545-556, "Intestinal microbiota signatures of
clinical response and immune-related adverse events in melanoma patients
treated with anti-PD-1", doi:10.1038/s41591-022-01698-2.
- Data availability (verbatim): "All sequencing (human and microbiome) data
  and de-identified metadata that support the findings have been deposited in
  NCBI databases and are all accessible via BioProject accession no.
  PRJNA762360." https://europepmc.org/article/PMC/PMC10246505
- ENA: 94 WGS runs (397.5 GB, NovaSeq 6000) + 80 amplicon runs (1.8 GB):
  https://www.ebi.ac.uk/ena/browser/view/PRJNA762360
- Treatment: "single-agent anti-PD-1 immunotherapy (nivolumab, pembrolizumab
  or investigational anti-PD-1) or pembrolizumab in combination with peg-IFN
  in the context of a separate clinical trial" (HCC 13-105); results text
  gives 49 anti-PD-1-alone vs 14 + peg-IFN patients (that analysis's cohort).
  Same article link.
- FMT: the 8 keyword hits are background citations of two prior FMT trials,
  not an intervention in this study.
- Also: the statement lists the Houston cohort as PRJEB22893 (public).

**Spencer 2021** — Science 374:1632-1640, "Dietary fiber and probiotics
influence the gut microbiome and melanoma immunotherapy response",
doi:10.1126/science.aaz7015.
- Data availability (verbatim): "Raw sequencing data and all relevant human
  data necessary for reproducing results are available in the NCBI Sequence
  Read Archive under BioProject ID PRJNA770295." https://europepmc.org/article/PMC/PMC8970537
- NCBI BioProject: MD Anderson, "Raw sequence reads", 807 SRA experiments,
  0.43 TB: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA770295
  ENA split: 309 WGS runs (557.4 GB; HiSeq 2000 / HiSeq 3000 / HiSeq X) and
  498 amplicon runs (7.5 GB): https://www.ebi.ac.uk/ena/browser/view/PRJNA770295
- Design: observational study of diet and probiotic use in patients on ICB
  (128 patients in the fiber analysis; 132 in the new anti-PD-1 cohort).
  FMT appears only in germ-free-mouse experiments and in re-analysis of two
  *other* groups' published FMT trials. Whether any patients got combination
  ICB regimens is **not confirmed** from the text read.

**Tsakmaklis 2023** — *BMC Cancer* 23:1160, "TIGIT+ NK cells in combination with
specific gut microbiota features predict response to immune checkpoint
inhibitor therapy in melanoma patients", doi:10.1186/s12885-023-11551-5.
- Data availability (verbatim): "The datasets generated and analyzed during the
  current study are available in the Sequence Read Archive (SRA) repository
  (NCBI), [PRJNA1011235]." https://europepmc.org/article/PMC/PMC10685659
- NCBI BioProject: University Hospital Frankfurt (Goethe University),
  title "ANIMA", 29 SRA experiments, 79,237 Mbytes:
  https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1011235
  ENA: 29 WGS runs, 92.8 GB, NovaSeq 6000, first public 2023-09-02:
  https://www.ebi.ac.uk/ena/browser/view/PRJNA1011235
- Design: "prospective, non-interventional" study, 29 cutaneous melanoma
  patients; regimens anti-PD-1 15 (51.7%), anti-PD-1 + anti-CTLA-4 13
  (44.8%), anti-CTLA-4 1 (3.5%); no FMT mentioned.
  Same article link.
- Discrepancy: the BioProject description says "40 melanoma patients" but
  the project holds 29 runs and the paper analyses 29 patients.

### Gopalakrishnan 2018 — is the patient-level shotgun data public?

**Yes. The shotgun data are public on ENA. EGA holds exome data only.**
- Science 359:97-103, "Gut microbiome modulates response to anti-PD-1
  immunotherapy in melanoma patients", doi:10.1126/science.aan4236. Data
  statement (verbatim): "Fecal, oral and murine 16S, and fecal WGS data are
  available from the European Nucelotide Archive under accession numbers
  PRJEB22894, PRJEB22874, PRJEB22895 and PRJEB22893 respectively. Human WES
  data are available from the European Genome-phenome Archive under accession
  number EGAS00001002698." https://europepmc.org/article/PMC/PMC5827966
- So the accessions map as: PRJEB22893 = fecal **WGS (shotgun)**; PRJEB22894 =
  fecal 16S; PRJEB22874 = oral 16S; PRJEB22895 = murine 16S.
- ENA, all first public 2017-11-06
  (https://www.ebi.ac.uk/ena/browser/view/PRJEB22893):
  - PRJEB22893: **25 paired WGS runs, 56.5 GB, HiSeq 2000**, FASTQ FTP links
    present for every run.
  - PRJEB22894 / PRJEB22874 / PRJEB22895: **one amplicon run each** (0.09 /
    0.13 / 0.01 GB, MiSeq). The paper reports 43 fecal 16S samples, so one
    run per project is odd; this was not resolved (possibly pooled/barcoded
    submissions). Not shotgun, so it does not affect the question asked.
- EGA was not checked directly; this reports what the paper's statement says
  for EGAS00001002698 (exome, not microbiome).
- The same-lab public record is also in McCulloch's statement, which lists
  PRJEB22893 as the public "Houston" cohort:
  https://europepmc.org/article/PMC/PMC10246505

### Something that affects the project: the compilation's "Gopalakrishnan_2018" is in Spencer's BioProject

All 22 `Gopalakrishnan_2018` runs in `metadata/response_sheet.xlsx` resolve in
ENA to **PRJNA770295, the Spencer 2021 project** (HiSeq 2000, 46.4 GB, sample
names like `W21.203`), not to Gopalakrishnan's own PRJEB22893. See
`ena_per_run_response_sheet_studies.tsv`. So the 22 runs are public, but
whether they are the same patients as the 25 runs in PRJEB22893 **cannot be
checked from anything available** (no patient-ID crosswalk). Also relevant:
Spencer's paper excluded patients from the earlier Gopalakrishnan cohort from
its new-cohort analysis, so overlap between "Gopalakrishnan_2018" (22) and
"Spencer_2021" (134) in the compilation is not ruled out either. Worth
resolving before any future pooling that includes those studies.

### What could not be confirmed
- Liu 2022: which immunotherapy agents the 14 patients received, and whether
  they had concurrent chemotherapy.
- Spencer 2021: whether any patients received combination ICB.
- Why only 11 / 27 / 134 / 22 runs of Heshiki / McCulloch / Spencer /
  Gopalakrishnan are in the compilation (not documented in the sources read).
- EGA record contents (not queried).
