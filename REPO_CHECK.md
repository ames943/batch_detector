# REPO_CHECK

Checks run on this folder before it is pushed. Nothing has been pushed and no
remote is configured. Size: **56.9 MB** of working files (about
61.4 MB with `.git`).

## 1. Tree

```
batch_detector_repo/
├── .gitignore
├── .zenodo.json
├── CITATION.cff
├── LICENSE
├── README.md
├── REPO_CHECK.md
├── batch_detector.py
├── examples/   (4 files, 0.0 MB)
│   ├── example_features.tsv
│   ├── example_metadata.tsv
│   ├── make_example_data.py
│   ├── run_example.sh
├── paper/   (223 files, 52.3 MB)
│   ├── README.md
│   ├── data/   (2 files, 0.0 MB)
│   ├── metadata/   (5 files, 0.0 MB)
│   ├── requirements.txt
│   ├── results/   (152 files, 51.7 MB)
│   ├── scripts/   (62 files, 0.6 MB)
├── requirements-dev.txt
├── requirements.txt
├── tests/   (1 files, 0.0 MB)
│   ├── test_batch_detector.py
```

`paper/` keeps the original project layout (`scripts/`, `metadata/`, `results/`)
because the scripts use paths relative to the project root. Input tables
therefore sit under `paper/results/ml/` and `paper/metadata/`; `paper/data/`
holds the data-source notes and the download helper.

## 2. batch_detector.py

Changes relative to the original v2.0.0 tool:

- Comments, docstrings and log wording rewritten; type annotations removed;
  `__version__ = "2.0.0"`, a `--version` flag and `main(argv=None)` added.
- `_pooling_recommendation_v2` docstring now says the 0.25 and 4 cutoffs are
  conservative rules of thumb and that the calibrated simulation found no clean
  ratio-based threshold. The same wrong claim was in the text the tool prints
  and saves (`recommendation_reason`) and in the `--legacy-recommendation`
  help; those now say "heuristic cutoff". Decision logic and numbers are unchanged.
- Figure, right panel in the default mode: the message now reads
  "Power curve drawn only with --legacy-power. / Default power grid:
  batch_detector_power_grid.tsv" (string only; the old text said
  "Power analysis unavailable (non-binary response)").
- The second loop of `_find_n_80pct` was removed; it could never run.
- Defaults are unchanged. No hard-coded paths (a test checks this).

Checks:

1. Old and new file compared as syntax trees (docstrings and string literals
   blanked): 27 of 30 functions identical; the other three are the parser
   (`--version`), `_find_n_80pct` (dead loop) and `main` (argv, log lines).
2. Old and new run on the example data in default mode after the figure-text
   change: all 36 JSON fields (except the timestamp) and the power table are
   identical. Earlier, the same comparison was also run in legacy, single-cohort
   and multiplicative-replacement modes (identical).
3. New file on the paper's n = 283 data, clean Python 3.9 virtualenv, pinned
   versions: all 33 fields of `batch_detector_report.json`, the recommendation
   text and the power grid match `paper/results/revision/task9_n283/`
   (198 s). That run was done before the figure-text change, which only affects
   the PDF.
4. The rendered figure was inspected: the new message fits inside panel B.

## 3. Data: licenses and what is included

Sources and wording are in `paper/data/DATA_SOURCES.md`.

| input | file(s) and size | original source | safe to redistribute? | included |
|---|---|---|---|---|
| Genus-level count/CLR tables, 4 cohorts | `results/ml/{n283_4cohort,lee2022,n118_3cohort}/X_genus_{raw,clr}.tsv` (9.1 + 3.7, 5.0 + 2.0, 3.4 + 1.4 MB); phylum/species CLR tables for n = 118 (0.15 and 19.8 MB) | derived from public INSDC reads: PRJNA397906, PRJNA615114, PRJNA399742, PRJEB43119 | Yes. INSDC places no restrictions on use of deposited data; the tables contain no reads. The Peng 2020 article's own data statement is paywalled and was not read. | **yes** (decision: keep) |
| Orletskaia response labels | `metadata/response_labels*.tsv`, `response_labels_n283.tsv`, `lee2022_labels.tsv` (3 to 6 KB; run, response, cohort only) | Orletskaia & Olekhnovich 2026 (CC BY 4.0); authors' repo JeniaOle13/cancer-biomarkers (MIT) | Yes, with attribution. `response_sheet.xlsx` (with a `log_ratio` column) could not be matched to any file the authors publish. | subset yes; `response_sheet.xlsx` **no** |
| Lee 2022 site metadata | `metadata/lee2022_sites.tsv` (3 KB) | ENA sample attribute `Cohort` of PRJEB43119 (matches the file for 165 of 165 runs) | Yes, with attribution | yes |
| cBioPortal tumor tables | `combined_tumor_features.tsv` and per-study tables (0.1 MB); raw JSON (95 MB) | cBioPortal: Riaz 2017, Liu 2019, Hugo 2016 | Not confirmed: ODbL share-alike; no license file found for the Riaz study | **no**; rebuilt by `cbioportal_download.py` + `process_cbioportal.py` |
| Duvallet CRC data | `input_feature_matrix.tsv` (0.55 MB), `input_labels.tsv` | microbiomeHD (Zenodo 10.5281/zenodo.840333) | No: CC BY-NC 4.0, and the GitHub repo has no license | **no**; `data/get_duvallet_crc.sh` rebuilds them (byte-identical to the paper's files, checked by md5) |
| Kraken2 reports, fastp JSON | 173 MB and 168 MB | derived | would be fine; left out for size | no |
| Raw FASTQ | about 1 TB | INSDC | public | no |

`paper/data/DATA_SOURCES.md` lists the four accessions with their original
papers and carries the line "Genus tables were derived by the author from
public INSDC reads; please cite the original studies if you use them."

## 4. Tests and example

`tests/test_batch_detector.py`, 17 tests (E0 formula, chance level against
shuffled labels, PERMANOVA on separated groups, reproducible example data,
cohort effect with p < 0.05, response excess about 0, NO-GO call,
`--zero-handling multiplicative_replacement`, single-cohort mode, the GO /
CAUTION / NO-GO branches, `--help` and `--version`, no hard-coded paths).

- Clean Python 3.9 virtualenv, pinned versions: `17 passed, 14 warnings in 5.79s`
- Python 3.13, newer libraries (numpy 2.1, pandas 2.2, scipy 1.17, matplotlib
  3.10): `17 passed in 7.13s`
- The warnings in the 3.9 run come from matplotlib's own dependencies.

`bash examples/run_example.sh`, from a fresh copy in the clean virtualenv
(exit 0, 15 s):

```
Response PERMANOVA : R2 = 0.0055   p = 0.660 (ns)
Batch PERMANOVA    : R2 = 0.4077   p = 0.001 (**)
Response delta_R2  : -0.00289  (E0=0.00840)
Batch delta_R2     : +0.39093  (E0=0.01681)
Batch:signal ratio : 73.99x
rho null percentile: 100.0
Recommendation     : NO-GO  (v2_chance_calibrated)
Min n (80% power)  : 240 same-protocol samples
```

Only `--n-grid 60,120,240` differs from the defaults (it leaves out n = 480);
the full default run takes about 40 s.

## 5. Does paper/ reproduce the paper?

Each script was run once from a fresh copy (clean Python 3.9 environment,
pinned versions, 90 s limit); details are in `paper/README.md`.

- 29 scripts finished (`phase2_simulation.py` only printed its usage because it
  needs `--mode`), 9 are long jobs that started fine, 15 stopped at once: 14
  because they need data that is not shipped (Kraken2 reports, tumor and CRC
  tables, raw cBioPortal JSON) and 1 because it needs `--dataset`.
- 101 regenerated TSV/JSON files are identical to the shipped ones.
- Elastic-net results differ in the third decimal on rerun (Table X: up to
  0.006 AUC and 0.016 permutation p, for example C2 to C1 0.550 shipped and
  0.544 rerun; n = 39 LOOCV AUC 0.447 and 0.450). Random-forest rows are
  identical. The cause was not found; no conclusion changes.
- `steps_1_4_c4_integration.py`, `review_task4*.py` and
  `review2_task2_rarefaction_5M.py` overwrite their output tables with empty
  ones when the Kraken2 reports are missing (warned about in `paper/README.md`;
  scripts not changed).
- Run order matters for `bootstrap_cis.tsv`, `permdisp_results.tsv` and
  `permanova_comparison.tsv` (documented in `paper/README.md`).
- Two items had no script: `permanova_n281_dropped2.tsv` (now reproduced exactly
  by the new `review2_task3_drop_two_samples.py`) and
  `task7_cross_cohort_full_table.tsv` (a hand-made reformatting of
  `holdout_summary.tsv`).
- Not run: the multi-hour jobs (`task5`, `phase2 --mode full`, nested CV to
  completion) and everything that needs the Kraken2 reports.

## 6. Privacy grep

Run over every file except this one (text files by line; the 16 PDF/PNG files
were searched as binary and have no hits for the username, email address, home
path or "Claude"). Patterns: `/Users/`, `/home/`, `ameygarg52`, `token`, `key`,
`password`.

| pattern | lines | files | what they are |
|---|---|---|---|
| `/Users/` | 1 | 1 | `tests/test_batch_detector.py:152`, a string literal in the test that checks the tool has no hard-coded paths |
| `/home/` | 2 | 2 | the same test line, and `paper/data/DATA_SOURCES.md:12`, part of the NCBI URL `.../home/about/policies/` |
| `ameygarg52` | 0 | 0 | none |
| `password` | 0 | 0 | none |
| `token` | 4 | 1 | `paper/scripts/build_matrix.py` lines 16, 28, 29, 54: comments about the "last token" of a genus name |
| `key` | 78 | 40 | all benign: dictionary idioms (`.keys()`, `sort_keys`, `dict.fromkeys`, variables `dkey`, `mkey`, `key`), the word "keyword(s)", and bacterial genus names in the abundance tables (Dickeya, Starkeya, Hawkeyevirus, Monkeypox virus) |

No API keys, tokens, `.env` files, credentials, usernames or home-directory
paths were found. A further search of all text files for "Claude", "Anthropic",
"assistant", "per instructions", "project memory" and "your suspicion" found
nothing apart from the reword log below.

## 7. Reworded lines in the copied notes and scripts

Wording only; no number, file name or result was changed. 37 changes
(" ⏎ " marks a line break inside the quoted text). Row 33 is `paper/README.md`
(rewritten from first person); the last three rows are comments or a printed
message in figure scripts.

| # | file | old | new |
|---|---|---|---|
| 1 | `paper/results/revision/SUMMARY.md` | not melanoma; see project ⏎ memory), consistent | not melanoma; see project ⏎ notes), consistent |
| 2 | `paper/results/revision/SUMMARY.md` | Per your instruction to match "the same prevalence | The aim was to match "the same prevalence |
| 3 | `paper/results/revision/SUMMARY.md` | lookup on `SRR11413606` this session, | lookup on `SRR11413606` during this review, |
| 4 | `paper/results/revision/SUMMARY.md` | One update worth acting on: Orletskaia & | One update: Orletskaia & |
| 5 | `paper/results/revision/SUMMARY.md` | Journal* (doi:10.34133/csbj.0065) — consider citing that version instead of, ⏎ or alongside, the preprint. | Journal* (doi:10.34133/csbj.0065); that version can be cited instead of, ⏎ or alongside, the preprint. |
| 6 | `paper/results/revision/review/REVIEW_TASKS.md` | (verified against project metadata files ⏎ this session, not just recalled from memory): | (verified against project metadata files ⏎ during this review, not just taken from earlier notes): |
| 7 | `paper/results/revision/review/REVIEW_TASKS.md` | I queried the ENA API directly this session (`api.crossref.org` is for ⏎ papers; | The ENA API was queried directly during this review (`api.crossref.org` is for ⏎ papers; |
| 8 | `paper/results/revision/review/REVIEW_TASKS.md` | If you want this pursued further, the concrete next ⏎ step is: | If this is pursued further, the concrete next ⏎ step is: |
| 9 | `paper/results/revision/review/REVIEW_TASKS.md` | The machine used throughout this session — including | The machine used throughout this work — including |
| 10 | `paper/results/revision/review/REVIEW_TASKS.md` | queried directly this session: | queried directly during this review: |
| 11 | `paper/results/revision/review/REVIEW_TASKS.md` | Two ⏎ requested changes, plus one fix found while checking the PNG: | Two ⏎ changes, plus one fix found while checking the PNG: |
| 12 | `paper/results/revision/figures/FIGURES.md` | **Note per instructions:** | **Note:** |
| 13 | `paper/results/revision/figures/FIGURES.md` | **percentile_norm excluded per instructions**: | **percentile_norm excluded as planned**: |
| 14 | `paper/results/revision/figures/FIGURES.md` | Two layout fixes, requested after review of v2: | Two layout fixes, made after review of v2: |
| 15 | `paper/results/revision/figures/FIGURES.md` | Two requested changes: (1) axes relabeled | Two changes: (1) axes relabeled |
| 16 | `paper/results/revision/final/FINAL_CHECKS.md` | No code changed (per instructions) — this is a | No code was changed — this is a |
| 17 | `paper/results/revision/final/FINAL_CHECKS.md` | Your suspicion is confirmed: **0.10 | Confirmed: **0.10 |
| 18 | `paper/results/revision/final/FINAL_CHECKS.md` | **No.** Recommend replacing the paper's sentence. Two honest options: | **No.** The paper's sentence should be replaced. Two options: |
| 19 | `paper/results/revision/final/FINAL_CHECKS.md` | - Purely honest, no false precision: | - Plain version, no false precision: |
| 20 | `paper/results/revision/final/FINAL_CHECKS.md` | - Slightly shorter variant matching your suggested phrasing: | - Slightly shorter variant: |
| 21 | `paper/results/revision/final/FINAL_CHECKS.md` | Only actionable ⏎ item: consider citing the now-published CSBJ version of Orletskaia & ⏎ Olekhnovich alongside or instead of the bioRxiv preprint. | Only actionable ⏎ item: the now-published CSBJ version of Orletskaia & ⏎ Olekhnovich can be cited alongside or instead of the bioRxiv preprint. |
| 22 | `paper/results/revision/review2/REVIEW2_TASKS.md` | (agents not stated in the text I read) | (agents not stated in the text read) |
| 23 | `paper/results/revision/review2/REVIEW2_TASKS.md` | the six studies asked about, and I did not check EGA itself. | the six studies asked about, and EGA itself was not checked. |
| 24 | `paper/results/revision/review2/REVIEW2_TASKS.md` | stated in anything I could access. | stated in anything accessible. |
| 25 | `paper/results/revision/review2/REVIEW2_TASKS.md` | in the text I extracted, so I can't confirm either.** The paper's second | in the extracted text, so neither can be confirmed.** The paper's second |
| 26 | `paper/results/revision/review2/REVIEW2_TASKS.md` | other public SRA projects that I did ⏎   not identify; | other public SRA projects that were ⏎   not identified; |
| 27 | `paper/results/revision/review2/REVIEW2_TASKS.md` | is **not confirmed** from the text I read. | is **not confirmed** from the text read. |
| 28 | `paper/results/revision/review2/REVIEW2_TASKS.md` | run per project is odd; I did not resolve why (possibly | run per project is odd; this was not resolved (possibly |
| 29 | `paper/results/revision/review2/REVIEW2_TASKS.md` | - I did not check EGA directly; I'm reporting what the paper's statement says | - EGA was not checked directly; this reports what the paper's statement says |
| 30 | `paper/results/revision/review2/REVIEW2_TASKS.md` | checked from anything I have** (no patient-ID crosswalk). | checked from anything available** (no patient-ID crosswalk). |
| 31 | `paper/results/revision/review2/REVIEW2_TASKS.md` | ### What I could not confirm | ### What could not be confirmed |
| 32 | `paper/results/revision/review2/REVIEW2_TASKS.md` | (not documented in the sources I read). | (not documented in the sources read). |
| 33 | `paper/README.md` | I compared the files regenerated by those runs with the shipped ones. | The files regenerated by those runs were compared with the shipped ones. |
| 34 | `paper/results/revision/SUMMARY.md` | call is "there IS real signal (p=0.001) but pool at your own risk" (rho | call is "there IS real signal (p=0.001) but pooling is risky" (rho |
| 35 | `paper/scripts/revision/figures/fig7_simulation_grid_v3.py` | artifact, now shown rather than omitted, per instructions. | artifact, now shown rather than omitted. |
| 36 | `paper/scripts/revision/figures/fig4_power_curves.py` | were run) -- per instructions, only the single (n,power) marker for that bound is | were run) -- only the single (n,power) marker for that bound is |
| 37 | `paper/scripts/revision/figures/fig4_power_curves.py` | marker is drawn, per instructions." | marker is drawn." |

## 8. Left out, and open items

- Tumor tables, Duvallet data and `response_sheet.xlsx` are left out (section
  3); `task1_*` and `task2_*` cannot run until the tumor and CRC tables are
  rebuilt, and `review2_task4_ena_lookup.py` needs the label sheet.
- Placeholders still to fill: `GITHUB_USERNAME` in `CITATION.cff`;
  `arXiv:XXXX.XXXXX` in `README.md` and `paper/README.md`; the Zenodo DOI,
  which the README says will be added after the first release.
- Commit author is `Amey Garg <noreply@example.com>` until the steps below
  are run.
- Not copied from the original project: the HUMAnN3 scripts, the first drafts
  of the LOOCV, meta-analysis and early figure code, the TCGA scripts, logs, and
  the raw and trimmed read directories.

## 9. Commands to run (not run here)

Fill in the placeholders and set the author. Replace `YOURNAME` with the
GitHub username and `GITHUB_NOREPLY_EMAIL` with the address from
GitHub > Settings > Emails (it looks like `12345678+YOURNAME@users.noreply.github.com`).
These are for macOS (`sed -i ''`).

```bash
cd ~/Downloads/batch_detector_repo

# 1. fill in the placeholders (add the arXiv ID once it exists)
sed -i '' 's/GITHUB_USERNAME/YOURNAME/' CITATION.cff
# sed -i '' 's/arXiv:XXXX.XXXXX/arXiv:2610.01234/' README.md paper/README.md

# 2. set the author and fold the edits into the single commit
git config user.name  "Amey Garg"
git config user.email "GITHUB_NOREPLY_EMAIL"
git add -A
git commit --amend --reset-author --no-edit
git log -1 --format='%an <%ae>%n%s'      # check author and message

# 3. create an EMPTY repo named batch_detector on GitHub (no README, no license),
#    then connect it and push
git remote add origin git@github.com:YOURNAME/batch_detector.git
#   (or over HTTPS: https://github.com/YOURNAME/batch_detector.git)
git branch -M main
git push -u origin main
```

For the Zenodo archive: sign in to zenodo.org with GitHub, switch on the
`batch_detector` repository under GitHub settings there, then publish a GitHub
release tagged `v2.0.0`. Zenodo reads `.zenodo.json` and creates the DOI;
afterwards add the DOI to `README.md` and `CITATION.cff` and commit.
