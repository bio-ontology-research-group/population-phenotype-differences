# Population differences in Mendelian disease phenotypes: analysis code and data

Reproducibility package for the case-to-case analysis of population-specific
phenotype differences in published and clinical Mendelian disease cases
(Althagafi, Abdelhakim, Schofield, Hoehndorf; manuscript in preparation).
Cases are compared with each other, within gene, never with a disease profile
from OMIM or HPOA. Population differences are separated from provenance
(curation corpus, reporting group, author region) by permuting between
reporting groups within gene x corpus strata.

`RESULTS.md` is the running analysis log, with every number reported and the
command that produced it. `docs/feasibility_2026-10-24.md` is the initial
feasibility assessment on PAVS and the GA4GH Phenopacket Store.

## Data

### Phenopackets (`results/phenopackets/`)

All 66,198 patients analysed, written as GA4GH Phenopacket v2 JSON
(`to_phenopackets.py`; validated with `validate_phenopackets.py`):

| Directory | Cases | Source |
|---|---|---|
| `store/` | 9,517 | GA4GH Phenopacket Store v0.1.26, with patient origin re-extracted from the source articles |
| `pmcp/` | 6,792 | PMC-Patients v2 case reports, Qwen extraction, RAG-HPO-style HPO mapping |
| `pubmed_abs/` | 45,911 | PubMed case-report abstracts, Qwen extraction |
| `pavs_clin/` | 2,576 | PAVS Saudi clinical cases |
| `pavs_lit/` | 1,278 | PAVS curated Saudi literature cases |
| `pavs_us/` | 124 | PAVS source Q (Ziats et al. 2020), a US cohort |

`results/phenopackets/index.tsv` lists every packet with corpus, PMID, gene,
population, ancestry x place-of-care setting, sex, age, origin and place of
care. Population and provenance labels are subject extensions
(`populationRegion`, `populationRegionSource`, `reportingAuthorRegion`,
`curationCorpus`, `reportingGroup`, `ancestryCareSetting`, `geographicOrigin`,
`ancestry`, `placeOfCare`). Excluded (negated) findings carry `excluded: true`.

Article text is not redistributed; phenotype descriptions are the short
extracted phrases only.

### Case tables and results

- `results/cases_*.tsv`: one row per case (id, corpus, pmid, gene, HPO terms,
  population, author region, reporting group); `results/v4/fine/` splits the
  OTHER bin into SASIA, LATAM and AFR, `results/v4/fine2/` adds age band and
  sex, `results/v4/setting2/` adds ancestry x place-of-care labels.
- `results/termclass.jsonl`: Qwen classification of each HPO term into
  laboratory, imaging, physiological test, pathology, clinical examination or
  symptom/history.
- `results/groups_lastauthor.tsv`, `results/all_authors.tsv`: reporting groups.
- `results/v4/`: all v4 outputs (all-pairs similarity contrasts, gene recovery,
  term and outcome tests per population, organ system, age and sex, World Bank
  covariates, figures) and `results/v4/qvalues.tsv` (Benjamini-Hochberg
  q-values per test family).

### Not in git

- HPO release 2026-01-08: `data/get_hpo.sh` downloads it into `data/hpo/`.
- The OntoFUNC term database (`term.txt`, `term2term.txt`, `graph_path.txt`)
  built from that release with ELK; expected in `data/ontofunc_hp/`.
- Similarity caches (`*.npy`, up to 17 GB); recomputed by `case_pairs.py`.
- Raw LLM responses that quote source text, and the open-access lookup cache.

## Pipeline

1. **Corpora.** `store_origin_extract.py`, `tag_store_pmids_country.py`
   (Phenopacket Store origins), `pmcp_extract.py` + `pmcp_to_cases.py`
   (PMC-Patients), `pubmed_abstracts.py` + `abs_to_cases.py` (PubMed
   abstracts), `tdm_fetch.py` (publisher TDM full text; needs institutional
   keys in an env file, never in the repository), `hpo_map2.py` +
   `qwen_tasks.py` (HPO normalisation, term classes, ancestry and care
   setting), `validate_extraction.py` (extraction audit against Store
   curation). LLM: Qwen3.8-27B-FP8 served with vLLM; job scripts in `jobs/`
   and `jobs/orix/`.
2. **Labels.** `fine_regions.py` (six populations), `add_age_sex.py`,
   `setting_labels.py`, `to_phenopackets.py`.
3. **Similarity decomposition.** `case_pairs.py` (symmetric best-match-average
   Resnik with corpus information content, hubness z-normalisation) and
   `population_model.py pairs|allpairs [--vs-rest]` (gene fixed effects,
   group-cluster bootstrap, total and patient-carried permutations).
4. **Gene recovery.** `recovery_main.py`, `retrieval_perm.py`.
5. **Term and outcome tests.** `population_model.py terms` (Mantel-Haenszel
   per term, organ-system branches, test-based vs bedside outcomes;
   `--where`, `--key-extra age_band|sex`).
6. **Country covariates.** `healthspend.py` (World Bank health expenditure,
   GDP, physician density; country-level permutation).
7. **Multiple testing and figures.** `qvalues.py`, `network_figs.py`.

Earlier stages of the project (PAVS benchmark re-analysis, Store-only origin
models) are `run_analyses.py`, `followup_analyses.py`, `rank_engine.py`,
`origin_analyses.py`, `pooled_origin_model.py`, `merged_model.py`,
`phenotype_specificity.py`.

The exact commands for each result are in `RESULTS.md` and `jobs/`
(`terms_tasks.txt`, `setting_tasks.txt` hold the array-task arguments). Job
scripts contain the cluster paths they were run with; adjust the working
directory before reuse.

## Licence

Code (`*.py`, `jobs/`, `data/*.sh`): BSD 3-Clause, see `LICENSE`.
Data (`results/`, including the phenopackets): CC BY 4.0, see `LICENSE-DATA`.
Records derived from PMC-Patients (`results/phenopackets/pmcp/`,
`results/cases_pmcp*.tsv`, `results/pmcp/`) are also subject to the
PMC-Patients licence (CC BY-NC-SA 4.0).

## Software

Python 3.10 or later; `pip install -r requirements.txt`. Statistics were run
with numpy 1.26 / scipy 1.11 (unimatrix) and numpy 2.2 / scipy 1.15 (local);
results agree to floating-point precision. Figures need graphviz (`twopi`).
