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
- `results/v4/`: all v4 outputs and `results/v4/qvalues.tsv` (Benjamini-Hochberg
  q-values per test family, `qvalues.py`). Directories, the script that fills
  each, and where the paper uses them:

  | Directory | Script (job) | Paper |
  |---|---|---|
  | `allpairs/` | `population_model.py allpairs` (`jobs/allpairs_um.sbatch`); `collect_tables.py` builds `contrasts.tsv` | Figure 4, Table S11, Figure S5 |
  | `within/` | `population_model.py allpairs --vs-rest` (`jobs/vsrest_um.sbatch`) | Table S12 |
  | `r2/` | `population_model.py pairs --r2-only` (`jobs/r2_um.sbatch`); `r2_summary.py` builds `summary.csv` | Table S17 (variance decomposition) |
  | `recovery/` | `recovery_main.py` (`jobs/rec_um.sbatch`) | Tables S14 and S15, Figure S6 |
  | `terms_all/` | `population_model.py terms` (`jobs/terms_um.sbatch`, `terms_tasks.txt`); `collect_tables.py` builds `all_outcomes.tsv` | Figures 1 and 2, Tables 1, S6, S7, Figures S1 and S3 |
  | `setting_all/` | `population_model.py terms` on setting labels (`jobs/setting_um.sbatch`, `setting_tasks.txt`); `collect_tables.py` builds `test_share_q.tsv` | Figure 3a, Table S9 |
  | `setting_pooled/` | `setting_pooled.py` | Table S9 (pooled rows) |
  | `worldbank/` | `healthspend.py` | Figure 3b and c, Table S10 |
  | `gda/` | `gda_impact.py rank` (`jobs/gda_um.sbatch`) and `gda_impact.py stats`; `collect_tables.py` adds q | Table 1, Table S13 |
  | `gda_rand/` | `gda_impact.py rank --random-control` (`jobs/gda_rand_um.sbatch`); `random_control.py` | Table S13 (continued) |
  | `../v5/tool_restriction/` | `export_profiles.py` (`jobs/tr_export_um.sbatch`); `tool_restriction_run.py` (`jobs/tool_restriction_um.sbatch`, configuration in `jobs/tool_restriction/`); `tool_restriction_universe.py`; `tool_restriction.py` | Table S13 (continued, Exomiser and LIRICAL) |
  | `terms_noclass/` | `population_model.py terms --require-classified` (`jobs/zeroshare_um.sbatch`, `zeroshare_tasks.txt`); `zero_share_sensitivity.py` | Table S6 (continued) |
  | `symptom_sign/` | `symptom_sign.py` | Table S5 |
  | `age_sex/` | `age_sex_interaction.py` | Table S7, Figure S2 |
  | `comorbidity/` | `comorbidity_net.py` | Table S8, Figure S4 |
  | `validation/` | `validate_against_store.py`, `corpus_replication.py` | Table S2 (c to e) |
  | `audit_abs/`, `audit_pmcp/` | language-model judgements of the extraction audits (`docs/extraction_audit.md`) | Table S2 (a, b) |
  | `sensitivity/` | `sensitivity_pakistan.py` | Table S16 |
  | `figures/` | `gim_figures.py`, `network_figs.py` | all figures |
  | `summaries.txt` | `summaries.py` | descriptive counts |

  Files at the top of `results/v4/` other than `qvalues.tsv` (`terms_abs_*`,
  `terms_age_*`, `terms_sex_*`, `terms_setting_*`, `terms_all_*`,
  `cases_*_age_*`, `cases_*_sex_*`, `cases_*_setting.tsv`,
  `recovery_v3data_*`) are superseded intermediate runs kept for traceability.

### Not in git: inputs to obtain

`jobs/run_local.sh` expects them under `$DATA` (default `data/`); the cluster
job scripts in `jobs/` used the same layout under the job directory.

- HPO release 2026-01-08 (`hp.obo`, `hp.owl`, `phenotype.hpoa`,
  `genes_to_phenotype.txt`): `data/get_hpo.sh` downloads it into `$DATA/hpo/`.
- The OntoFUNC term database (`term.txt`, `term2term.txt`, `graph_path.txt`) in
  `$DATA/ontofunc_hp/`: built from `hp.owl` of that release by classifying it
  with ELK and exporting the inferred is-a closure in the OntoFUNC/GO-style
  table format (the Bio-Ontology Research Group's standard pipeline).
- GA4GH Phenopacket Store release 0.1.26
  (https://github.com/monarch-initiative/phenopacket-store/releases/tag/0.1.26,
  `all_phenopackets.zip`) in `$DATA/ppktstore/0.1.26/`.
- PMC-Patients V2 (https://huggingface.co/datasets/zhengyun21/PMC-Patients,
  CC BY-NC-SA 4.0) for `pmcp_extract.py`.
- PubMed case-report abstracts: `pubmed_abstracts.py` queries E-utilities and
  writes `$DATA/abs/abstracts.jsonl` (titles and abstracts; not redistributed).
  The same script fetches the abstracts of the Store publications
  (`$DATA/storeabs/cands.jsonl`) for the validation.
- PAVS phenopackets (`PAVS_phenopackets.json`) from the PAVS database export
  (https://pavs.borg.kaust.edu.sa); `to_phenopackets.py --pavs-json`.
- HGNC complete set (`hgnc_complete_set.txt`,
  https://www.genenames.org/download/archive/) in `$DATA/hgnc/`.
- BGE-small embedding model (`BAAI/bge-small-en-v1.5`, downloaded by
  sentence-transformers on first use) for `hpo_map2.py`.
- World Bank indicators (SH.XPD.CHEX.PC.CD, NY.GDP.PCAP.CD, SH.MED.PHYS.ZS):
  `healthspend.py` downloads them from https://api.worldbank.org/v2 into
  `results/v4/worldbank/`.
- Similarity caches (`*.npy`, up to 17 GB); recomputed by `case_pairs.py`, which
  stores the case order beside the cache and refuses a cache built for another
  case set.
- Natural Earth 1:110m country shapes (`ne_110m_admin_0_countries.shp`, public
  domain, https://www.naturalearthdata.com) in `$DATA/naturalearth/` for the
  country map (`gim_figures.py --shapes`).
- Raw LLM responses that quote source text (the tracked origin and setting files
  were stripped of their evidence fields with `tools/strip_evidence.py`), and
  the open-access lookup cache.

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
7. **Gene identification against annotations.** `gda_impact.py rank` (cluster
   chunks) and `gda_impact.py stats` (leave-publication-out rank, coverage,
   bedside-only and test-only rankings).
8. **Signs versus symptoms, age and sex, pooled place of care.**
   `symptom_sign.py`, `age_sex_interaction.py`, `setting_pooled.py`.
9. **Organ-system co-reporting networks.** `comorbidity_net.py` (`--split` for
   test-based versus bedside findings).
10. **Extraction performance and replication.** `validate_against_store.py`
    (exact and lenient precision and recall against Store curation),
    `corpus_replication.py` (test-share gap per corpus and abstract-versus-full-text
    bias); the language-model audits are described in `docs/extraction_audit.md`.
    `sensitivity_pakistan.py` repeats the main contrast with Pakistan and
    Afghanistan in South Asia.
11. **Multiple testing, tables and figures.** `collect_tables.py` (contrast,
    outcome, setting and gene-ranking tables from the cluster outputs),
    `qvalues.py` (Benjamini-Hochberg q per family; the q-values of the
    symptom_sign, age_sex, setting_pooled, gda, comorbidity and validation
    families are computed inside those scripts), `summaries.py` (descriptive
    counts), `network_figs.py` (`shift`, `dag`, `spend`, `popnet`;
    `popnet --print-layout` reports the node layout), `gim_figures.py` (the three main figures, including corpus replication
    and paired bedside-restriction ranking, and supplementary plots), and
    `network_figs.py popnet` / `dag` (the supplementary population network
    and HPO subgraphs). The main ranking figure reads the archived paired
    estimates and intervals in `results/v4/gda/gda_stats.csv`; no ranking
    computations are repeated when plotting.

Earlier stages of the project (PAVS benchmark re-analysis, Store-only origin
models) are `run_analyses.py`, `followup_analyses.py`, `rank_engine.py`,
`origin_analyses.py`, `pooled_origin_model.py`, `merged_model.py`,
`phenotype_specificity.py`.

`jobs/run_local.sh` records every local stage with the exact arguments, in
order; the cluster stages are the sbatch files in `jobs/` (`terms_tasks.txt`
and `setting_tasks.txt` hold the array-task arguments) and `jobs/orix/` holds
the vLLM extraction jobs. Job scripts contain the cluster paths they were run
with; adjust the working directory before reuse. `RESULTS.md` is the running
log of every result with the command that produced it. Language-model
extraction was run with temperature 0, but GPU inference is not
bit-reproducible, so re-running the extraction can change individual phrases.

## Licence

Code (`*.py`, `jobs/`, `data/*.sh`): BSD 3-Clause, see `LICENSE`.
Data (`results/`, including the phenopackets): CC BY 4.0, see `LICENSE-DATA`.
Records derived from PMC-Patients (`results/phenopackets/pmcp/`,
`results/cases_pmcp*.tsv`, `results/pmcp/`) are also subject to the
PMC-Patients licence (CC BY-NC-SA 4.0).

## Software

System tools: `pdftotext` (poppler-utils; `store_origin_extract.py`,
`tdm_fetch.py`), graphviz `twopi` (`network_figs.py dag`), vLLM 0.28.0 with
`Qwen/Qwen3.8-27B-FP8` for the extraction (an OpenAI-compatible endpoint given
by `$LLM_ENDPOINT` or `--endpoint`).

Python 3.10 or later; `pip install -r requirements.txt`. Statistics were run
with numpy 1.26 / scipy 1.11 (unimatrix) and numpy 2.2 / scipy 1.15 (local);
results agree to floating-point precision. Figures need graphviz (`twopi`).
