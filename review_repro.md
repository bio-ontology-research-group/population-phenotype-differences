# Reproducibility review of the analysis repository

Reviewed 2026-10-28 (snapshot around 21:45 local time). Nothing in the repository was modified except this file. All
dry runs wrote only under `/home/leechuck/tmp/repro_check/`. Note: README.md, requirements.txt, gim_figures.py,
population_model.py, corpus_replication.py and jobs/r2_um.sbatch were edited while the review ran (mtimes 15:28 to
21:41). The findings below refer to the files as they stood at the end of the review.

## Summary

The core statistics are in good shape. Every permutation and bootstrap script has a fixed seed, the permutation
counts that can be checked from p-value floors match the paper (2,000 for terms, outcomes, recovery, GDA, symptoms,
age/sex, setting and co-reporting; 1,000 for similarity; 500 country bootstraps; 1,000 case bootstraps for the GDA
loss), and all scripts parse under Python 3.10 and 3.11. Two dry runs reproduced exactly: `qvalues.py` regenerated
`results/v4/qvalues.tsv` (2,482 rows, 0 differences in q), and `age_sex_interaction.py` gave the same output under
two `PYTHONHASHSEED` values. `gim_figures.py` rebuilt 7 of 8 figures pixel-identical. The repository contains no keys
or tokens.

The main gaps are:

1. **Undocumented commands.** For most analyses added on 2026-10-28, and for the aggregation steps that build the
   tables the paper cites, the exact command is recorded nowhere.
2. **Audit without code or inputs.** The extraction audits reported in the paper (Table S2a/b) have no code, and the
   PMC-Patients audit inputs are not in the repository.
3. **Quoted source text in tracked files.** Tracked JSONL files contain verbatim quotes from the source articles,
   including Wiley TDM full text. This contradicts the README, and the licence metadata does not cover the
   PMC-Patients share-alike terms.
4. **Inputs a reader cannot obtain.** Several inputs come with no instructions: the OntoFUNC termdb build, publication
   years for abstracts, the PAVS phenopackets, HGNC, the Phenopacket Store release and the PMC-Patients download.
5. **Stale or mismatched outputs.** RESULTS.md is out of date in one place. The current `gim_figures.py` would label
   the Figure 3 panels wrongly.

## Findings by severity

### Critical

**C1. The command that produced most v4 results is not recorded** (README.md, RESULTS.md, jobs/)
- Problem: README says "The exact commands for each result are in RESULTS.md and jobs/". That holds only for the
  terms (jobs/terms_um.sbatch + terms_tasks.txt), setting (setting_um.sbatch + setting_tasks.txt), allpairs,
  vs-rest, recovery and GDA rank jobs. The following have only a generic usage template in a docstring, or nothing
  at all:
  - `fine_regions.py` (results/v4/fine), `add_age_sex.py` (results/v4/fine2) and `setting_labels.py`
    (results/v4/setting2).
  - `symptom_sign.py` (symptom_sign/), `age_sex_interaction.py` (age_sex/), `setting_pooled.py` (setting_pooled/),
    `comorbidity_net.py` (comorbidity/organ*, including the `--split` run and the undocumented pilot_organ_*), and
    `healthspend.py` (worldbank/healthspend_*).
  - `validate_against_store.py` (validation/abs_vs_store_*, pmcp_vs_store_*), `corpus_replication.py`
    (validation/replication_*), `abs_to_cases.py` on the Store abstracts (validation/cases_storeabs.tsv), and
    `validate_extraction.py` (validation/pmcp_map2_validation.csv).
  - `gda_impact.py stats` (gda/gda_stats.csv, stats.log).
  - `network_figs.py` shift/dag/spend/popnet, including the `--print-layout` runs that made gim_fig4_popnet and
    gim_fig3_shift_nodes.csv. `gim_figures.py` itself is also unrecorded.
- Problem: several files the paper relies on are built by ad hoc code that is not in the repository:
  - `allpairs/contrasts.tsv` and `within/contrasts.tsv`: flattened from popmodel_*_allpairs.json, with q_BH added.
  - `terms_all/all_outcomes.tsv`: the 50 outcome CSVs concatenated, with an `excess` column.
  - `setting_all/test_share_q.tsv`: BH over 15 contrasts.
  - `gda/gda_stats_q.csv`: BH within version added to gda_stats.csv.

  `gim_figures.py`, `network_figs.py popnet` and `qvalues.py` read these files, so Figures 1, 2 and 4 and the
  similarity and setting q-values cannot be rebuilt from code.
- Problem: the copy step from the cluster layout (terms/, gda/, results/fine/...) into results/v4/* is not
  recorded. The same goes for job 10488_1 (setting task 1; setting_um.sbatch only covers array 2-15) and for the
  PMC-Patients audit inputs.
- Fix: add `jobs/run_v4_local.sh` (or a Makefile) with the exact command line of every local step, in order, using
  repository-relative paths. Add a small `collect_v4.py` that builds contrasts.tsv, all_outcomes.tsv,
  test_share_q.tsv and gda_stats_q.csv from the raw outputs, and call it from that script. Each RESULTS.md section
  from "All populations against all" onward should name the script line or job that produced it. Have every script
  print `sys.argv` and the git commit at the top of its log.

**C2. The extraction audits in the paper have no code, and part of their input is missing** (results/v4/audit_abs/,
RESULTS.md lines 285-305 and 605-627; paper main.tex line 96, supplement Table S2)
- Problem: the abstract audit numbers (precision 93.7%, recall 92.2%, mapping 80.3/14.5/5.1%, 94.2% finding-class
  agreement on 450 phrases, per-population table) come from `batch*_judged.json` plus `regions.json`. No script
  draws the 120-case sample (seed 2026), builds the region-blind batches, holds the judge prompt or scores the
  judgements. The PMC-Patients audit (60 cases, 887 phrases, precision 98.4%, recall 90.4%; reported in main.tex)
  points to `scratchpad/audit/batch*_judged.json`, which is not in the repository.
- Fix: add `audit_sample.py` (sampling with the seed, batch shuffling), the judge prompt as a text file, and
  `audit_score.py` that turns the judged JSON into the Table S2 numbers. Copy the PMC-Patients audit batches into
  `results/v4/audit_pmcp/` (or `results/audit_pmcp/`). State in README that the judges were Claude Sonnet agents, so
  re-judging is model-dependent but scoring is exact.

**C3. Tracked files contain quoted source text, and the licence metadata does not match the data** (.gitignore,
README "Not in git", .zenodo.json)
- Problem: README says "Raw LLM responses that quote source text" are not in git. These tracked files carry
  verbatim `evidence` quotes from the articles:
  - `results/origin.jsonl` (2,541 quoted evidence strings) and `results/origin_merged.jsonl` (3,142);
    `results/origin_merged_v3.jsonl` is likely similar.
  - `results/pmcp/origin.jsonl` (12,162) and `results/pmcp/setting.jsonl` (9,321).
  - `results/origin_tdm_wiley.jsonl` (174 records, evidence quoted from Wiley full texts obtained under the KAUST
    TDM subscription).

  `results/pmcp/fetch.jsonl` and `results/abs/abs_*.jsonl` are tracked too; they hold metadata and extracted
  phrases only. Redistributing quotes from subscription TDM texts may breach the publisher's TDM terms.
- Problem: `.zenodo.json` declares `"license": "cc-by-4.0"` for the whole deposit, which includes BSD code and the
  PMC-Patients-derived records. Those records are CC BY-NC-SA 4.0, and share-alike means they cannot be relicensed
  under CC BY. README and the paper describe this correctly; the deposit metadata does not.
- Problem (to check): `results/phenopackets/pavs_clin/` (2,576 PAVS Saudi clinical cases) is redistributed under
  CC BY 4.0. Confirm that the PAVS licence and the data-use terms allow this.
- Fix:
  - Strip the `evidence`, `notes` and `consanguinity_evidence` fields from the tracked origin/setting JSONL, or move
    those files to .gitignore and ship a derived TSV (pmid, labels) instead. Drop `origin_tdm_wiley.jsonl` from the
    repository in any case.
  - Because the repository was pushed, purge these files from history before making it public.
  - In .zenodo.json, either split the PMC-Patients-derived files into a separate deposit under cc-by-nc-sa-4.0, or
    state both licences in the description and set the record licence to the most restrictive one that applies.

### High

**H1. Inputs that are not in the repository and have no instructions** (README "Not in git", data/get_hpo.sh)

| Input | Used by | In repo? | README says how to obtain? |
|---|---|---|---|
| HPO 2026-01-08 (hp.obo, hp.owl, phenotype.hpoa, genes_to_phenotype.txt, phenotype_to_genes.txt) | most scripts via `--hpo-dir`, `--obo` | no | yes: `data/get_hpo.sh`. All five URLs resolve (HTTP 200). |
| OntoFUNC termdb (term.txt, term2term.txt, graph_path.txt) | `population_model.py terms --termdb`, comorbidity_net, phenotype_specificity, hpo_map2 | no | partly: "built with ELK", but no tool, version or command |
| Similarity cache simcache_v4.npy | population_model, recovery_main, retrieval_perm | no | yes: recomputed by case_pairs.py. Needs about 17 GB and RAM to match. |
| World Bank JSON | healthspend.py | yes (results/v4/worldbank) | no: download URL and date not recorded |
| Publication year per abstract PMID (`abstracts.jsonl` from pubmed_abstracts.py) | healthspend.py `--years-abs` | no (contains abstract text) | no. Without it the healthspend analysis cannot be rerun. |
| HGNC `hgnc_complete_set.txt` | abs_to_cases, pmcp_to_cases, merged_model, to_phenopackets | no | no |
| Phenopacket Store 0.1.26 JSON (`<sp>/ppktstore/0.1.26/*/*.json`) | to_phenopackets, case_pairs `--store` | no | no |
| PAVS_phenopackets.json | to_phenopackets (hardcoded path), merged_model, run_analyses | no | no |
| PMC-Patients V2 | pmcp_extract | no | no |
| BGE embedding model | hpo_map2.py | no | no (model name only in code) |
| Natural Earth 110m shapes | gim_figures.py `--shapes` | no | yes (added during review) |
| LLM endpoint (Qwen3.8-27B-FP8 via vLLM) | qwen_tasks, pmcp_extract, store_origin_extract | n/a | partly. The default `--endpoint` is the internal host http://unimatrix01.kaust.edu.sa:8000; the jobs/orix scripts show a self-hosted vLLM recipe. |
| API keys (BORG LLM key, Elsevier/Wiley TDM) | the scripts above, tdm_fetch | correctly absent | `--key-file` and the env file are documented in docstrings; README mentions TDM keys only |

- Problem: README says the HPO goes to `data/hpo/` and the termdb to `data/ontofunc_hp/`, but every job uses
  `--hpo-dir hpo`, `--obo hpo/hp.obo` and `--termdb ontofunc_hp` relative to the working directory.
- Fix:
  - Add `data/build_ontofunc.sh`, or a README paragraph giving the ELK/OntoFUNC version and command.
  - Ship `results/abs/pmid_year.tsv` (pmid, year, journal; no abstract text) and point healthspend at it.
  - Add a "How to obtain" line per input, with URL and version: HGNC download date, the Phenopacket Store release
    zip, PAVS, PMC-Patients V2 (Hugging Face/Figshare) and the BGE model id.
  - Record the World Bank API URL and the download date.
  - Use the same layout (`data/hpo`, `data/ontofunc_hp`) in the jobs, or say in README to symlink them.

**H2. RESULTS.md and code are out of step with the result files** (RESULTS.md lines 648-652; population_model.py,
corpus_replication.py, gim_figures.py)
- Problem: the abstract-bias numbers in RESULTS.md (MENA minus EUR 95% CI -10.7 to +1.8, p=0.16; LATAM p=0.11) do
  not match `validation/replication_abstract_bias.csv` (CI -11.13 to +1.98, p=0.152; LATAM p=0.126), which matches
  the paper (main.tex line 190). corpus_replication.py was modified at 15:28 and rerun at 15:29, after RESULTS.md
  was last written at 15:20.
- Problem: population_model.py (modified 20:17) is newer than every output it produced (allpairs 10-27, terms
  10-28 morning, vs-rest 10-28). gim_figures.py (21:41) is also newer than the committed figures. The repository
  does not record which code version produced which output.
- Fix: update RESULTS.md from the CSV. Before release, rerun the cheap local steps (qvalues, gim_figures,
  corpus_replication, setting_pooled, symptom_sign, age_sex_interaction, comorbidity) with the final code and check
  that nothing changes. For the cluster outputs, confirm that the population_model.py diff since job 10204 touches
  only new options (for example `--r2-only`). Write the git commit hash into every log.

**H3. Regenerating Figure 3 now gives the wrong panel letters** (gim_figures.py lines 228 and 263)
- Problem: a dry run of `gim_figures.py` (inputs symlinked, output under repro_check) gave pixel-identical
  gim_fig1, fig2_setting, figS_decomp, figS_organ, figS_recovery, figS_agesex and figS_coreport. gim_fig3_shift
  differs only in its panel letters: the code now labels them "b" and "c", while the committed figure (identical to
  paper/figures/gim/gim_fig3_shift.pdf) and the main.tex caption use (a) and (b). The relabelling for the new map in
  Figure 2 seems to have reached fig3 as well.
- Fix: in `fig3`, change `panel_label(a, "b", ...)` to "a" and `panel_label(b, "c", ...)` to "b".

**H4. Job scripts do not match the repository layout, and two of them point at a Claude scratchpad**
(jobs/rerun_v3.sh, jobs/rerun_v4_local.sh, jobs/*_um.sbatch)
- Problem: rerun_v3.sh and rerun_v4_local.sh hardcode
  `/home/leechuck/tmp/claude-1000/.../818691a1-.../scratchpad/{statsenv,hpo,ontofunc_hp,simcache_merged_v3.npy}`
  and `cd /home/leechuck/Documents/papers/...`. The `_um` jobs read `results/fine`, `results/fine2` and
  `results/setting2`, and write `terms/`, `gda/` and `results/fine/`. The repository instead has
  `results/v4/fine`, `results/v4/fine2`, `results/v4/setting2`, `results/v4/terms_all` and so on.
- Fix: make both shell scripts take `ROOT`, `PY`, `HPO` and `TERMDB` from environment variables with
  repository-relative defaults. In the sbatch files, replace `W=/data/hohndor/...` with `W=${W:-$(pwd)}` and use
  the repository paths, or add one README line that maps the cluster layout to the repository layout.

### Medium

**M1. Dependencies** (requirements.txt)
- Imports missing from requirements: `pandas` (pooled_origin_model.py) and `scikit-learn` (followup_analyses.py,
  imported lazily). `statsmodels` is only used by the legacy pooled_origin_model.py. External tools that are not
  documented: `pdftotext` (poppler-utils; store_origin_extract.py, tdm_fetch.py) and vLLM 0.28.0 (pinned only in a
  comment in jobs/orix).
- No lock file. README states numpy 1.26/scipy 1.11 on unimatrix and numpy 2.2/scipy 1.15 locally. The local check
  here ran numpy 2.2.4/scipy 1.15.3 on Python 3.13.5.
- No numpy-2-incompatible API is used (scanned for np.float_, np.in1d, np.trapz, row_stack, copy=False and
  similar). All 37 scripts parse under Python 3.10, 3.11 and 3.13.
- Fix: add `pandas` and `scikit-learn` (or mark them as legacy-only extras). Add a pinned `requirements-lock.txt`
  from the environment used for the final figures and statistics. List poppler-utils and graphviz under "System
  tools" in README.

**M2. Absolute or site-specific paths in Python scripts**

| File:line | Path | Suggested fix |
|---|---|---|
| to_phenopackets.py:278 | `~/Public/software/pavs/data/PAVS_phenopackets.json` | `--pavs` argument (required) |
| to_phenopackets.py:83, 86, 187, 253 | implicit layout under `--sp`: `hpo/hp.obo`, `hgnc/hgnc_complete_set.txt`, `abs/abstracts.jsonl`, `ppktstore/0.1.26/` | separate `--obo`, `--hgnc`, `--abstracts` and `--store-dir` arguments, or document the `--sp` layout in the docstring |
| pmcp_extract.py:218, 222; qwen_tasks.py:243; store_origin_extract.py:306 (and docstring line 13) | default `--endpoint http://unimatrix01.kaust.edu.sa:8000` | default to `$LLM_ENDPOINT` or `http://127.0.0.1:8000`, and document an OpenAI-compatible endpoint |
| tdm_fetch.py:52 | default `--keys ~/.config/pavs-tdm/tdm.env` | acceptable (user config outside the repository); mention in README |
| docstrings of pmcp_extract.py and store_origin_extract.py | `--key-file ~/.borg_llm_key` | acceptable as an example |

  No Python script contains `/data/hohndor`, `/ibex`, `/mnt` or scratchpad paths; these appear only in jobs/ (see
  the appendix). `gim_figures.py --res` defaults to a path relative to the script, which is good.

**M3. README does not document recent scripts and result directories** (README.md)
- Scripts not mentioned: gda_impact.py, symptom_sign.py, age_sex_interaction.py, setting_pooled.py,
  comorbidity_net.py, validate_against_store.py, corpus_replication.py and oa_lookup.py. gim_figures.py is
  mentioned only for `--shapes`, not as the figure generator. The `network_figs.py` subcommands (shift, dag, spend,
  popnet) and `--print-layout` are not mentioned; Figure 4 comes from `popnet --print-layout`.
- Result directories not named: allpairs, within, recovery, terms_all, setting_all, setting_pooled, worldbank, gda,
  symptom_sign, age_sex, comorbidity, validation, audit_abs and figures. Superseded files at the top of
  `results/v4/` are not explained: terms_abs_*, terms_age_*, terms_sex_*, terms_setting_*, the top-level
  terms_all_*, cases_abs_age_*/sex_*/setting.tsv, cases_pmcp_map2_setting.tsv and recovery_v3data_MENA_*.
- Fix: extend the Pipeline list with steps 8 to 11: gene-disease annotations (gda_impact), symptoms vs signs, age
  and sex interaction, pooled setting, co-reporting networks, and validation/replication (validate_against_store,
  corpus_replication). Add a table that maps each `results/v4/<dir>` to its script, command and paper
  table/figure. Move the superseded top-level files to `results/v4/legacy/` or note that they are superseded.

**M4. The similarity cache is not checked against the case set** (case_pairs.py lines 81-82)
- Problem: if `--cache` exists, it is memory-mapped without checking its shape or case order. A stale cache from
  another case set is either used silently or fails later with an index error. The same simcache_v4.npy serves
  results/cases_* and results/v4/fine/* (same order, relabelled regions), which works only because the order is
  unchanged.
- Fix: assert `S.shape == (N, N)` and store a sidecar hash of the ordered case ids (`cache + ".ids.sha1"`); refuse
  the cache on mismatch.

**M5. Permutation counts: defaults and docstrings disagree with what was run**
- comorbidity_net.py: the default and docstring say `--perms 1000`; the paper says 2,000. The logs (p=0.0015 =
  3/2001) show 2,000 were used. Change the default to 2000.
- population_model.py: the default `--perms 1000` serves the similarity analyses; the terms jobs pass 2000 (correct).
  Consider a per-subcommand default.
- corpus_replication.py: run with 1,000 group permutations (p=0.1079 = 108/1001). The supplement (line 202) states
  a floor of 1/2001 for "class, outcome, recovery and country-level tests", and the paper does not give the count
  for the full-text replication p-values. Either rerun with `--perms 2000` or state 1,000 in the supplement. Its
  5,000 label permutations and 2,000 bootstraps for the abstract bias are hardcoded.
- healthspend.py: 500 country bootstraps are hardcoded (line 114). This matches the paper; expose it as
  `--boots 500`. gda_impact.py hardcodes 1,000 case bootstraps, which matches main.tex line 167.
- All 16 stochastic analysis scripts set seeds (`default_rng(3..43)`, `random.seed`). gim_figures.py has no random
  elements. LLM calls use temperature 0. GPU inference is still not bit-reproducible; say so in README.

**M6. qvalues.py writes into the repository on any invocation** (qvalues.py)
- Problem: there is no argparse, so `python3 qvalues.py --help` overwrites results/v4/qvalues.tsv. The docstring
  lists the organ families as total/pc/adjage/adjsex, but the output also has organ:age_* and organ:sex_* families
  (the filter on line 31 never fires because `fam` is already "organ:..."). q-values for gda, symptom_sign,
  setting_all, setting_pooled, age_sex and comorbidity are computed inside those scripts, not here.
- Check: run in an isolated copy, the script reproduced the committed qvalues.tsv exactly (2,482 rows, 18 families,
  0 q differences).
- Fix: add `--res` and `--out` arguments with `--help`. Correct the docstring, or apply the intended filter. Add a
  README sentence on where each family's BH correction is done.

**M7. Licence and citation metadata** (CITATION.cff, .zenodo.json, LICENSE)
- Althagafi's affiliation is KAUST in CITATION.cff and .zenodo.json but Taif University in main.tex.
- ORCIDs are missing for Althagafi, Schofield and Hoehndorf. Verify that 0000-0001-6816-2119 belongs to
  Abdelhakim.
- CITATION.cff: `type: software` with `license: BSD-3-Clause` only; no `version`, `date-released`, `doi` or
  `preferred-citation` (the paper). Use `license: [BSD-3-Clause, CC-BY-4.0]` and add a preferred-citation block once
  the paper DOI exists.
- .zenodo.json: `upload_type: dataset` (fine) with `license: cc-by-4.0` only (see C3). `related_identifiers` uses
  `isSupplementTo` for the GitHub URL; that relation should point to the paper DOI, and the GitHub link should use
  `isDerivedFrom` or `isSupplementedBy`. Add a version.
- LICENSE: "Copyright (c) 2026, The Authors (Bio-Ontology Research Group, KAUST)". Two authors are not at KAUST;
  use "Copyright (c) 2026, Azza Althagafi, Marwa Abdelhakim, Paul N. Schofield, Robert Hoehndorf". LICENSE-DATA is
  the standard CC BY 4.0 legal code. The author list is consistent across all three files and the paper.
- The repository is private, while main.tex (Data Availability) says the code "is available". Make it public (after
  the C3 history clean-up) and fill in the Zenodo DOI before submission.

### Low

**L1. Smoke test (`python3 <script> --help`, run outside the repository with bytecode writing off)**
- 33 of 37 scripts exit 0.
- Failures:
  - `tag_store_pmids_country.py` (no argparse; opens `--help` as a file).
  - `validate_phenopackets.py` (no argparse; `phenopackets` not installed locally, and with it installed `--help`
    would be read as a directory).
  - `pooled_origin_model.py` (`statsmodels` not installed locally; argparse but no main guard).
- Not run with `--help`: `qvalues.py` (would overwrite results; see M6) and `rank_engine.py` (a module; imports
  cleanly).
- Fix: give tag_store_pmids_country.py, validate_phenopackets.py and qvalues.py argparse with `--help` and a
  `__main__` guard; add a main guard to pooled_origin_model.py and validate_extraction.py.

**L2. Orphaned and unrecorded figure outputs** (results/v4/figures/)
- `gim_fig3_shift_nodes.csv` comes from a `network_figs.py shift --print-layout` run whose PDF gim_figures.py later
  overwrote; gim_figures.py reads fig_shift_network_nodes.csv instead. `fig_popnet`, `fig_shift_network`,
  `fig_healthspend` and `fig_dag_AFR`/`fig_dag_LATAM` (empty, no q<0.05 terms) are not in the paper.
- The gim_figures.py docstring lists neither gim_figS_agesex nor gim_figS_coreport. gim_fig4_popnet comes from
  network_figs.py, not from gim_figures.py.
- Fix: record the `network_figs.py popnet --print-layout --out results/v4/figures/gim_fig4_popnet` command (or
  whatever was run), delete or label the unused outputs, and complete the docstring.

**L3. RESULTS.md numbers with no result file**
- The audit tables (see C2) and "Hits@1: Saudi clinical 1.8%, US clinical 1.7%; MENA vs EUR by corpus" (RESULTS
  line 393; derivable from recovery/*_percase.csv, but the computation is not recorded).
- The age and sex composition per population (RESULTS lines 581-583; not in age_sex/log.txt or the CSV).
- "termclass: 2,627 new terms (24 unclassified)" and the corpus counts in the v4 section (derivable from case
  tables).
- Everything else checked has a file: allpairs and within contrasts, the recovery table, the 50 terms runs, the
  healthspend slopes, setting_all test_share_q, the GDA tables, symptom_sign, age_sex, comorbidity global p,
  setting_pooled, abs_vs_store, replication_corpus, and the paper figures (all 13 PDFs used by the paper are
  byte-identical to results/v4/figures).
- Fix: write a small `summaries.py` that prints these descriptive numbers into `results/v4/summaries.txt`.

## Appendix A: traceability of results/v4 directories

| Directory | Script | Recorded command | Gaps |
|---|---|---|---|
| allpairs | population_model.py allpairs | jobs/allpairs_um.sbatch (job 10204; writes results/fine/popmodel_fine) | contrasts.tsv built by unrecorded code; copy step |
| within | population_model.py allpairs --vs-rest | jobs/vsrest_um.sbatch (job 10491) | contrasts.tsv built by unrecorded code |
| recovery | recovery_main.py | jobs/rec_um.sbatch (job 10236) | copy from results/fine; recovery_v3data_* at top level unrecorded |
| terms_all | population_model.py terms | jobs/terms_um.sbatch + terms_tasks.txt (job 10420) | all_outcomes.tsv unrecorded; copy from terms/ |
| setting_all | population_model.py terms | jobs/setting_um.sbatch + setting_tasks.txt (array 2-15; task 1 = job 10488_1 not recorded) | test_share_q.tsv unrecorded |
| setting_pooled | setting_pooled.py | docstring template only | exact command missing |
| worldbank | healthspend.py (+ manual World Bank download) | docstring template only | command, download URL/date and years input missing |
| gda | gda_impact.py rank / stats | rank: jobs/gda_um.sbatch (job 10537); stats: docstring template only (job 10577 not in jobs/) | stats command and the gda_stats_q.csv step missing |
| symptom_sign | symptom_sign.py | docstring template only | exact command missing |
| age_sex | age_sex_interaction.py | docstring template only | exact command missing |
| comorbidity | comorbidity_net.py | docstring template only (says --perms 1000) | commands for organ, organ_split and pilot missing |
| validation | validate_against_store.py, corpus_replication.py, validate_extraction.py, abs_to_cases.py; Qwen on jobs/storeabs_um.sbatch (jobs 10600, 10602) | Qwen stages recorded; the rest only as docstring templates | analysis commands missing |
| audit_abs | none | none | sampling, prompt and scoring code missing (C2) |
| figures | gim_figures.py, network_figs.py | `gim_figures.py [--res]` usage only; network_figs templates | exact network_figs commands missing; fig3 labels (H3) |
| fine, fine2, setting2 | fine_regions.py, add_age_sex.py, setting_labels.py | docstring templates only | exact commands missing |
| qvalues.tsv | qvalues.py | implicit (no arguments) | reproduces exactly |

## Appendix B: absolute paths in job scripts

README says job scripts keep their cluster paths, which is acceptable for the sbatch files. They are:
- `/home/hohndor/pavs-origin`: abs, origin, pmcp, qtasks, recover, tdm_llm and map2/select. Also
  `/home/hohndor/.borg_llm_key`.
- `/data/hohndor/pavs-popmodel`: allpairs_um, gda_um, r2_um, rec_um, setting_um, terms_um, vsrest_um and
  storeabs_um.
- `/ibex/scratch/projects/c2014/pavs-popmodel`: allpairs, popv4 and recv4.
- `/mnt/data/u/hohndor/{qwen-serve,hf,qwen-serve-cu129}`: orix/*.
- `/scratch/$SLURM_JOB_ID`: fallback TMPDIR.

Not acceptable: `jobs/rerun_v3.sh` and `jobs/rerun_v4_local.sh` hardcode a Claude Code session scratchpad
(`/home/leechuck/tmp/claude-1000/.../scratchpad/...`) and the local checkout path (H4).

## Appendix C: secrets scan

The repository has no keys or tokens. The scan covered:
- `sk-`, AKIA, ghp_, xox, Bearer and PEM headers;
- `api_key=` and `token=` literals and 40+ character hex strings in code, jobs and metadata;
- key-like strings in results/, docs/ and data/;
- files named *key*, *token*, *secret* or .env.

Key material is read only from `--key-file` paths outside the repository or from the TDM env file. The ORIX jobs
generate a per-job random key and delete it afterwards.
