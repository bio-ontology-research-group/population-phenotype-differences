# Reproducing the source, genetics and frequency extension

Run from the analysis repository. Frozen ontology inputs are the original HPO
2026-01-08 release in `data/hpo/`. Existing v4 ranks and case tables are not rewritten.
All public results are under `results/v5/`; source XML, quotations, model responses
and blinded review packets are intentionally under ignored `data/v5_private/`.
No credentials are stored in these results or commands.

## Source and genetic analyses

```sh
python3 source_robustness.py --out results/v5/source_robustness
python3 hpoa_provenance_audit.py --out results/v5/hpoa_audit
python3 genetic_context.py --metadata results/v5/hpoa_audit/case_metadata.csv --out results/v5/genetic_context
python3 gene_regional_analysis.py --metadata results/v5/hpoa_audit/case_metadata.csv --out results/v5/gene_regional
python3 gene_regional_analysis.py --metadata results/v5/genetic_context/case_metadata.csv --versions disease_inheritance,protein_consequence_allelic --out results/v5/gene_genetic_context
python3 gene_influence.py --metadata results/v5/hpoa_audit/case_metadata.csv --out results/v5/gene_influence
python3 paired_source_analysis.py --termclass-extra data/v5_private/paired_mapping_final/termclass_extra.jsonl --out results/v5/paired_source_analysis
python3 frequency_overlap.py --hpoa data/hpo/phenotype.hpoa --out results/v5/frequency_overlap
```

The paired comparison requires the archived successful extraction and mapping
outputs, not a new stochastic model run. Their hashes, prompts, model, generation
settings, mapping cache reuse and missing raw-selection-response limitation are in
`results/v5/paired_source_analysis/pipeline_provenance.json`. The source selection
manifest preserves all 110 sampled publications, including inaccessible full texts.
The proposal extraction pass on 35 frequency sources is an aid to source review;
its 294 proposals are not a validated frequency dataset.

## Frequency curation and ranking

```sh
python3 curate_frequency_tables.py --raw data/v5_private/frequency_sources --obo data/hpo/hp.obo --out results/v5/frequency_curation
python3 curate_frequency_cohorts.py --tables results/v5/frequency_curation --raw data/v5_private/frequency_sources --obo data/hpo/hp.obo --out results/v5/frequency_curation_combined
python3 clinical_evaluation_audit.py --pavs /home/leechuck/Public/software/pavs --out results/v5/clinical_evaluation_audit
python3 frequency_benchmark.py --out results/v5/frequency_benchmark
python3 frequency_source_contrast.py --out results/v5/frequency_source_contrast
python3 extension_tables.py --out ../paper/extension_tables.tex
```

Replace the PAVS path with the local checkout containing the original
`combined_normalized-gpt4oss.tsv`; its hash is saved in the eligibility audit.
Original publisher spreadsheets for PMID 31618753 established the US within-gene
identity discrepancy. Do not assume its phenotype/variant records can be repaired
by rearranging similar cases. The 12-query primary and eight-query stringent
identity sets preserve the frozen phenotypes and gene targets.

The scoring design was specified before examining ranks. Main frequency results
use four source cohorts; patient/family units, unknown assessments and source
discrepancies remain explicit. The later German RPGRIP1 provenance comparison is
separate and was not added retrospectively to that ranking annotation set. Results
do not establish diagnostic benefit or a population-specific frequency effect.

The per-query ranking file includes every sensitivity arm, best/worst tied ranks,
the highest-scoring disease profile for the target gene and the intended-disorder
score. Publication bootstrap intervals use only three query publications and are
descriptive. Source overlap checks use publications, known reporting groups,
stated cohorts/care locations and reported coding variants; they do not replace
private-identifier linkage or normalized cross-transcript equivalence.

## Verification and manuscript build

```sh
python3 -m unittest discover -s tests -p test_extension_statistics.py
python3 -m unittest discover -s tests -p test_frequency_benchmark.py
```

Checks cover clustered resampling, source confounding, lack of care overlap,
unknown versus absent assessment, variant-syntax scope, paired missingness,
multi-reference leakage, frequency conflicts, family weighting, tie midranks and
a hand-calculated frequency-weighted Resnik example.

The manuscript imports `extension_methods.tex`, `extension_results.tex` and the
generated `extension_tables.tex`; the two extension figures are copied from
`results/v5/figures/`. From the manuscript repository:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../review-work/extension-build main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../review-work/extension-build gim_supplement.tex
```

Compiled outputs and visual-review renders stay in the project-local review
directory. The manuscript remains a draft while independent clinician reviews,
adjudication and scores are pending. See `human_validation_analysis_plan.md` for
the reference-freezing and scoring procedure. Empty return files must never be
interpreted as completed negative assessments.
