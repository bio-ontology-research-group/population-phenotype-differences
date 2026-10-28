# Extraction audits

Two region-blind audits of the language-model extraction (Qwen3.8-27B, FP8) and
HPO mapping were run, one on PMC-Patients summaries and one on PubMed
abstracts. In both, a second language model (Claude Sonnet, run as an agent
with the case text and the extracted phrases) acted as judge. The judgements
are model-versus-model agreement; a clinician validation is pending.

## Sampling

- PMC-Patients: 60 cases, 20 each from EASIA, EUR and MENA (the three
  populations of the earlier analysis), shuffled into three batches of 20.
  Judgements: `results/v4/audit_pmcp/batch{1,2,3}_judged.json`. The batch
  inputs held the PMC-Patients summary text (CC BY-NC-SA) and are not
  redistributed.
- PubMed abstracts: 120 cases, 20 per population (EUR, MENA, SASIA, LATAM,
  AFR, EASIA), one patient per publication, drawn with random seed 2026 and
  shuffled into four batches of 30. The judge saw the title and abstract.
  Judgements: `results/v4/audit_abs/batch{1,2,3,4}_judged.json`; the
  population of each case: `results/v4/audit_abs/regions.json`.

The judge never saw the population label, the author list or the affiliation.

## Judgement scheme

For every extracted phrase the judge recorded

- extraction: `correct`, `not in text`, `negated` (the text excludes the
  finding), `other person` (a relative, not the proband), `not a phenotype`;
- mapping of the phrase to its HPO class: `exact`, `broader` (a correct
  ancestor), `too specific`, `wrong`;
- how the finding is established in the text: test-based (laboratory, imaging,
  physiological test, pathology) or bedside (examination sign, symptom or
  history); this is compared with the per-class finding class in
  `results/termclass.jsonl`;

and, per case, the findings present in the text that the extraction missed,
each with the same test-based or bedside label.

## Scoring

Precision = correct phrases / all extracted phrases. Recall = correct phrases /
(correct phrases + missed findings). Recall by finding class uses the judge's
per-mention label. Mapping shares are over correct phrases. The "test share
extracted" is the test-based share among the correctly extracted phrases and
the "test share in text" is the same share after adding the missed findings.
The summary tables are in `RESULTS.md` ("Direct audit" for PMC-Patients;
"Extraction audit on PubMed abstracts" for the abstracts) and in
Supplementary Table S2 of the paper. The agreement between the judge's
per-mention label and the per-class finding class is 94.2% of 450 classified
phrases (abstract audit).

## Benchmark against curated full text

Independently of the audits, the abstract pipeline was run on the abstracts of
the 1,488 Phenopacket Store publications and compared with the curators'
full-text HPO classes (`validate_against_store.py`; outputs in
`results/v4/validation/abs_vs_store_*.csv`), and the PMC-Patients pipeline with
the same curation where publications overlap (`pmcp_vs_store_*.csv`).
