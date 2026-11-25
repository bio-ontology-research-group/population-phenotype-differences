# Independent review of phenotype extraction from 20 rare disease case reports

Version 1.0.0. Licence: CC BY 4.0 for the annotations (see "Licence" for quoted text).

## What this is

Two reviewers (Paul N. Schofield, code `PNS`; Marwa Abdelhakim, code `MA`) independently read the
abstract and the full text of 20 published case reports of Mendelian disorders and recorded, for the
patients with the named causative gene:

- every phenotypic finding, with the supporting phrase, whether it was present or explicitly absent,
  its Human Phenotype Ontology (HPO) class, and how the finding was established in that report
  (history, examination, laboratory, imaging, physiological test, pathology, mixed, unclear);
- each patient or group of patients, with diagnosis, causal gene, variant, allelic state, stated
  ancestry and country of care;
- explicit frequency statements (reviewer `MA` only);
- whether the document was complete.

The dataset also contains the output of the automated extraction and HPO mapping procedure of the
accompanying study for the same 40 documents, so that the published agreement figures can be
recomputed. It can serve as a small reference set for phenotype extraction, HPO normalisation and
inter-annotator agreement.

## How it was made

- Sample: 20 publications drawn at random within six regional groups (EUR 4, MENA 4, East Asia 4,
  Latin America 4, South Asia 2, sub-Saharan Africa 2; `prepare_human_review.py`, seed 20261031)
  from publications with both a PubMed abstract and an open-access PubMed Central full text.
- Presentation: abstract and full text as 40 separate documents in random order with neutral
  identifiers (D001 to D040). Model output and regional group were withheld. Geography that is
  stated in the text remained visible.
- Instructions given to the reviewers: `protocol.txt` (verbatim).
- The reviewers worked independently. Disagreements were not adjudicated; this is two parallel
  references, not a consensus.
- Reviewer `MA` returned no findings for document D004. Reviewer `PNS` returned no frequency or
  document-assessment sheet.
- HPO release: 2026-01-08 (`http://purl.obolibrary.org/obo/hp/releases/2026-01-08/hp.obo`).

## Files

All tables are UTF-8, tab-separated, one header row, empty string for missing values.
`datapackage.json` (Frictionless Data Package) gives every column with type, description, allowed
values, keys and a SHA-256 checksum per file.

| File | Rows | Content |
|---|---|---|
| `documents.tsv` | 40 | Document key: PMID, PMCID, DOI, title, licence, target gene, regional group |
| `findings.tsv` | 1,788 | Findings recorded by the reviewers |
| `patients.tsv` | 90 | Patients and aggregates recorded by the reviewers |
| `frequencies.tsv` | 91 | Frequency statements (reviewer MA) |
| `document_assessments.tsv` | 40 | Completeness of each document (reviewer MA) |
| `pipeline_findings.tsv` | 614 | Automated extraction and HPO mapping for the same documents |
| `finding_classes.tsv` | 758 | Finding class assigned to each HPO class in the study |
| `protocol.txt` | | Reviewer instructions |
| `datapackage.json` | | Machine-readable description of all tables |

`document_id` links every table to `documents.tsv`. Publications are identified by PMID
(`https://pubmed.ncbi.nlm.nih.gov/<pmid>`), PMCID and DOI; HPO classes by CURIE (`HP:0001250`,
resolving to `http://purl.obolibrary.org/obo/HP_0001250`); genes by HGNC symbol.

## Raw and normalised values

Reviewer entries are kept as entered in the columns ending in `_as_entered`. The columns `hpo_id`,
`hpo_label`, `status` and `observation_basis` are normalised:

- cells listing several HPO identifiers give one row per identifier (same `row` number);
- identifiers are resolved in the HPO release above (alternative identifiers replaced by the primary
  one; `hpo_id` is empty if the identifier is not in that release, which affects 45 rows);
- a basis cell with several values becomes that value if the others are "unclear", otherwise `mixed`.

Free-text columns of `patients.tsv` and `frequencies.tsv` are as entered, with line breaks removed.

## Reproducing the agreement figures

The scripts are in the code repository
(https://github.com/bio-ontology-research-group/population-phenotype-differences):
`clinician_review_normalize.py`, `clinician_validation.py`, `clinician_validation_patients.py`,
`clinician_agreement.py`, and `build_review_dataset.py`, which wrote this dataset. Summary outputs
are in `results/v5/clinician_validation/` of the repository.

## Licence

The annotations and tables are released under CC BY 4.0. The source documents are not included.
Phrases quoted from them (`source_phrase`, `source_locator`, evidence columns) remain under the
licence of the publication given in `documents.tsv` (CC BY, CC BY-NC, CC BY-NC-SA or CC BY-NC-ND) and
are attributed by PMID and DOI.

## Citation

Please cite the accompanying article and this dataset (DOI assigned by Zenodo on deposit).
