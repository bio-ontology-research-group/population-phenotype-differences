# Clinician reference and extraction evaluation

Paul and Marwa will independently annotate 40 source documents from 20 publications.
The sample was drawn before inspecting the new model predictions: four publications
per region when available, with two each for South Asia and sub-Saharan Africa.
Both abstract and full-text views are included. Model outputs, pairing keys and
assigned analysis regions are withheld. Source type and geography can remain evident
in the text. The sample supports an audit of this selected accessible set; it does
not estimate performance across all publications without additional assumptions.

Each reviewer first marks document review completion and records eligible primary
patients, findings and explicit assessed denominators. Empty files are pending work.
A reviewed document with no eligible patients requires an explicit completed status
and a note; this is distinct from a missing return. The reference must distinguish
individuals, families and aggregates. Previously published comparison cases and
unaffected relatives are ineligible. Normal tests support only the specific absence
the source and reviewers establish. Unmentioned features remain unknown.

Before adjudication, report agreement on patient eligibility, exact HPO sets and
present/absent status, plus the number of disagreements in observation basis,
ancestry, care and frequency denominators. Set precision/recall/F1 and Jaccard are
appropriate for independently generated finding lists; an artificial universe of
all unmentioned HPO terms must not supply true negatives. Identity alignment is
reviewed before counting agreement. Preserve both original returns, their hashes
and the adjudication record. Resolve disagreements jointly, recording rationale
and unresolved items; model output is revealed only after the reference is frozen.

The primary extraction metrics are precision, recall and F1 for distinct
patient–finding–status records against the adjudicated reference, separately for
abstracts and full texts. Count incorrect patient attribution as an error. Report
patient/aggregate eligibility separately so that literature-table contamination
cannot disappear through phenotype deduplication. Report both pooled counts and
document macro averages, with explicit conventions for reviewed empty documents.
An empty prediction for a nonempty reference has recall zero; when no positive
prediction exists, precision is undefined and its denominator is reported.
Incomplete and unresolved documents are excluded from the primary score with their
counts and reasons shown; their effect belongs in a sensitivity analysis.

Separate extraction from mapping. After reference freezing, reviewers match
predicted source phrases to reference findings without requiring an exact HPO
identifier. This establishes source-supported extraction precision and recall.
For those matched phrases, score exact HPO mapping and record excessive specificity,
overly broad mappings, incorrect anatomy, incorrect finding type and null mappings.
Ancestor-compatible matching is a secondary result, never relabelled as exact
accuracy. Do not award multiple true positives for repeated descriptions of the
same finding. Evaluate the bedside/test classification against the clinician's
observation basis for matched findings, with mixed and unclear findings reported
separately.

For ancestry and place of care, report supported field accuracy and missingness
separately, including unsupported inference from affiliation or care location.
For frequency extraction, score exact affected/assessed count agreement only on
matched cohort–phenotype statements. Wrong cohort membership, pooling of historical
patients, mistaken assessment denominators and source contradictions are separate
error categories. The absence of explicit frequency statements is a valid reference
outcome, not permission to fabricate a total-cohort denominator.

Resample publications, retaining both views and all patients from each publication,
for 95% bootstrap intervals and paired abstract/full-text accuracy differences.
Report regional counts and error rates descriptively; two publications in a region
cannot support precise regional performance claims. The single-individual
machine-output subset is only a sensitivity analysis until identity is checked.
The existing model-judge audit and these clinician measurements must remain
separately labelled in Methods, Results and the supplement.

Status: packets prepared; independent returns, adjudication and scores pending.
