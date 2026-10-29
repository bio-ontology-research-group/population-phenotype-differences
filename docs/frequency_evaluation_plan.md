# Frequency annotation and ranking evaluation

The four disorder targets and the 18 available clinical queries were frozen in
`results/v5/frequency_sources/design.json` before any frequency-ranking result.
The evaluation is exploratory because the clinical sample is small. SPINK5 is a
provenance case study unless an additional independent evaluation cohort is secured.

## Annotation eligibility

Maintain one record per source cohort, disorder, phenotype, assessment method and
age/time window. Record recruitment, genetic confirmation, observed population,
place of care, numerator, assessed denominator, unknown count, source locator and
patient overlap. Treat the source paper and its supplements as one source. Preserve
conflicting statements for resolution; an automated quotation match establishes
text provenance only, not correctness of the denominator or cohort membership.

Binary patient tables can supply phenotype-specific denominators when their symbol
definitions distinguish presence, absence and unknown/unassessed. Positive-only
clinical summaries supply positive observations but generally leave assessment
coverage unresolved. Single-case 1/1 and explicit 0/1 observations are retained
separately: pooling only features selected for mention would create frequencies
near one. The primary frequency update therefore requires an explicit cohort
assessment scheme or a complete, interpretable phenotype-by-patient table.
Case-report mention frequencies may be assessed separately and labelled as such.

Check genetic subtype directly, especially KMT2D versus KDM6A and clinical Kabuki
diagnoses without a molecular assignment. Do not assign a multi-reference HPOA
frequency to one paper or infer a population from an author's country. Keep
recruitment geography, stated ancestry/origin and care setting separate. A single
mixed-ancestry cohort has no population-specific denominator until its patient-level
membership and assessments can be linked.

## Within-disorder comparison

Compare the same phenotype and assessment definition within the same genetic
disorder, with age and ascertainment described. Report exact binomial uncertainty
for genuinely assessed denominators and study-level heterogeneity. Relatives and
repeated publications require family/cohort clustering or deduplication. Pool only
nonoverlapping, compatible cohorts. If provenance is unknown or study and population
are inseparable, report the missingness or nonidentifiability instead of fitting a
population coefficient. These are selected clinical cohorts; frequencies are
conditional on recruitment and assessment.

## Independent ranking

Verify each clinical query's causal gene, eligible disorder subtype and independence
from source patients. Remove evaluation publications and reporting-group overlap
from both baseline and updated reference construction. A pooled HPOA annotation
containing an excluded reference must be excluded when its contributions cannot be
separated. Audit residual overlap by family, variant and case description where
available. Freeze query phenotypes and use the same candidate set and scoring
method in all comparisons.

The benchmark must separate the effect of adding phenotype associations from that
of adding frequency information. Compare an unweighted reference, the same support
with existing HPOA frequencies, and the same support with verified updated
frequencies. Evaluate any new association support in a distinct comparison. Freeze
the frequency-aware scoring rule and handling of counts, categories, missingness
and uncertainty before examining rank outcomes. Categories retain interval bounds;
missing frequency is never set to zero. Any representative point value used in
scoring requires interval-endpoint sensitivity analyses.

Report every eligible query's rank, reciprocal rank and top-10 status under each
reference, with paired changes and cohort/gene breakdown. Resample independent
families or publications rather than treating repeated patients as independent.
With 18 candidate queries, estimates will be imprecise; a negative or inconclusive
result is retained. Improve annotation provenance regardless of whether ranking
improves. Avoid tuning the update or scoring rule on these test outcomes.

Current status: source retrieval complete for accessible articles; source curation,
patient-overlap verification, frequency model specification and evaluation pending.


## Frozen scoring specification (before first rank evaluation)

Use disease profiles for positive phenotypic HPO annotations attached to OMIM
identifiers with gene mappings in the frozen HPO release. Remove annotations whose
references include a held-out clinical publication or a known matching last-author
reporting group. Exclude NOT annotations and the excluded-frequency category.
Preserve the remaining phenotype support in every frequency-only comparison.
Use the verified KMT2D mapping for Kabuki syndrome 1 (OMIM:147920); the frozen
cross-gene mapping of that identifier to KDM6A is removed in every arm. Other
mapped OMIM disease profiles remain available as candidates.

Reuse Resnik best-match-average similarity. The query-to-disease direction remains
an unweighted mean of each query term's best match. In the reverse direction,
weight each disease term's best query match by its phenotype-frequency weight and
normalize by the sum of weights. Average the two directions. A gene receives the
maximum score of its mapped disease profiles. Rank all genes with eligible
profiles; use midranks for ties and also report the best and worst rank in a tie.
Information content is computed from gene frequencies in the filtered, unweighted
reference and held fixed across all arms. With equal weights this reduces to the
unweighted BMA comparator.

For existing HPOA frequencies, use unqualified age/sex/modifier-independent
annotations only. A single count value uses (affected + 0.5)/(assessed + 1), a
Jeffreys-smoothed proportion. Duplicate identical frequency values count once;
conflicting values are left unresolved rather than pooled without cohort identity.
Category intervals are those defined in the frozen ontology (100%, 80–99%, 30–79%,
5–29%, 1–4%, 0%); their midpoints are scoring representatives, not inferred exact
frequencies. Percentage values retain their stated values. Missing or unresolved
frequency has a neutral weight of one. Sensitivities use both category endpoints
and a missing-frequency weight of 0.5. These settings are all reported, never
selected by the test ranking.

Updated frequencies replace the existing weight for the same disorder and exact
HPO identifier using the source-checked cohort counts, after overlap checks.
Within the update, pool assessed counts only for compatible, nonoverlapping source
cohorts. Do not add the HPOA counts to the new counts because their source patients
may overlap. A sensitivity gives each family equal contribution to affected and
assessed totals, preserving unknown assessments. Another uses only the source's
reported truncating KMT2D subset from PMID 25896430. Report association-only and
association-plus-frequency expansion separately, with unchanged information
content and candidates, so a new term cannot masquerade as a frequency benefit.

The eligibility audit retains 12 of the 18 frozen candidate queries. Five have
unsolved/undetermined source status; one is a solved MPV17 neuropathy case outside
the hepatocerebral annotation target. The five US queries have consistent KMT2D
gene targets but discrepancies in within-gene patient/variant linkage across the
original supplementary tables and prose. Report a prespecified sensitivity
restricted to the one US query with concordant ID, age and gene, leaving eight
queries overall. Do not repair patient identity by guessing from phenotype
similarity. Source-derived query phenotypes remain frozen.
