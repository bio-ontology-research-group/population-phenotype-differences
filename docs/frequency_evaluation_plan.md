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
