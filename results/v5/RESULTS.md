# Extension results: source, genetics and phenotype frequency

New outputs use the original HPO release (2026-01-08). Source, genetic and exploratory frequency analyses are complete at the available evidence resolution and are incorporated in the manuscript draft. Independent review by Paul and Marwa, adjudication and human-validation scores remain pending. Model-output concordance and source checks are not labelled clinician validation.

## Source-specific bedside restriction

The v4 per-case ranks were reanalysed separately by corpus. Eligible cases have both bedside and test-based findings. Rank ratios compare bedside-only with full profiles for the same case. Intervals resample publications (2,000 replicates); reporting-group resampling is supplied as a sensitivity. These are source-specific descriptive effects, not comparisons adjusted to the same gene distribution across sources.

| Corpus | Region | Paired cases | Rank ratio (95% publication-bootstrap interval) | Top 10: full → bedside (%) |
|---|---|---:|---:|---:|
| store | EUR | 2563 | 2.66 (2.15–3.48) | 48.0 → 31.0 |
| store | MENA | 405 | 2.06 (1.54–2.94) | 46.2 → 28.9 |
| store | SASIA | 33 | 5.11 (2.04–11.52) | 66.7 → 33.3 |
| store | LATAM | 60 | 1.83 (1.30–3.57) | 58.3 → 45.0 |
| store | AFR | 22 | 1.05 (0.76–1.49) | 59.1 → 59.1 |
| store | EASIA | 556 | 3.61 (2.55–5.10) | 53.6 → 33.3 |
| pmcp | EUR | 1670 | 3.84 (3.48–4.22) | 40.3 → 19.9 |
| pmcp | MENA | 527 | 3.99 (3.31–4.86) | 46.5 → 24.1 |
| pmcp | SASIA | 191 | 3.69 (2.88–4.78) | 50.8 → 29.3 |
| pmcp | LATAM | 185 | 3.49 (2.54–4.63) | 47.0 → 31.4 |
| pmcp | AFR | 26 | 3.71 (2.23–6.23) | 53.8 → 26.9 |
| pmcp | EASIA | 1856 | 4.02 (3.63–4.42) | 47.6 → 24.9 |
| pubmed_abs | EUR | 8170 | 4.74 (4.53–4.96) | 40.6 → 16.2 |
| pubmed_abs | MENA | 2315 | 4.46 (4.08–4.87) | 41.6 → 17.7 |
| pubmed_abs | SASIA | 732 | 4.84 (4.19–5.67) | 43.4 → 20.4 |
| pubmed_abs | LATAM | 620 | 4.29 (3.62–5.07) | 48.4 → 24.4 |
| pubmed_abs | AFR | 79 | 5.55 (3.69–8.53) | 46.8 → 15.2 |
| pubmed_abs | EASIA | 5224 | 4.70 (4.45–4.97) | 44.5 → 19.0 |

All 18 corpus–region estimates are reported, including the near-null African Store estimate (22 paired cases). The pooled four-fold effect is not uniform across corpora. Source-specific gene distributions, profile lengths and annotation overlap can contribute to this heterogeneity; the corpus contrast is not itself a causal effect of extraction source.

Regional comparisons are in `source_robustness/regional_contrasts.csv`: three phenotype outcomes and three ranking outcomes, each separately within Store, PMC-Patients and abstracts, in total and author-region-adjusted versions. Each outcome family/version controls BH across 45 tests. Null intervals describe permutations, not confidence intervals. New primary phenotype analyses exclude cases with no classified findings; exclusions are recorded in `denominators.csv`.

## HPOA frequency and genetic metadata audit

The frozen annotation file has 280,730 rows: 97,804 with affected/assessed counts, 118,685 frequency categories, 109 percentages and 64,132 missing frequencies.

Count annotations cite 7,711 distinct PMIDs. 20,497 count annotations overlap a publication and explicit disease identifier in the existing case collection. This is a candidate linkage, not verified membership of those patients in the denominator. Ancestry, origin, care and regional assignment remain separate fields. Annotation rows cannot be summed as patient counts.

Detailed genetic metadata coverage and disorder candidates are in `hpoa_audit/`. Store has explicit disease identifiers for all 9,517 entries and variant expressions for 9,088. PMC-Patients and abstract phenopackets carry no normalized disease identifiers or explicit variant expressions; their allelic-state fields do not establish a variant mechanism. Further source extraction/curation is needed for subtype and variant matching.

## Paired-source experiment

A deterministic sample was drawn from the 400 archived single-patient comparison units, with one unit per publication and a target of 30 publications per region before access checks. Of 110 selected publications, 59 currently have usable open-access paired texts, 41 lack a PMC full text, and 10 fetch attempts failed. Accessible pairs: EUR 14, MENA 13, EASIA 22, LATAM 6, SASIA 2, AFR 2. The limited coverage of underrepresented regions must be reported. Fetch failures are distinct from absent full text and remain candidates for retry.

All 118 document views have now been extracted with the same source-neutral prompt; all 1,407 phrase mappings and nine additional observation-type classifications are complete. Existing abstract/full-text mapping results mix extraction and curation processes and are not being relabelled as this new controlled experiment. Human identity alignment, source-text extraction review and denominator verification remain outstanding.

## Outstanding validation

Paul/Marwa independent returns, patient alignment, adjudication and source-text extraction/mapping scores remain outstanding. The manuscript explicitly distinguishes the completed machine analyses from this pending evaluation. The original model-judge audit also lacks archived exact judge settings/prompts, which remain marked as missing in the draft rather than reconstructed by assumption.

## Gene-specific exploratory analysis

With at least 10 cases and three publications per regional group after matching, 342 gene–region contrasts were eligible within corpus, three after conditioning on author region, and six after matching disease identifier and allelic state. No care-country-conditioned contrast with stated origin met these support thresholds. No composition or individual-HPO-term test reached BH q < 0.05 across the complete tested families (4,999 permutations). This reflects limited overlap and multiplicity as well as uncertainty; it does not establish absence of genetic effects. The lowest composition q was 0.0684 for CYP17A1/EASIA test share, an exploratory candidate rather than a confirmed finding. The disease/allelic-state subset includes FBN1, KMT2D, RPGRIP1, SPINK5, TGFBR2 and UMOD; explicit-variant, inheritance and protein-consequence matching are reported below.

The matching audit is in `gene_regional/support.csv`, all tested effects in `contrasts.csv`, and settings and test-family summaries in `summary.json`. No thresholds were relaxed after observing these results. The strong reduction in supported comparisons on conditioning is directly relevant to the scope of any genetics claim.

## Frequency evidence scope

Among count annotations, 88,951 have exactly one reference that is a PMID, 8,262 have multiple or mixed references including a PMID, and 591 have no PMID. A frequency supported by multiple references cannot be allocated to a particular study or population from HPOA alone. The audit now records this explicitly in `frequency_scope`; source cohorts and their denominators must be recovered before a population comparison or pooled frequency update.

Blinded reference packets for Paul and Marwa have been prepared for 20 publications (40 documents each). They contain empty annotation return sheets, source texts, instructions to distinguish unknown from assessed absence, and a separate adjudication template. Their status remains `prepared_not_reviewed`.

## Influence and exact-variant overlap

The deletion analysis recomputes each uncentered, within-gene/corpus Mantel–Haenszel test-share contrast after removing one gene. Largest changes (percentage points): MENA: HBB -0.14; SASIA: ATP7B +0.26; LATAM: CYP17A1 -0.24; AFR: HBB +0.94; EASIA: CYP17A1 -0.11. These are descriptive influence measures, not significance tests.

There are 36 exact variant/allelic-state profiles represented across regional groups, but 0 regional comparisons meet the original support thresholds after matching corpus, exact variant expressions and allelic state. Matching uses genomic expressions where available, otherwise coding expressions; variants have not been normalized across builds/transcripts. No relaxed-threshold confirmatory claim is made.


## Inheritance and protein consequence sensitivity

Seven gene–region contrasts met the support thresholds after matching corpus,
disease identifier and disease-linked inheritance annotations. Four met them after
matching corpus, descriptive protein consequence and allelic state. No composition
or HPO-term test reached BH q < 0.05 in either analysis. Complete outputs are in
`gene_genetic_context/`. Disease inheritance is linked from HPO annotations to an
explicit disorder identifier; it is not a segregation analysis. Protein classes
are recovered from HGVS protein syntax and describe substitutions, frameshifts,
stops and other stated consequences. They do not assign a molecular mechanism.
Store supplied 7,479 records with a classifiable consequence for every available
protein expression. Exact-variant matching remains too sparse for the prespecified
regional comparisons.

## Frequency curation sample and independent evaluation set

Four disorders were selected before frequency-model evaluation: Netherton syndrome
(SPINK5), mitochondrial DNA depletion syndrome 6 (MPV17), Kabuki syndrome 1 (KMT2D)
and Leber congenital amaurosis 6 (RPGRIP1). The frozen list contains 97 candidate
source publications after exclusion of clinical evaluation publications and their
reporting groups. Thirty-five full texts were retrieved, 51 lacked PMC full text,
and 11 fetches failed. Source availability is recorded separately from eligibility.
The initial clinical evaluation list had 18 candidates: three RPGRIP1 and three
MPV17 Saudi cases, plus seven Saudi and five US KMT2D cases. Source diagnosis and
subtype checks retained 12 for the primary benchmark, as detailed below. SPINK5 has
no eligible clinical query and remains a provenance case study.

The frequency extraction development audit found two issues in its initial pilot:
older Kabuki reports using MLL2 were missed, and explicit binary phenotype tables
or single-patient findings were rejected unnecessarily. The pilot is preserved in
private outputs and its revised prompt explicitly handles the historical gene name
and assessment-based denominators. These development papers cannot serve as an
independent validation of that prompt. Source verification of every accepted
frequency remains required; model proposals are not curated frequency annotations.


## Paired extraction output comparison

Full texts yielded 13.92 additional distinct HPO terms per publication on average
(95% paired-publication bootstrap interval 11.61–16.64), including 4.34 additional
test-based terms (3.37–5.44) and 8.88 bedside terms (7.14–10.83). Among the 53 pairs
with a defined test share in both views, the mean change was -0.65 percentage points
(-7.08 to 5.69). For the 48 publications where each output contained one individual,
the term-count difference was 14.48 (11.88–17.75); test-share change was -1.85 points
(-8.30 to 4.54), based on 42 pairs with defined shares. Unknown shares were excluded,
not set to zero. These estimates use 5,000 publication bootstrap draws. No eligible
regional source-interaction test passed BH q < 0.05; groups needed at least five
finite paired observations each, leaving South Asia and Africa descriptive only.

These outcomes compare raw model-generated publication-level phenotype unions.
The one-individual subset is selected from output cardinality; clinician identity
alignment remains pending. At least one full-text output included previously
published comparison patients against the prompt instruction. That output has been
retained for the extraction audit rather than silently corrected. Four abstract
views yielded empty mapped phenotype profiles. Exact abstract-term retention in
full-text output averaged 76.6%; ancestor-compatible retention averaged 85.4%, over
55 nonempty abstract profiles. These are source concordance measures, not accuracy
against a clinical reference. Tables, denominators, mappings/provenance hashes and
all regional comparisons are in `paired_source_analysis/`.


## Study-specific frequency overlap

Across HPOA, 296 disease–phenotype pairs have count annotations attached to multiple
individual PMIDs. Only 79 annotation groups remain comparable after matching onset,
sex and modifiers and restricting to positive phenotypic annotations. Among the four
selected disorders there is one such group, for cryptorchidism in Kabuki syndrome;
all three recorded proportions are 1.0. Source cohort independence and population
provenance are still unverified. HPOA therefore supplies little ready-made
study-specific frequency variation for these selected disorders. Recovering
source-level assessment denominators is necessary before a population explanation
can be tested. Full overlap outputs are in `frequency_overlap/`.

Two KMT2D tables contribute 39 rows and 144 patient–phenotype assessments. PMID
29482518 supplies 17 rows for two boys, with source age/variant discrepancies
retained. PMID 25896430 supplies 22 rows restricted to five author-assigned causal
cases; the two computationally assigned missense cases are excluded in a specified
sensitivity. Blank cells remain unassessed. The primary missense cohort in PMID
31949313 concerns a disorder distinct from Kabuki and was not pooled into the update.

The model proposal pass has completed for all 35 accessible frequency-source papers:
294 candidate statements satisfy integer count bounds, 285 resolve to a proposed
cohort identifier, and 256 contain a whitespace/case-normalized quotation found in
the supplied source. These are technical checks only. Source checking identified
incorrect cohort inclusion and subtype mixing despite valid count syntax. Every
proposal still requires source review before entering an annotation update. The
proposal audit is retained separately in `frequency_proposal_audit/`.

## Frozen assessment-based frequency update

`frequency_curation_combined/` contains 62 rows, 236 patient–phenotype assessments
and 15 patients from four source cohorts: two Chinese KMT2D cohorts above, four
Korean MPV17 hepatocerebral cases in three families (PMID 27536553), and four Japanese
RPGRIP1 LCA cases in two families (PMID 31666973). Explicit absence and unknown are
separate; contradictory assessment claims were deferred. The Japanese cohort is a
follow-up and was counted once. These are selected, source-checked clinical cohorts,
not prevalence samples or independent clinician annotations.

The accepted rows collapse to 52 exact disorder–HPO pairs: 40 existing reference
associations and 12 new positive associations. Cohort ages and ascertainment differ;
pooling is exploratory, with source-deletion, family-balanced and truncating-KMT2D
sensitivities. No curated source shared an evaluation publication or reporting
group, stated care location or exact reported coding variant. Identity checks are
source-based and do not supply private-identifier linkage or HGVS equivalence.

## Clinical query eligibility and correction to the clinical benchmark

Five of the 18 frozen frequency queries were source-reported unsolved or undetermined.
One solved MPV17 neuropathy case fell outside the hepatocerebral target. The 12
remaining queries are KMT2D eight (three Saudi/five US), RPGRIP1 three Saudi and
MPV17 one Saudi, from three publications. Original US supplementary tables give
consistent KMT2D gene targets but discordant within-gene phenotype/age/variant links.
The prespecified strict identity sensitivity retains the concordant US ID 78,
leaving eight queries. No patient reassignment was guessed.

The whole Saudi clinical collection includes 1,601 solved, 313 undetermined, 661
unsolved and one unknown-status record. Restricting the existing nearest-case
clinical benchmark to source-reported solved diagnoses leaves 1,260 rankable Saudi
queries: top 1 is 2.78%, top 10 is 15.63%, median rank 143.5 and MRR 0.07663.
All 124 US records are source-reported solved, with 121 rankable queries unchanged
(top 1 1.65%, top 10 14.05%). Clinical records were never in that reference library,
so filtering query eligibility does not require recomputing reference scores.
Main text and Supplementary Table S15 use the solved-only counts; the full original
and restricted summaries are retained in `clinical_evaluation_audit/`.

## Frequency-ranking results

`frequency_benchmark/` evaluates the frozen rule in `docs/frequency_evaluation_plan.md`.
The filtered reference has 6,550 OMIM disease profiles mapped to 4,891 genes; 799 rows
citing evaluation publications or their known reporting groups were removed in all
arms. The Kabuki1-to-KDM6A mapping was removed in all arms. Query profiles, candidates
and gene-derived IC are fixed; weights enter the disease-to-query BMA direction.

| Reference | MRR | Median rank | Top 10 / 12 |
|---|---:|---:|---:|
| Unweighted | 0.004807 | 275.5 | 0 |
| Existing HPOA frequencies | 0.005040 | 307.0 | 0 |
| Updated frequencies, fixed support | 0.005050 | 305.0 | 0 |
| Added associations, unweighted | 0.004849 | 270.5 | 0 |
| Added associations and updated frequencies | 0.005027 | 300.5 | 0 |

The primary update improved four ranks, worsened one and left seven unchanged.
Two RPGRIP1 queries and one KMT2D query attained their maximum gene score through
another mapped disorder profile; their gene scores therefore did not respond to
the intended-disorder frequency update. Winning profiles and intended-profile
scores are included in the query-level outputs.
Paired MRR change is 0.0000101; a publication bootstrap interval is [0, 0.0000182],
but only three publications contribute and this interval is descriptive. All
category-bound, missing-weight, family-balanced, truncating-only, source-deletion
and strict-identity analyses are released; none yielded top-10 recovery. This is
an inconclusive test of diagnostic benefit for this weighting rule, with extremely
sparse clinical profiles and no prospective outcomes. It is retained in the
supplement, not promoted to the paper's central claim.

## Within-disorder provenance comparison

A further source review found four unrelated German RPGRIP1 LCA patients with
explicit yes/no assessments (PMID 30576320, Tables 1 and 3). This source was checked
after the first ranking run and was not retrospectively added to that frozen update.
`frequency_source_contrast/` compares nystagmus and cataract with the Japanese cohort:
nystagmus is 4/4 in each; cataract is 1/4 German versus 0/4 Japanese. The positive
German case is aged 47; the Japanese table gives ages 1–9 with later tests at 7–15.
Patient counts are descriptive; exact intervals use the six independent families
whose observed status is concordant for these traits. Study, population and care
country are perfectly confounded, while ages and variants differ. Thus source
denominators support a comparison but no adjusted population effect is identifiable.
This is Supplementary Table S24.

## Manuscript placement and verification

The main paper retains care-setting associations and sensitivity to phenotype
restriction as its argument. It now includes source-specific heterogeneity and
gene-level support limits. Frequency provenance and ranking remain exploratory
supplementary results. Tables S18–S24 and two figures document the extension;
new Methods specify extraction, planned human review, genetic matching, assessment
denominators, leakage exclusions, scoring and sensitivities. All pending clinical
validation is labelled. Reproduction commands are in `docs/extension_reproduction.md`.
