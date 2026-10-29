# Extension results: source, genetics and phenotype frequency

Work in progress. New outputs use the original HPO release (2026-01-08). Human review by Paul and Marwa is pending.

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

Both views will use the same source-neutral prompt and mapping pipeline. Existing abstract/full-text mapping results mix extraction and curation processes and are not being relabelled as this new controlled experiment. Human identity alignment, source-text extraction review and denominator verification remain outstanding.

## Pending completion

Gene-specific adjusted analyses, matched-source extraction and mapping, Paul/Marwa review and adjudication, verified frequency-denominator curation, independent frequency-aware ranking evaluation, and manuscript integration remain in progress. No manuscript conclusion has yet been changed on the basis of the unfinished extension.

## Gene-specific exploratory analysis

With at least 10 cases and three publications per regional group after matching, 342 gene–region contrasts were eligible within corpus, three after conditioning on author region, and six after matching disease identifier and allelic state. No care-country-conditioned contrast with stated origin met these support thresholds. No composition or individual-HPO-term test reached BH q < 0.05 across the complete tested families (4,999 permutations). This reflects limited overlap and multiplicity as well as uncertainty; it does not establish absence of genetic effects. The lowest composition q was 0.0684 for CYP17A1/EASIA test share, an exploratory candidate rather than a confirmed finding. The disease/allelic-state subset includes FBN1, KMT2D, RPGRIP1, SPINK5, TGFBR2 and UMOD; explicit variant matching and inheritance classification remain to be assessed.

The matching audit is in `gene_regional/support.csv`, all tested effects in `contrasts.csv`, and settings and test-family summaries in `summary.json`. No thresholds were relaxed after observing these results. The strong reduction in supported comparisons on conditioning is directly relevant to the scope of any genetics claim.

## Frequency evidence scope

Among count annotations, 88,951 have exactly one reference that is a PMID, 8,262 have multiple or mixed references including a PMID, and 591 have no PMID. A frequency supported by multiple references cannot be allocated to a particular study or population from HPOA alone. The audit now records this explicitly in `frequency_scope`; source cohorts and their denominators must be recovered before a population comparison or pooled frequency update.

Blinded reference packets for Paul and Marwa have been prepared for 20 publications (40 documents each). They contain empty annotation return sheets, source texts, instructions to distinguish unknown from assessed absence, and a separate adjudication template. Their status remains `prepared_not_reviewed`.

## Influence and exact-variant overlap

The deletion analysis recomputes each uncentered, within-gene/corpus Mantel–Haenszel test-share contrast after removing one gene. Largest changes (percentage points): MENA: HBB -0.14; SASIA: ATP7B +0.26; LATAM: CYP17A1 -0.24; AFR: HBB +0.94; EASIA: CYP17A1 -0.11. These are descriptive influence measures, not significance tests.

There are 36 exact variant/allelic-state profiles represented across regional groups, but 0 regional comparisons meet the original support thresholds after matching corpus, exact variant expressions and allelic state. Matching uses genomic expressions where available, otherwise coding expressions; variants have not been normalized across builds/transcripts. No relaxed-threshold confirmatory claim is made.
