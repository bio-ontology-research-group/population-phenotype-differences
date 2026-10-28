# Re-analysis results (2026-10-24)

Engine: `rank_engine.py` (vectorised Resnik BMA, extrinsic/intrinsic IC exactly as
`evaluation/similarity/scripts/phenotype_similarity_analysis.py` in the PAVS repo).
Reference: HPO release 2026-01-08 (hp.obo, genes_to_phenotype.txt with 5,193 genes,
phenotype.hpoa). Cases: `data/PAVS_phenopackets.json` (sources A,B,F,P = Saudi-clinical,
M = Saudi-literature, Q = US-clinical) and Phenopacket Store 0.1.26.
Validation: 299/300 sampled cases reproduce the repo's per-case rank exactly when
negated features are counted as present (`--keep-excluded`), Spearman 1.000
(`results/validate.csv`). Runtime: ~100 cases/s.

## 1. Negated features were counted as present in the published numbers

`phenotype_similarity_analysis.py` takes every `phenotypicFeatures[].type.id`
without checking `excluded`. In the Store, 92% of phenopackets carry excluded
features, typically the disease-typical features recorded as absent. Dropping
them (extrinsic Resnik, `results/baseline_summary.csv`):

| Cohort | n | terms/case (paper) | terms/case (present only) | Hits@1 (paper) | Hits@1 (present only) | MRR |
|---|---|---|---|---|---|---|
| Literature (Store) | 8,886 | 21.5 | 8.48 | 57.78 | 40.46 [39.4, 41.5] | 0.498 |
| Saudi-literature (M) | 1,018 | 3.8 | 3.79 | 8.35 | 8.35 [6.7, 10.0] | 0.129 |
| Saudi-clinical (A,B,F,P) | 1,948 | 3.8 | 3.84 | 1.39 | 1.39 [0.9, 1.9] | 0.038 |
| US-clinical (Q) | 122 | 3.8 | 3.82 | 1.64 | 1.64 [0.0, 4.1] | 0.034 |

Saudi and US cohorts are unaffected (44 negations in M, none elsewhere). The
8.5 present terms per Store case agrees with Danis et al. 2025 (8.2).

## 2. Leave-publication-out (LPO) circularity control (`results/lpo_percase.csv`)

For each Store case with source PMID p, every (gene, HPO term) pair in
genes_to_phenotype whose only supporting reference in phenotype.hpoa is PMID p is
removed before ranking (IC kept fixed). 5,124 of 8,886 cases have at least one
such term in their causative gene's profile.

| Subset | n | Hits@1 full -> LPO | MRR full -> LPO | Hits@10 full -> LPO |
|---|---|---|---|---|
| All Store cases | 8,886 | 40.46 -> 22.90 | 0.498 -> 0.317 | 66.6 -> 48.1 |
| PMID in HPOA | 6,025 | 53.28 -> 27.39 | 0.632 -> 0.365 | 80.7 -> 53.4 |
| PMID solely supports >=1 term of the gene | 5,124 | 56.28 -> 25.84 | 0.653 -> 0.339 | 81.9 -> 49.8 |
| PMID not in HPOA | 2,861 | 13.46 (unchanged) | 0.214 | 36.8 |

Roughly half of the top-1 performance on HPOA-sourced cases is attributable to
the case's own publication being in the reference. The residual gap between
LPO-corrected HPOA cases (27%) and non-HPOA cases (13%) is consistent with
indirect circularity (OMIM-referenced rows citing the same reports, companion
papers) and with HPOA-covered diseases being better characterised; it cannot be
separated with this data.

## 3. Depth does not explain the literature-vs-clinical gap (`results/depth_percase.csv`, `results/matched_summary.csv`)

Store cases were subsampled to term counts drawn from the Saudi-clinical
distribution (mean 3.84; 3 replicates), either at random or keeping the
highest-IC terms.

| Subset (extrinsic Resnik) | Hits@1 | Hits@10 | MRR |
|---|---|---|---|
| Store, full depth | 40.46 | 66.6 | 0.498 |
| Store, random subset at clinical depth | 21.24 | 45.7 | 0.297 |
| Store, top-IC subset at clinical depth | 31.85 | 59.9 | 0.415 |
| Store, PMID not in HPOA, random subset at clinical depth | 6.16 [5.7, 6.7] | 21.6 | 0.114 |
| Store, PMID not in HPOA, top-IC subset at clinical depth | 9.48 [8.9, 10.1] | 29.7 | 0.162 |
| Saudi-literature as recorded (3.79 terms) | 8.35 [6.7, 10.0] | 21.9 | 0.129 |
| Saudi-clinical as recorded (3.84 terms) | 1.39 [0.9, 1.9] | 8.0 | 0.038 |
| US-clinical as recorded (3.82 terms) | 1.64 [0.0, 4.1] | 5.7 | 0.034 |

At matched depth and with circularity removed, worldwide literature cases
(6-9% Hits@1) and Saudi literature cases (8%) perform alike, and both remain
4-7 times above routine clinical records (1.4-1.6%). Term choice (which
features a curator extracts from a diagnosis-aware publication versus which a
clinician writes on a referral form), not term count and not circularity,
carries the remaining difference. Mean IC: Store random subset 3.46, Saudi
literature 3.01, Saudi clinical 2.45.

## 4. Source distinguishability: population is not identifiable (`results/source_distinguishability.csv`)

Logistic regression on ancestor-propagated HPO vectors, 5-fold CV AUC, balanced
samples; "shared-gene" restricts to cases whose causative gene occurs in both
groups (controls disease mix).

| Pair | What differs | AUC all | AUC shared-gene |
|---|---|---|---|
| F vs P (KFSHRC 2017 vs 2019) | same institution, same population | 0.59 | 0.59 |
| A+B vs P (KAMC vs KFSHRC) | institution, same population | 0.79 | 0.76 |
| A+B vs F | institution, same population | 0.76 | 0.72 |
| P vs Q (KFSHRC vs Michigan) | institution AND population | 0.82 | 0.82 |
| F vs Q | institution AND population | 0.85 | 0.75 |
| A+B vs Q | institution AND population | 0.69 | (n<40) |
| M vs P (Saudi literature vs Saudi clinical) | source type, curator, same population | 0.91 | 0.87 |
| M vs Store (Saudi vs world literature) | population AND curator | 0.95 | 0.92 |
| P vs Store | everything | 0.96 | 0.94 |
| placebo: random split within P / within Store | nothing | 0.52 / 0.50 | |

Two Saudi institutions are as distinguishable from each other (0.72-0.79) as
Saudi clinical records are from US clinical records (0.69-0.85). Source type
and curator (0.87-0.95) dominate everything. With population perfectly nested
in source, no population effect can be estimated from PAVS.

## Files

- `results/baseline_percase.csv`, `results/baseline_summary.csv`, `results/baseline_ci.csv`
- `results/lpo_percase.csv`
- `results/depth_percase.csv`, `results/matched_summary.csv`
- `results/source_distinguishability.csv`
- `results/validate.csv`

## Phenopacket Store stratified by patient origin (LLM-extracted labels, 2026-10-24)

Pipeline: `store_origin_extract.py fetch` (Europe PMC full text when open, abstract otherwise; 64% of papers abstract-only by token count) then `store_origin_extract.py llm` (Qwen3.8-27B-FP8 on unimatrix, temperature 0, JSON schema) over the union of 0.1.26 and 0.1.27 source PMIDs; 1,731 of 1,733 papers returned a result. Labels: `origin_analyses.py make-labels-llm` -> `results/labels_llm.tsv` (majority region of patient_origins weighted by n_patients; recruitment country fallback for 137 papers; MIXED if the majority <70%). Regions: MENA (Arab countries, Turkey, Iran, Pakistan, Afghanistan), EASIA, EUR (Europe, North America, Australia/NZ, Israel, Russia), OTHER.

Paper-level labels (1,731 papers): UNKNOWN 909, EUR 317, EASIA 248, MENA 151, OTHER 65, MIXED 41. Origin stated explicitly in 690 papers; among the 822 labelled papers confidence high 746 / medium 71 / low 5. Agreement with the first-author-country proxy where both are known: 574/669 (85.8%); the main disagreement is European-affiliated authors reporting MENA patients (54 papers), i.e. the proxy undercounts MENA. Spot checks of 13 MENA/EUR papers: evidence quotes support the label in all.

### Ranking by region (Store, extrinsic IC; `results/origin_rank_llm.csv`, `results/origin_rank_llm_ci.csv`)

| region | cases | papers | cases with PMID in HPOA | Hits@1 full [95% CI] | Hits@1 LPO [95% CI] | Hits@1, PMID not in HPOA (n) |
|---|---|---|---|---|---|---|
| EUR | 2,590 | 259 | 57.7% | 31.1 [29.3, 32.9] | 20.4 [18.9, 22.0] | 8.6 [7.0, 10.3] (1,095) |
| EASIA | 784 | 227 | 57.7% | 40.2 [36.7, 43.6] | 27.6 [24.5, 30.6] | 23.2 [18.7, 27.7] (332) |
| MENA | 500 | 127 | 75.0% | 52.0 [47.6, 56.0] | 23.0 [19.4, 26.6] | 15.2 [9.6, 21.6] (125) |
| OTHER | 228 | 58 | 69.7% | 53.5 [46.5, 59.7] | 32.5 [26.3, 38.6] | 33.3 [21.7, 43.5] (69) |
| MIXED | 201 | 31 | 85.6% | 52.7 [45.8, 59.7] | 26.4 [20.4, 32.8] | 17.2 [3.4, 31.0] (29) |
| UNKNOWN | 4,583 | 791 | 73.6% | 43.4 [41.9, 44.8] | 22.9 [21.7, 24.0] | 13.8 [11.8, 15.8] (1,211) |

The apparent MENA advantage under the full reference (52.0 vs 31.1) is mostly circularity: MENA publications are more often cited in HPOA (75% vs 58%), and under LPO the MENA and EUR intervals overlap. Region also differs in disease/gene mix, so these marginal differences are not population effects.

### Within-gene contrasts (curator fixed)

Genes with >=5 MENA and >=5 EUR cases, each from >=2 papers: 5 (47 MENA cases); EASIA vs EUR: 7. Paper-level permutation test of profile divergence (JSD, 2,000 permutations): ATP13A2 p=0.47, CLDN16 p=0.08, ITPR1 p=0.48, MPV17 p=0.53, RPGRIP1 p=0.90; none significant (`results/origin_profile_llm.csv`).

Population-specific reference test (leave-one-paper-out; gene profile replaced by union of terms from other MENA papers vs a size-matched EUR sample; `results/origin_popref_llm.csv`, threshold 5 cases):

| held-out | cases (genes) | Hits@1 MENA-ref | Hits@1 EUR-ref | Hits@1 HPOA |
|---|---|---|---|---|
| MENA | 55 (6) | 18.2 | 36.4 | 49.1 |
| EUR | 39 (5) | 51.3 | 35.9 | 51.3 |

The same-region reference was better in 27% and worse in 58% of held-out MENA cases; for held-out EUR cases better in 23%, worse in 31%. Results are dominated by single genes (RPGRIP1: 16 MENA cases from 3 papers). The small reference profiles built from a handful of papers underperform the full HPOA in both directions; no advantage of a population-matched reference is detectable. With 5-7 testable genes the test has little power: it rules out large effects for these genes only.

### Pooled within-gene model (`pooled_origin_model.py`, `results/origin_pooled_model.csv`; exploratory, 2026-10-25)

OLS with gene fixed effects, covariates log(n_terms), mean IC, HPOA citation; SEs clustered by publication. Only genes with cases from >=2 regions contribute: 1,838 cases, 362 papers, 79 genes (MENA 148, EASIA 383, EUR 1,120).

- MENA vs EUR, LPO: log2 rank +0.18 [-0.45, +0.80]; Hits@1 -9.1 points [-22.5, +4.2]. Not significant; the CI excludes only large effects.
- EASIA vs EUR, LPO: log2 rank +0.44 [+0.00, +0.88], Hits@1 -7.2 [-14.3, -0.2] (p~0.05).
- Adding first-author region (affiliation proxy): MENA patient Hits@1 -11.8 [-23.0, -0.7] (p=0.038, one of ~20 tests, log-rank not significant); EASIA author effect -16.4 [-21.5, -11.2] (p=5e-10) while the EASIA patient effect turns +9.9. Recorder/author effect is larger than patient effect where they can be separated.
- MENA patients reported by European vs MENA authors (no gene FE, disease mix confounded): log2 rank +4.0, not interpretable.

Power: SE for MENA Hits@1 ~6.8 points with 148 informative cases; detecting a 5-point difference (SE ~2.5) needs roughly 7x more within-gene MENA cases (~1,000).

Europe PMC (title/abstract nationality terms + sequencing case-report terms, 2026-10-25): MENA 8,738 papers (5,680 open access), European 16,160 (10,911 OA). Per test gene, 30-160 OA case-report papers vs 5-10 Store papers.

## Case-to-case design (2026-10-25)

Per Robert: population differences must be tested by comparing cases with cases, not cases with an aggregated (OMIM/HPOA) disease profile. `case_pairs.py`:
- similarity: symmetric BMA Resnik between the observed terms of two cases; IC from the case corpus itself (no HPOA); hubness-normalised (per-case z-scores averaged over both cases); pairs from the same publication never compared;
- `pairs`: within gene, mean similarity of same-region pairs minus cross-region pairs, depth-adjusted, averaged over genes; three permutation nulls over publications: within gene (combined population+recorder effect), within (gene, author region) (patient effect beyond recorder), within (gene, patient region) (author effect beyond population);
- `retrieval`: nearest-case gene ranking (max similarity to library cases of each gene), the query's own gene represented by k cases from the same vs the other region; paired per query, bootstrap over publications.

Store results (LLM patient labels, first-author proxy as recorder; 1,000 permutations):

| contrast | pairs (genes) | same-minus-cross | combined p | patient excess (p) | author excess (p) |
|---|---|---|---|---|---|
| MENA vs EUR | 7,299 (27; 12 with both pair types) | -0.13 | 0.60 | -0.02 (0.82) | -0.12 (0.26) |
| EASIA vs EUR | 36,133 (45; 32) | +0.17 | 0.26 | -0.01 (0.75) | +0.01 (0.77) |

Retrieval, MENA queries (n=64, 21 publications, 9 genes): log2(rank_other/rank_same) -0.33 [-1.06, +0.77]; EASIA queries (n=276): +0.36 [-0.06, +0.78]. Raw-similarity retrieval was biased by library-case depth (hubs); z-normalisation reduces but does not remove it.
In the Store patient and author region are nearly collinear, so the two nulls preserve each other's structure: no separation is possible; the combined effect is not significant either.

## Expansion: PMC-Patients + Qwen (running 2026-10-25)

PubCaseFinder checked: its API lists case-report PMIDs per MONDO disease but returns no phenotypes; its per-article HPO annotations (PubCasesHPO on PubAnnotation, dictionary-based, abstracts only) stop around 2017 (0/36 Store papers from 2018 on covered, 7/24 older) with 2-3 terms per abstract: unusable as phenotype source. RareArena: rephrased, Orphanet diagnoses, no genes/HPO. Used instead: PMC-Patients V2 (Zhao et al. 2023; 250,294 patient summaries from 210,069 PMC OA articles, CC BY-NC-SA). Filter (HPOA gene symbol + variant vocabulary): 14,500 candidates / 12,979 articles; first-author regions EUR 5,626, EASIA 4,737, MENA 781, OTHER 690 articles.
- `pmcp_extract.py patient` (Qwen): germline Mendelian triage, causal gene(s), observed/excluded phenotype phrases with proposed HPO label, origin in summary, consanguinity. Pilot: ~47% Mendelian; somatic/aneuploidy/inconclusive correctly rejected.
- `pmcp_extract.py map`: unmatched labels -> token-retrieval candidates -> Qwen picks (no more specific than the finding). Pilot coverage: exact 64.7%, Qwen 28.2%, unmapped 6.9%.
- article-level origin: `store_origin_extract.py` fetch+llm on Mendelian articles.
- `pmcp_to_cases.py` -> case table; `validate_extraction.py` compares with Store curation on shared publications, by region.
Jobs: unimatrix 9666 (extract + origin), 9669 (map, after 9666). Store unknown-origin recovery: 9664 (OpenAlex OA copies, ~7% retrievable by script; publishers block).

### Store case-to-case, v2 origin labels (after full-text recovery, 2026-10-26)

z-similarity, publication-level permutations (2,000), pairs within the same publication excluded.

| Contrast | Pairs (genes) | Combined excess, p | Patient-beyond-author excess, p | Author-beyond-patient excess, p |
|---|---|---|---|---|
| MENA vs EUR | 8,801 (30; 16 with both types) | -0.29, 0.32 | -0.03, 0.87 | -0.10, 0.31 |
| EASIA vs EUR | 41,436 (54; 38 with both types) | +0.08, 0.57 | +0.01, 0.89 | -0.05, 0.70 |

Retrieval (k same-region vs other-region library cases for the query's gene): MENA queries (74, 25 pubs) +0.08 log2 [-0.53, +0.70]; EUR queries -0.10 [-0.62, +0.25] (MENA contrast); EASIA queries +0.84 [+0.34, +1.39], but EUR queries in the same contrast -0.50 [-0.96, +0.03]. The EASIA asymmetry (EASIA library cases help both EASIA and EUR queries) points to a library-composition effect, not population matching; no evidence for a population-matched advantage in the Store.

## PMC-Patients case-to-case analysis (2026-10-26)

Corpus: PMC-Patients V2 genetic candidates (14,500 summaries) run through Qwen3.8-27B on unimatrix (pmcp_extract.py patient + map; article origin via full text). 7,013 Mendelian germline, 6,828 single-gene cases (1,652 genes, 6,027 publications; results/cases_pmcp.tsv). Regions: EASIA 2,139, EUR 1,945, MENA 600, OTHER 540, MIXED 59, UNKNOWN 1,545. Phrase mapping: exact 67.4%, Qwen candidate choice 25.5%, fuzzy 0.2%, unmapped 6.8%.

Patient x author region: MENA patients with EUR-affiliated authors 113 cases (MENA-authored 371); EASIA patients almost always EASIA-authored (2,008/2,109), so EASIA cannot separate patient from author.

Validation against Store curation (51 shared publications; results/pmcp_validation.csv): precision exact 0.44 / lenient 0.61, recall exact 0.56 / lenient 0.77; Qwen records more terms (15.7 vs 11.0). Similar across regions (MENA lenient P/R 0.77/0.75, EASIA 0.63/0.78, EUR 0.61/0.69).

All tests below exclude pairs and library cases from the same publication or from publications sharing any author (results/pmcp_authors.tsv; shared-author pairs were 1-4% of within-gene pairs).

Pair test (mean within-gene similarity, z-similarity, 2,000 publication-level permutations):
- MENA vs EUR: 3,661 pairs, 136 genes. Combined excess -0.14 (p=0.14); MENA-MENA pairs are the least similar (1.87 vs 2.17-2.18).
- EASIA vs EUR: 17,034 pairs, 318 genes. Combined excess +0.12 (p=0.038); patient beyond author +0.004 (p=0.72); author beyond patient -0.03 (p=0.09). Not separable from author region.

Retrieval test. Naive per-region averages are dominated by a library main effect (MENA cases are better prototypes for EUR queries too). Population matching is the interaction: within genes with queries from both regions, mean of (d_r1 + d_r2)/2 where d = log2(E[rank other-region library] / E[rank same-region library]), expectation exact over k-subsets (retrieval_perm.py). This statistic is biased upward under the null (selection of genes/libraries), so it is only interpretable against the permutation null (1,000 perms; labels permuted between publications within gene, or within gene x author region).

| Contrast | Genes | Statistic | Null (gene) excess, p | Null (gene x author) excess, p |
|---|---|---|---|---|
| MENA vs EUR, mean, k=5 | 44 | +0.88 | +0.50, 0.007 | +0.25, 0.027 |
| MENA vs EUR, median, k=5 | 44 | +0.73 | +0.48, 0.012 | +0.38, 0.028 |
| MENA vs EUR, raw similarity | 44 | +0.85 | +0.46, 0.008 | +0.22, 0.035 |
| MENA vs EUR, k=3 | 44 | +0.85 | +0.48, 0.008 | +0.22, 0.040 |
| MENA vs EUR, k=1 | 44 | +0.71 | +0.29, 0.050 | +0.11, 0.17 |
| EASIA vs EUR, mean | 152 | +0.44 | +0.18, 0.023 | -0.03, 0.79 |
| EASIA vs EUR, median | 152 | +0.26 | +0.15, 0.029 | -0.03, 0.79 |
| Store MENA vs EUR (replication) | 7 | +0.26 | -0.06, 0.58 | -0.12, 0.69 |
| Store EASIA vs EUR | 20 | +0.20 | -0.01, 0.51 | +0.02, 0.34 |

Descriptive: MENA queries reported by EUR-affiliated authors (n=41) have mean d +1.10, MENA-authored +0.55, so the MENA advantage is not carried by MENA authors. Per gene, 26/44 genes positive; most genes contribute 2 MENA and 2-4 EUR queries; the mean is pulled by a few genes (NOD2, XPC, ABCD1, VPS13B, AAAS), but the median-based test gives the same conclusion.

Reading: one-sided, two contrasts, small gene count; MENA p-values 0.03-0.04 with author fixed would not survive a strict multiplicity correction across all variants tried. EASIA: any effect is fully explained by author region. MENA: modest evidence that a patient's nearest same-gene case is more informative when it comes from a MENA patient, beyond the reporting authors' region, while average similarity is not higher; consistent with founder alleles giving near-identical presentations in a subset of families (not tested). Store does not replicate (7 genes; underpowered).

## Is the population difference phenotype-, gene- or disease-specific? (2026-10-26)

Script: phenotype_specificity.py (termdb / terms / genes). FUNC 0.4.8 (OntoFUNC) termdb built from hp.obo 2026-01-08 (19,408 terms). PMC-Patients cases, same-publication and shared-author pairs excluded for the gene-level test.

Term level (MH-weighted within-gene difference in proportion of cases with each propagated HPO term, publication-level permutations, 2,000):
- MENA vs EUR, gene strata: 968 cases (319 MENA), 185 genes, 977 terms tested. 75 terms p<0.05 (49 expected); none FDR<0.1 or FWER<0.05 (top: eosinophilia +3% vs 0%, immune-system physiology -9.5 points; q=0.12). Global sum z^2 p=0.013. No depth difference (log #terms -0.05, null [-0.09, +0.09]).
- Gene x author strata (separates patient from author): 315 cases (92 MENA), 60 genes; global p=0.48; nothing significant.
- OntoFUNC func_binom on gene-balanced annotation counts (5 random balancings, root HP:0000118): global tests p 0.07-0.59; smallest FWER 0.073. Nominal and consistent across balancings: MENA more abdominal symptoms, skin hyperpigmentation, tendon/contracture terms; EUR more arrhythmia, lymphocyte, motor delay, immune physiology terms. Not significant after correction.

Gene level (per-gene same-region minus cross-region depth-adjusted pair similarity, z against publication-level permutation null):
- MENA vs EUR: 122 genes; mean z -0.07 (null [-0.18, +0.18]); heterogeneity sum z^2 p=0.39 (gene x author strata: 44 genes, p=0.17). No gene individually significant (top AAAS z +2.7, ALPL +2.4, MEFV +1.9). AR-only vs AD-only genes: Mann-Whitney p=0.34.
- EASIA vs EUR: 300 genes; mean z +0.19 (null [-0.11, +0.11]); heterogeneity p=0.052; LMNA (z +4.0), HNF1B (+3.3) top. Gene x author strata: 32 genes, p=0.39. AR vs AD p=0.61.
- OntoFUNC func_wilcoxon (gene -> HPOA phenotype and inheritance terms, value = per-gene z, 1,000 random sets):
  - MENA, phenotypic abnormality root: global high-ranks p=0.018. Abnormal corneal epithelium morphology (11 genes) FWER 0.008, Keratitis (10 genes) FWER 0.027; genes AAAS, MEFV, IKBKG, NOD2, IFNGR1, XPC, WAS, LMNA, FBN1, AIRE, NTRK1, mostly immune/autoinflammatory. Broad immune categories not enriched (Abnormality of the immune system p_high 0.43, Autoimmunity 0.18). Mode of inheritance root: global p 0.11 / 1.
  - EASIA: global p 0.87 / 0.63; Abnormality of head or neck p_high 0.003, FWER 0.95.

Reading: the data do not support a broad shift affecting all phenotypes or all genes, nor an organ-system or inheritance-mode pattern; the within-gene phenotype-profile difference is detectable only globally and only without the author control. The one category hit (a small cluster of immune/autoinflammatory genes annotated to keratitis, led by AAAS and MEFV) overlaps the genes driving the retrieval interaction (AAAS, NOD2, XPC), but rests on ~10 genes with 2-4 MENA cases each; hypothesis-generating only.

## OntoFUNC with ELK, and the merged multi-provenance model (2026-10-26)

OntoFUNC: FUNC termdb rebuilt from hp.owl 2026-01-08 with ~/Public/software/ontofunc/MakeFuncStructure.groovy (ELK 0.4.2; copy patched in the scratchpad only to replace CliBuilder for Groovy 4). The inferred closure adds 30,956 HP-HP ancestor links to the asserted 197,919 (+15.6%), removes none. Re-running everything on it changes nothing material: term test global p=0.016 (gene strata), 0.41 (gene x author); func_binom smallest FWER 0.069 to 0.45 over balancings; func_wilcoxon MENA high ranks global p=0.039, Abnormal corneal epithelium morphology FWER 0.009, Keratitis 0.022; EASIA nothing.

Merged case set (merged_model.py; similarity over all 19,827 cases, IC from the combined corpus): Store 9,517, PMC-Patients 6,828, PAVS Saudi clinical (A/B/F/P) 2,576, PAVS Saudi literature (M) 1,278, PAVS US clinical (Q) 124. Pairs and libraries exclude same-publication and shared-author publications (author lists for 7,654 PMIDs; the clinical cohorts' cohort papers are used for their authors, so literature from the same groups is excluded).

Pair model (37,375 within-gene MENA/EUR pairs, 498 genes; OLS with gene FE, depth, provenance-pair FE, author-region terms; 1,000 publication-level permutations):
- Provenance dominates: provenance-pair effects span ~2.4 z-units (Store-Store +2.56, PMCP-PMCP +1.52, PAVS clinical vs Store +0.70, vs PMCP +0.34, relative to the baseline pair class); same author region +0.47; patient region same-region -0.20.
- Same patient region beyond provenance and author region: all pairs excess -0.01 (p=0.89); literature-literature -0.03 (p=0.73); Saudi-clinical vs literature +0.05 (p=0.42). Without the author control the same-region coefficient is negative (p=0.003), i.e. what looks like a patient-region effect is author-region structure.
- Calibration: the same publication and gene recorded by two processes (Store vs PMCP) has mean z 4.97, versus 2.78 for different Store publications of the same gene.

Retrieval (retrieval_perm.py, library restricted to literature corpora, permutation strata include corpus):
- Saudi clinical queries (887 queries, 121 genes): MENA literature cases rank them better, gene-averaged log2 rank ratio +0.68; excess over gene null +0.34 (p=0.003), beyond author region +0.09 (p=0.032). Without matched EUR clinical queries this statistic includes any general prototype advantage of MENA literature cases.
- US clinical queries (68, 20 genes): -0.30, no preference (p=0.61 / 0.86). Clinical interaction (10 genes with both): nothing.
- Literature queries, all three literature corpora (64 genes): excess +0.33 (p=0.021), beyond author +0.13 (p=0.12): weaker than PMC-Patients alone (p=0.027).

Reading: pooling sources does not strengthen the population signal. Recording provenance and the region of the reporting group explain far more case-to-case similarity than patient origin; the only effect beyond author region in the merged data is a small advantage (+0.09 log2 rank) of MENA literature cases as nearest neighbours for Saudi clinical cases.

## Literature on the leads (2026-10-26; Sonnet agents, PubMed-verified, abstracts read)

Files: literature/leads_immune_genes.{md,bib} (48 refs), leads_founder_genes.{md,bib} (58), leads_easia_and_recording.{md,bib} (29).
- Documented population differences exist only for MEFV (Eurofever: eastern-Mediterranean residents more severe than migrants of the same ancestry in Europe; M694V homozygosity and amyloidosis) and AIRE (Iranian Jews vs Finns: rare candidiasis, no keratopathy), i.e. allele plus environment; the AIRE ocular and MEFV inflammation directions are opposite to our leads.
- Founder alleles are well documented for AAAS (c.1331+1G>A, Maghreb) and XPC (c.1643_1644delTG, Maghreb; Libyan founder with cutaneous-only disease), but the literature reports no genotype-phenotype correlation for AAAS and no evidence that founder populations have more homogeneous phenotypes. ALPL Saudi series severe and consanguineous (zygosity/ascertainment); ABCD1 and VPS13B: nothing population-specific.
- NOD2, IKBKG, IFNGR1, WAS: no population comparisons; plausible confounders are disease stage at report (neonatal eosinophilia in incontinentia pigmenti), BCG vaccination policy, consanguinity shifting to complete recessive forms, and mixing of distinct disease entities per gene (Blau/Crohn/Yao).
- EASIA leads (LMNA, HNF1B, CLCN5): presentation follows the recruiting specialty and screening programme (Japanese school urine screening for Dent disease; renal vs diabetes clinics for HNF1B; cardiology vs neuromuscular for LMNA): ascertainment, consistent with our finding that the EASIA effect is author/reporting region.
- No published study quantifies recording provenance against patient origin; the provenance result appears new.

## Provenance vs population, with reporting groups as replication units (2026-10-26)

Estimand: a population effect is any difference between patient populations that replicates across independent reporting groups once the curation layer is held fixed; genetic or not (screening, vaccination, referral all count). Layers: curation (corpus: Store curators, Qwen/PMC-Patients, PAVS curators, PAVS clinical records) held fixed by corpus FE / within-corpus permutation; reporting group = senior (last) author (5,991 groups; largest 37 publications; initials-only names may merge labs, which makes tests conservative); labels permuted between groups. 'Total' = author region not held fixed (includes local practice); 'patient-carried' = permutation within author region (rests on patients reported from outside their region). Script population_model.py; retrieval_perm.py --groups.

Decomposition of within-gene pair similarity (z units, gene FE, all layers in one OLS, 200 group-cluster bootstraps):
- MENA vs EUR (62,756 pairs, 503 genes): curation-pair effects relative to Store-Store -1.2 to -3.0 (e.g. PMCP-PMCP -1.45 [-2.40, -0.53], PAVS clinical-Store -2.44 [-3.14, -1.39], PAVS clinical-clinical -3.02); same group +0.32 [-0.05, +1.02]; same author region +0.11 [-0.59, +0.45]; same population +0.00 [-0.25, +0.25]. Population terms: total excess -0.004 (p=0.96), patient-carried -0.02 (p=0.77).
- EASIA vs EUR (80,382 pairs, 497 genes): curation -1.4 to -2.4; same group +0.31; same population: total excess +0.10 (p=0.12), n_EASIA +0.14 (p=0.07); patient-carried +0.06 (p=0.16).

Term frequencies (MH within gene x corpus, inferred HPO closure, 2,000 group permutations):
- EASIA vs EUR, total: 4,510 cases (1,869 EASIA), 2,299 groups; 344 of 2,288 terms p<0.05 (114 expected), 95 FDR<0.1, 10 FWER<0.05; global p=0.0005. Largest: abnormal circulating nitrogen compound concentration +7.2 points (FWER 0.001), abdominal distention +2.0, cardiomegaly +1.4; obesity -1.9, dehydration -1.0, ataxia -3.1 points. Replication across curation processes: Spearman of term z Store vs PMCP +0.11 (p=0.001); signals nominal in both: circulating metabolite/nitrogen/organic compounds and azotemia higher in EASIA, ataxia/coordination and memory B cell terms lower. Branch level under group permutation: blood and blood-forming tissues higher (mean z +0.38, p=0.014), nervous system lower (-0.37, p=0.020), neither FWER<0.05 over 21 branches (FUNC's term-level random sets overstated these branches because nested terms are not independent). Patient-carried: 441 cases (112 EASIA), global p=0.18.
- MENA vs EUR, total: 1,475 cases (440 MENA), 935 groups; 84 of 1,248 p<0.05 (62 expected), 1 FWER<0.05 (abnormality of immune system physiology -8.5 points, FWER 0.027); global p=0.063. No replication across curation processes (Spearman -0.14; no term nominal in both). Branch: constitutional symptoms lower (p=0.008, FWER 0.11). Patient-carried: 614 cases (154 MENA), global p=0.90.

Nearest-case retrieval (library = literature corpora, group units):
- MENA literature queries: total excess +0.33 log2 (p=0.014); patient-carried +0.13 (p=0.13).
- Saudi clinical queries: total +0.33 (p=0.001); patient-carried +0.09 (p=0.035).
- EASIA literature queries: total +0.15 (p=0.035); patient-carried -0.01 (p=0.67).

Reading: population effects in the total sense exist and replicate across independent groups: EASIA vs EUR in which phenotypes are recorded (lab/metabolite findings more, neurological and immune less), and MENA vs EUR in nearest-case retrieval, including for Saudi clinical cases. They are an order of magnitude smaller than curation provenance (1-3 z units) and about the size of the reporting-group effect (+0.3). Only a small part travels with the patient when the reporting region is held fixed (Saudi clinical retrieval +0.09, p=0.035); the rest sits with the reporting region, where local practice (a population effect) and regional writing conventions (provenance) cannot be separated with ~110-150 patients reported from outside their region.

## Extraction and normalisation audit (2026-10-26)

Why the Store vs PMC-Patients term correlation is only 0.11: (1) only 14 genes have EASIA and EUR cases in both corpora (3 for MENA), so the two estimates mostly concern different genes (on the 14 shared genes: 0.05, n=180 terms); (2) split-half reliability of term-level EASIA-EUR differences (random halves of reporting groups, same process) is 0.27-0.31 within PMC-Patients and -0.11 to +0.14 within the Store (~400 EASIA cases in 56 gene strata). The attainable cross-corpus correlation is therefore ~0.15-0.2 even with identical recording; 0.11 is noise-limited, mainly by the Store.

Direct audit: 60 PMC-Patients cases (20 per region, shuffled into three region-blind batches), 887 extracted phrases judged against the source text by Sonnet agents (scratchpad/audit/batch*_judged.json).

| | EASIA | EUR | MENA | all |
|---|---|---|---|---|
| extraction precision | 98.5% | 97.2% | 99.3% | 98.4% |
| recall (vs findings in text) | 92.0% | 92.3% | 87.1% | 90.4% |
| mapping exact | 74.6% | 72.9% | 79.9% | 75.9% |
| mapping broader but correct | 8.1% | 11.7% | 8.1% | 9.0% |
| too specific | 5.4% | 2.5% | 2.7% | 3.7% |
| wrong | 5.7% | 4.6% | 6.4% | 5.6% |
| unmapped | 6.3% | 8.3% | 3.0% | 5.7% |

Recall by how a finding is established: physiological test 97%, clinical exam 94%, imaging 92%, symptom/history 91%, pathology 85%, laboratory 83%. Test-requiring findings: recall EASIA 90.8%, EUR 92.4%, MENA 79.8% (bedside/history 91-93% in all); share of extracted findings that are tests: EASIA 41.5%, EUR 30.4%, MENA 34.6%.

Mapping by route: exact label/synonym match (626 phrases) 85% exact + 7% broader, 8% errors; Qwen choice among 15 TF-IDF candidates (194) 68% exact + 17% broader, 11% wrong; unmapped 50 (47 had a fitting term). Weakest stage: candidate retrieval/selection for non-lexical phrases.

Reading: extraction is good and errors are not strongly region-dependent, except lower recall of lab findings in MENA cases (would bias MENA-EUR lab comparisons towards EUR; does not affect EASIA-EUR, where lab recall is equal and EASIA reports contain more test-based findings). Mapping errors (~15% of phrases wrong/too specific/unmapped) are non-differential noise that attenuates population effects; improving the candidate stage is the main lever.

### HPO normalisation: RAG-HPO-style mapper (2026-10-26)

Literature (literature/hpo_mapping_methods.md): RAG-HPO (Garcia et al., Genome Med 2025; MIT) fits a local OpenAI-compatible endpoint; dense retrieval over HPO labels/synonyms/definitions + LLM choice among candidates. Ported as hpo_map2.py + qwen_tasks.py select: index of HPO labels, EXACT/RELATED/NARROW synonyms, first definition sentence and RAG-HPO's 3,079 add-on phrases (62,925 entries); BGE-small (bge-small-en-v1.5) embeddings; candidates = dense top-15 for the verbatim phrase and for Qwen's proposed label, exact label/synonym hits, plus parents of the top hits; verbatim text that equals a label/exact synonym is accepted directly, otherwise Qwen chooses one listed candidate (same or closest more general; never more specific) or null.
- Retrieval recall of the audited gold term (859 phrases): BGE-small @1 0.81, @5 0.92, @15 0.95 (text+label); SapBERT @15 0.94.
- Against the audit gold (biased towards the old mapper, whose exact-judged outputs define the gold): old 84.6% acceptable, new 80.9%.
- Blind pairwise judgement of the 254 phrases where the mappers disagree (Sonnet, order randomised): acceptable old 51.6%, new 84.3%; too specific 54 vs 21, wrong 39 vs 9, missed 30 vs 10; preferred new 212, old 32, tie 10. New wins also where the old mapper had an exact label match (128 vs 15): Qwen's proposed label often matched a wrong HPO label exactly. Adopted for all 79,312 unique PMC-Patients phrases.

## Ancestry vs place of care (PMC-Patients, qwen_tasks.py setting; 2026-10-26)

6,267 of 7,013 Mendelian summaries parsed (746 to retry). Care basis: our_institution 4,304, explicit 1,458, none 338; migrant yes 203. Cases (cases_pmcp_setting.tsv): MENA ancestry treated in MENA 236, MENA ancestry treated in EUR/N America 103, ancestry not stated treated in MENA 289, EUR or unstated ancestry treated in EUR 2,335.
- Ancestry at fixed place of care (MENA-ancestry vs EUR/local patients, both treated in EUR; within gene, group permutation): 296 cases (69 MENA-ancestry), 55 genes: global term test p=0.40, nothing FDR<0.1; physiological-test findings -0.002 (p=0.98), share of test-based findings +0.001 (p=0.96).
- Place of care at fixed ancestry (MENA ancestry, MENA vs EUR care): 69 cases, 19 genes: global p=0.51; nothing.
Underpowered; the PubMed-abstract corpus extracts ancestry and place of care per patient and should enlarge these cells.

## Rerun on updated data (Store v3 labels incl. 174 Wiley TDM texts; PMC-Patients with RAG-HPO-style mapping; 2026-10-27)

Store v3 labels: unlabelled papers 775 -> 625; MENA papers 175 -> 192, EUR 389 -> 484 (cases MENA 693, EUR 4,411). PMC-Patients rebuilt with hpo_map2 (cases_pmcp_map2.tsv): 6,792 cases, 92.1% of phrases mapped.

| Test | MENA vs EUR total | MENA patient-carried | EASIA vs EUR total | EASIA patient-carried |
|---|---|---|---|---|
| Pair decomposition, same population excess | -0.03 (p=0.71) | 0.00 (p=1.0) | +0.16 (p=0.13); n_EASIA +0.19 (p=0.016) | +0.02 (p=0.81) |
| Term frequencies, global | p=0.013 (87/1,370 nominal; none FDR) | p=0.92 | p=0.0005 (382/2,449 nominal; 52 FDR<0.1; 11 FWER<0.05) | p=0.12 |
| Physiological-test findings (log count) | -0.100 (p=0.0015) | -0.169 (p=0.0065) | -0.034 (p=0.097) | -0.101 (p=0.24) |
| Share of test-based findings | -2.3 pts (p=0.17) | -2.7 (p=0.39) | +2.3 pts (p=0.041) | -0.9 (p=0.84) |
| Symptom/history findings | -0.075 (p=0.11) | -0.152 (p=0.063) | -0.049 (p=0.039) | +0.05 (p=0.54) |
| Retrieval, literature queries | +0.27 (p=0.031) | +0.04 (p=0.32) | +0.13 (p=0.057) | -0.02 (p=0.73) |
| Retrieval, Saudi clinical queries | +0.42 (p=0.001) | +0.11 (p=0.022) | | |

Old vs new mapping (total term tests): MENA global p 0.015 vs 0.013; EASIA 85 vs 52 FDR<0.1 terms, 9 vs 11 FWER<0.05, physiological-test MENA p 0.0020 vs 0.0015. The better mapping does not change the conclusions.

Most robust patient-level difference: MENA patients have fewer findings from instrument-based functional tests (ECG/EEG/EMG/audiometry/ophthalmological instruments), also when reported by groups outside the region (p=0.0065), not explained by extraction (audit recall of such findings ~97% in all regions).

## PubMed abstract corpus and v4 rerun (2026-10-27)

Corpus: 114,313 PubMed case reports with an HPOA gene; 64,102 candidate abstracts run through Qwen3.8-27B (qwen_tasks.py abstract, ORIX H200 x3; 4 failures); 42,173 abstracts with Mendelian patients. abs_to_cases.py -> results/cases_abs.tsv: 45,911 single-gene patients (corpus pubmed_abs), 94.3% of phrases mapped (RAG-style selection, ORIX; PMC-Patients selections reused). Region = ancestry if stated, else place of care: EUR 20,820, EASIA 12,076, MENA 4,732, OTHER 3,287, UNKNOWN 4,996. Setting: MENA ancestry cared for in MENA 1,685, in EUR 572; unstated ancestry cared for in MENA 2,285; EASIA ancestry in EASIA 5,267, in EUR 77. Sex stated 26,831; age in years 18,375, age band 23,006. No PMID overlap with Store, PMC-Patients or PAVS. Last-author groups: 46,709 papers. termclass: 2,627 new terms (24 unclassified).

Decomposition (IBEX, 66,198 cases, z-similarity, group permutations):
- MENA vs EUR (680,145 pairs, 1,920 genes): same population excess total +0.095 (p=0.002), patient-carried +0.053 (p=0.001); same group +0.72; Store-Store curation +3.0.
- EASIA vs EUR (1,140,700 pairs, 2,041 genes): total +0.076 (p=0.002); patient-carried +0.001 (p=0.92).

Term tests (store + pmcp_map2 + pubmed_abs; 2,000 group permutations; results/v4/):
- MENA vs EUR total: share of findings needing a test -5.5 points (p=0.0005); bedside/history findings +0.14 log-count (p=0.0005); 222 terms FDR<0.1, 32 FWER<0.05; integument, digestive, head/neck, blood, musculoskeletal, limbs higher (branch FWER<=0.002). Patient-carried: physiological-test findings -0.058 (p=0.0005), test share -2.5 (p=0.067), global p=0.15.
- EASIA vs EUR total: more findings of both kinds (test-based +0.034, bedside +0.059, both p=0.0005), test share unchanged; 415 FDR<0.1. Patient-carried: bedside +0.099 (p=0.013), global p=0.21.
- Ancestry vs place of care (pubmed_abs + PMC-Patients): MENA ancestry cared for in EUR vs European patients cared for in EUR: test share -0.3 points (p=0.87), global term p=0.0005 but no FDR term. MENA ancestry cared for in MENA vs cared for in EUR: test share -7.6 points (p=0.002), test-based findings -0.10 (p=0.021), bedside +0.11 (p=0.022). EASIA ancestry: cared for in EUR vs Europeans in EUR bedside +0.15 (p=0.015); in EASIA vs in EUR bedside -0.18 (p=0.020).
- Age and sex strata (pubmed_abs): the MENA lower test share holds in neonates/infants (-6.9, p=0.001), children (-5.3, p=0.003) and males (-7.2, p=0.0005); not significant in adolescents, adults or females; bedside excess in every stratum.

Reading: the lower share of test-based findings for MENA patients goes with being cared for in the MENA region, not with MENA ancestry: MENA patients cared for in Europe look like European patients on this outcome, and MENA patients cared for at home differ from the same ancestry cared for in Europe. This is a population (care-setting) difference in the broad sense, not a genetic one. The physiological-test deficit is the one outcome that also holds patient-carried.

Phenopackets: to_phenopackets.py -> results/phenopackets/ (66,198 packets; index.tsv; phenopackets.zip). All pass a strict GA4GH v2 protobuf parse (validate_phenopackets.py) apart from subject.extensions, which carry origin, place of care, migrant, onset, provenance and analysis labels (PAVS convention). Countries as GAZ terms (results/gaz_country_map.json, 178 names, OLS-checked) else text. Store and PAVS packets copied and augmented; PAVS copies fixed for two v2 violations in the source (bare-date metaData.created; modeOfInheritance on Disease, moved to an extension).

## All populations against all (2026-10-27, unimatrix job 10204)

OTHER split into SASIA, LATAM and AFR by `fine_regions.py` (patient origin, ancestry text, place of care, paper origins;
author region split the same way). Cases: EUR 27,285, EASIA 15,059, MENA 9,382, SASIA 1,724, LATAM 1,676, AFR 251.
`population_model.py allpairs`: the pairs decomposition for each of the 15 pairs, plus two direct contrasts per pair
(A-A pairs vs A-B pairs; B-B pairs vs A-B pairs), each with a group-permutation p (1000 perms), total and
patient-carried. 60 contrast tests; BH q in `results/v4/allpairs/contrasts.tsv`; full log `allpairs-10204.log`.

Contrasts with BH q < 0.05 (excess over the permutation null, z units; units = exchangeable groups):
- MENA-MENA vs MENA-EUR +0.17 (total, q=0.008; patient-carried +0.115, p=0.023, q=0.08). EUR-EUR vs MENA-EUR ~0.
- EASIA-EASIA vs EUR-EASIA +0.23 (total, q=0.008; patient-carried +0.06, ns).
- EASIA-EASIA vs MENA-EASIA +0.15 (total, q=0.008; patient-carried +0.085, ns). MENA-MENA vs MENA-EASIA -0.07, ns.
- EASIA-EASIA vs EASIA-SASIA +0.18 (total) and +0.28 (patient-carried, 188 groups).
- MENA-MENA vs MENA-SASIA +0.11 and SASIA-SASIA vs MENA-SASIA -0.19 (patient-carried only, 133 groups).
- Every A-A vs A-AFR contrast is large and positive (+0.6 to +0.7) and AFR-AFR vs A-AFR negative: AFR cases are
  less alike than any other pair type; 251 cases, 46-691 exchangeable groups; not interpretable as a population effect.
- Nothing for LATAM, EUR-SASIA, EUR-LATAM, EASIA-LATAM, SASIA-LATAM.
Reading: no general "same population, more alike" effect. EUR-EUR pairs are never more alike than EUR-X pairs.
MENA and EASIA cases are more alike among themselves; for EASIA this disappears when author region is held
fixed (care/reporting setting), for MENA it partly survives (vs EUR, vs SASIA). Effects are 0.1-0.3 z,
against +0.7 for two cases from the same reporting group.

## Gene recovery on all cases (2026-10-28, unimatrix job 10236)

`recovery_main.py` (vectorised; ranks identical to the previous version), 65,702 cases, library = store, pmcp,
pavs_lit, pubmed_abs; own publication and shared-author publications excluded; 2000 group permutations.
Statistic: within (gene, query corpus), mean log2 rank of population X minus EUR (negative = X recovered better).

| X vs EUR | strata | queries | diff | rank factor | total p | patient-carried excess | p |
|---|---|---|---|---|---|---|---|
| MENA | 1694 | 21375 | -0.104 | 0.93 | 0.083 | -0.052 | 0.016 |
| EASIA | 2371 | 34220 | -0.147 | 0.90 | 0.001 | +0.003 | 0.68 |
| SASIA | 835 | 12513 | -0.180 | 0.88 | 0.051 | -0.034 | 0.12 |
| LATAM | 737 | 12901 | -0.205 | 0.87 | 0.033 | -0.010 | 0.79 |
| AFR | 142 | 3562 | +0.042 | 1.03 | 0.75 | +0.005 | 0.98 |

Hits@1: Saudi clinical (pavs_clin) 1.8%, US clinical (pavs_us) 1.7%; MENA vs EUR literature: abstracts 13.1 vs
13.9, PMC-Patients 12.1 vs 10.8, Store 14.1 vs 18.9. No population is recovered worse than EUR. Non-EUR cases rank
slightly better in total; holding author region fixed removes this for all but MENA, where MENA patients rank
slightly better (rank factor ~0.96). Logs and per-case ranks in `results/v4/recovery/`.

## Test-based vs bedside findings for all populations, organ systems, age and sex (2026-10-28, unimatrix job 10420)

`population_model.py terms` on store + pmcp_map2 + pubmed_abs with fine population labels (results/v4/fine2, which add
`age_band` (infant <1 y, child 1-11, adolescent 12-17, adult) and `sex` from the phenopacket index via add_age_sex.py).
New: per-organ-system outcomes (findings, test-based, bedside in each child of HP:0000118), `--key-extra` (adds a column
to the stratum key) and `--where` (subset). Each population vs EUR; 2,000 group permutations. 50 runs; all outcomes in
results/v4/terms_all/all_outcomes.tsv.

Share of findings needing a test (X minus EUR, points; * p<0.05, ** p<0.005):

| X | total | patient-carried | age-adjusted | sex-adjusted | infant | child | adolescent | adult | male | female |
|---|---|---|---|---|---|---|---|---|---|---|
| MENA | -5.5** | -2.5 | -4.3** | -4.6** | -6.7** | -4.9** | -1.2 | -2.7 | -6.6** | -2.2 |
| SASIA | -7.0** | -9.0* | -6.6** | -8.1** | -7.3* | -8.1** | -11.9** | -3.0 | -9.8** | -5.8** |
| LATAM | -4.6** | -6.5* | -3.3* | -4.6** | +0.3 | -4.1 | -12.6* | -2.0 | -6.0** | -3.2 |
| AFR | -9.0** | -5.2 | -6.2* | -10.0** | +3.8 | -9.5* | +5.2 | -11.9* | -11.4** | -8.0 |
| EASIA | -0.6 | -3.7 | +0.3 | +0.5 | -5.5** | +0.6 | +2.9 | +1.3 | +0.4 | +0.7 |

- The lower test share in MENA, SASIA, LATAM and AFR comes from more bedside/history findings (+0.12 to +0.20
  log-count, all p=0.0005), not fewer test-based findings (-0.01 to -0.04, mostly ns). EASIA reports more findings
  of both kinds (+0.03 test, +0.06 bedside) and keeps the EUR share. Survives age and sex adjustment.
- Organ systems: the bedside excess is broad (13-16 of 21 systems for MENA and SASIA), largest in head/neck
  (+5.6 to +10 log-count x100), digestive, musculoskeletal, integument and nervous system, i.e. fuller examination
  and dysmorphology descriptions, not one organ. Test-based findings are higher only for blood (MENA, EASIA) and
  metabolism/respiratory/digestive (EASIA).
- Patient-carried: SASIA (-9.0) and LATAM (-6.5) keep the lower share; MENA does not (-2.5) but keeps fewer
  test-based findings (-0.053*), physiological tests (-0.059**) and imaging (-0.043*), concentrated in nervous
  system (-4.5*), head/neck (-3.9**) and eye (-3.2*) tests.
- Age: largest in infants and children; sex: larger in males for every population with a deficit.

## Health expenditure of the country of care (2026-10-28; healthspend.py)

World Bank indicators (API, results/v4/worldbank/): current health expenditure per capita (SH.XPD.CHEX.PC.CD),
GDP per capita PPP (NY.GDP.PCAP.PP.CD), physicians per 1,000 (SH.MED.PHYS.ZS); value in the publication year
(nearest available, 2000-2023); log10. Cases: PubMed abstracts + PMC-Patients of the six populations with one
place-of-care country matched to the World Bank (41,155 cases, 103 countries, 3,713 gene x corpus strata; Taiwan,
422 cases, has no World Bank data). Model: outcome ~ covariate with gene x corpus fixed effects; covariate permuted
between countries (2,000), 95% CI by bootstrap over countries.

| covariate (per 10-fold) | test share (points) | log(1+ test-based) | log(1+ bedside) |
|---|---|---|---|
| health expenditure per capita | +5.0 [4.3, 6.1], p=0.0005 | +0.019, p=0.059 | -0.112 [-0.145, -0.093], p=0.0005 |
| GDP per capita PPP | +9.6 [7.0, 11.9], p=0.0005 | +0.038, p=0.020 | -0.214, p=0.0005 |
| physicians per 1,000 | +9.2 [5.2, 12.7], p=0.0005 | +0.029, p=0.052 | -0.220, p=0.0005 |

Population differences in test share vs EUR without -> with health expenditure in the model (points):
MENA -5.4 -> -1.2; SASIA -7.8 -> +0.5; LATAM -4.6 -> -1.3; AFR -7.4 -> -3.7; EASIA -1.6 -> +1.0.
Bedside excess: MENA +0.140 -> +0.063; SASIA +0.178 -> +0.028; LATAM +0.123 -> +0.063; AFR +0.171 -> +0.104.
Health expenditure absorbs most of the population differences; physician density absorbs little. Higher-spending
countries report fewer bedside/history findings rather than more test-based findings. Outliers: Saudi Arabia and
other Gulf states spend like southern Europe but sit ~5 points below (fig_healthspend).

## Network figures (2026-10-28; network_figs.py; results/v4/figures/)

- fig_shift_network: 140 HPO terms 2-3 levels below the root, edges co-occurrence phi > 0.08 over 62,220 patients,
  colour = term z per population vs EUR, shape = finding class, black ring = q < 0.05; shared layout, one panel per
  population. Test-based cluster (circulating metabolites, liver/spleen, blood) shifts up in EASIA; bedside cluster
  (skin, face/skull, limbs, digestive, neurodevelopment) shifts up in MENA and SASIA.
- fig_dag_<POP>: terms with BH q < 0.05 and their is_a links after transitive reduction, radial from the root
  (MENA 150 terms, SASIA 107, EASIA 231, LATAM and AFR none).
- fig_healthspend: country means of within-stratum residual test share vs health expenditure per capita.

## Within-population vs outside similarity (2026-10-28, unimatrix job 10491; `population_model.py allpairs --vs-rest`)

Each population P against all other populations pooled (REST); contrast P-P vs P-REST pairs, same decomposition and
permutations as the all-pairs analysis. `results/v4/within/contrasts.tsv`; BH q over the 12 tests in
`results/v4/qvalues.tsv` (family within_vs_outside).

| P | total excess | q | patient-carried excess | q |
|---|---|---|---|---|
| EASIA | +0.216 | 0.012 | +0.071 | 0.35 |
| MENA | +0.058 | 0.19 | +0.078 | 0.19 |
| LATAM | +0.004 | 0.96 | +0.168 | 0.13 |
| SASIA | +0.040 | 0.72 | -0.208 | 0.19 |
| AFR | -0.500 | 0.46 | -0.747 | 0.066 |
| EUR | -0.038 | 0.46 | -0.024 | 0.60 |

Only East Asian cases are more alike among themselves than to cases from elsewhere, and only in the total version
(provenance included). Against the pooled rest, the MENA within-population similarity seen against EUR alone
(+0.17, q=0.008) shrinks to +0.06 to +0.08 and is not significant.

## Ancestry vs place of care, all populations (2026-10-28, unimatrix jobs 10488_1, 10492; `results/v4/setting_all/`)

Share of findings needing a test, difference in percentage points (r1 minus r2), total strata, 2,000 permutations;
BH q over the 15 contrasts (`results/v4/setting_all/test_share_q.tsv`).

| Contrast | diff | q | cases (r1) |
|---|---|---|---|
| MENA ancestry, MENA care vs EUR | -8.0 | 0.004 | 11,944 (1,673) |
| MENA ancestry, EUR care vs EUR | -0.3 | 0.81 | 5,017 (581) |
| MENA ancestry: MENA care vs EUR care | -7.2 | 0.019 | 824 (507) |
| SASIA ancestry, SASIA care vs EUR | -9.4 | 0.004 | 5,367 (512) |
| SASIA ancestry, EUR care vs EUR | -3.7 | 0.36 | 1,345 (83) |
| SASIA ancestry: SASIA care vs EUR care | -14.3 | 0.16 | 81 (51) |
| AFR ancestry, AFR care vs EUR | -12.5 | 0.015 | 1,528 (62) |
| AFR ancestry, EUR care vs EUR | -2.5 | 0.42 | 3,235 (202) |
| AFR ancestry: AFR care vs EUR care | -5.3 | 0.66 | 53 (20) |
| LATAM ancestry, LATAM care vs EUR | -2.0 | 0.36 | 7,164 (529) |
| LATAM ancestry, EUR care vs EUR | -6.1 | 0.029 | 2,299 (205) |
| LATAM ancestry: LATAM care vs EUR care | +2.3 | 0.81 | 110 (72) |
| EASIA ancestry, EASIA care vs EUR | -1.4 | 0.029 | 22,906 (5,810) |
| EASIA ancestry, EUR care vs EUR | -4.8 | 0.27 | 1,950 (91) |
| EASIA ancestry: EASIA care vs EUR care | +4.7 | 0.34 | 680 (600) |

The lower test share follows the place of care, not ancestry: patients of MENA, South Asian or African ancestry
cared for in Europe or North America are reported like European patients, while the same ancestries cared for at
home are reported with fewer test-based findings. The direct within-ancestry comparison is significant for MENA
(-7.2 points, q=0.019) and points the same way for SASIA (-14.3) and AFR (-5.3) with few cases. LATAM is the
exception: the gap appears for LATAM-ancestry patients cared for in Europe/North America.
