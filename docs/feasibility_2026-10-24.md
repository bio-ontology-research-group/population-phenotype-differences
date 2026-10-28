# Feasibility analysis (2026-10-24)

Derived data used to assess whether the population-specific-reference hypothesis
is testable with PAVS and the GA4GH Phenopacket Store (v0.1.26).

## Files

- `phenopacket_store_0.1.26_pmid_country.tsv`: every source PMID of the
  Phenopacket Store v0.1.26 (1,584 PMIDs) with year, journal, country of the
  first author's affiliation (regex over PubMed AffiliationInfo), all author
  countries with counts, whether the abstract mentions consanguinity, and
  geographic MeSH headings. Country is a proxy for patient origin; 118 PMIDs
  have no affiliation country in PubMed.

## Facts established from the local PAVS repo (`~/Public/software/pavs`)

- The PAVS "DDD cohort" (prefix D, 1,856 records) is the DDG2P gene panel
  (`submodules/knowledge-graph/data/phenotypes/ddd-diagnoses.tsv`: one row per
  gene/disease with a curated HPO list). The first 50 D phenopackets have HPO
  sets identical to their DDG2P row; 1,856 records cover 1,339 genes. These
  are not patients.
- Source Q (Ziats et al. 2020, Pediatr Res, PMID 31618753) is a US cohort
  (all 20 author affiliations USA); the PAVS phenopackets tag it
  `geographicOrigin = Saudi Arabia`.
- Per source: A 251 cases/238 genes, B 223/191, F 1,019 cases (595 solved)/441
  genes, P 2,217 (1,616 solved)/1,097 genes, M 1,422 (1,308 solved)/478 genes,
  Q 522 (142 solved)/114 genes.
- Cross-source gene coverage within Saudi clinical sources (A,B,F,P): 325 genes
  appear in >=2 sources with >=1 case, 15 genes with >=3 cases in each of >=2
  sources, 0 genes with >=5 cases in each of >=2 sources. Both institutions
  (KFSHRC = F,P; KAMC = A,B) with >=3 cases each: 6 genes.
- Saudi clinical vs Saudi literature (M): 41 genes with >=3 cases on each side,
  11 with >=5, 0 with >=10.
- M is drawn from 104 PMIDs; the five largest contribute 129, 113, 97, 90, 87
  cases; median 4 cases per PMID; 33 PMIDs with >=10 cases.
- M vs Phenopacket Store publication overlap: 5 shared PMIDs; 129 M cases vs
  14 Store cases from those PMIDs (same-paper, two-curator comparison is
  possible but tiny).
- Phenopacket Store by first-author region: Eur/NA/Oceania 7,039 cases,
  East Asia 1,223, MENA+ (Middle East, North Africa, Iran, Turkey, Pakistan)
  566, other 208, unknown 552. Mean HPO terms/case: 8.9, 7.2, 8.9. Saudi
  affiliation: 130 cases, 31 genes.
- Genes with >=3 cases in both MENA+ and Eur/NA within the Store: 23 (>=5: 10,
  >=10: 5). East Asia vs Eur/NA: 41 (>=5: 34, >=10: 18).
- Genes with >=3 PAVS Saudi cases and >=3 Store Eur/NA cases: 34 (>=5: 16).

## Sample sizes for the candidate designs (first-author-country proxy for Store patient origin)

- Design "Saudi-derived reference, Store test split by region": genes with >=3 M
  cases and >=3 Store cases in both MENA+ and Eur/NA: 2 (RPGRIP1, WWOX); with
  all Saudi PAVS sources >=3: 7 genes (216 MENA+ / 397 Eur/NA test cases, but
  dominated by LMNA, FBXL4, SPTAN1 where one MENA paper supplies most cases).
- Design "Store MENA+ vs Store Eur/NA reference (same curator), PAVS Saudi
  clinical test": 12 genes with >=3 Store cases on each side, 48 PAVS
  Saudi-clinical test cases (7 genes / 30 cases at >=5).
- Within-Store profile comparison (H1), genes with >=5 cases on each side: LMNA
  (115 MENA+ cases from 1 paper vs 135 from 7), FBXL4 (39/3 vs 54/10), SPTAN1
  (27/1 vs 55/5), RPGRIP1 (25/4 vs 36/7), SNUPN (18/1 vs 5/1), AIRE (13/1 vs
  35/3), ATP13A2 (9/4 vs 26/9), TBX5 (7/4 vs 121/12), SPINK5 (7/5 vs 63/14),
  NHEJ1 (5/1 vs 15/4). Where the MENA+ side comes from one paper, study and
  population are not separable for that gene.
- Trujillano et al. 2017 (Centogene, EJHG, PMID 27848944) Supplementary Table
  S3, fetched from Europe PMC: 320 solved probands, Geographic Region ME 244,
  Eu 33, SA 21, Am 8, Oceania 1, blank 13; 248 genes, 13 genes with >=3 cases;
  phenotypes given as comma-separated HPO labels (not IDs), one lab, one
  phenotyping pipeline.
