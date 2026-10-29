# Source, genetics and frequency extension

Active user objective: complete all five analyses below, document their results,
and then reconsider and revise the manuscript and supplement. Existing v4 outputs
remain immutable; new outputs go to `results/v5/`. Raw text and audit review packets
stay outside version control. Human reviewers: Paul Schofield and Marwa Abdelhakim
(confirmed by Robert on 2026-10-29). No annotations will be labelled human-validated
until their actual assessments have been returned and adjudicated.

## Completion evidence

1. Matched publication/source manifest, abstract and full-text extraction with the
   same checkpoint/prompt/mapping, eligibility and alignment checks, blinded human
   annotations, adjudication, and paired source-bias estimates with uncertainty.
2. Separate source-specific regional phenotype comparisons and paired gene-ranking
   changes, denominators/exclusions, publication-cluster uncertainty and multiplicity
   control; evaluated from the original case-level data.
3. Gene/disorder-specific regional contrasts, source/care-adjusted estimates where
   identifiable, explicit assessment of genetic metadata availability, subtype,
   inheritance/variant sensitivities where supported, and multiple-testing control.
4. HPOA frequency-format and reference audit; publication-to-cohort provenance
   mapping (stated ancestry separate from care and author location), overlap checks,
   and within-disorder frequency analyses where data support them.
5. Source-verifiable affected/assessed denominators for selected disorders, unknown
   separated from absent, independent validation cohorts/publications, and a direct
   comparison of rankings using original versus improved frequencies.
6. Complete result report including null/failed analyses and limitations; revised
   paper and supplement with verified numbers, references, figures and compilation.

## Design decisions before new analyses

- Retain v4 source boundaries (Store, PMC-Patients, abstracts). Shared publications
  do not constitute independent replication. Never count repeated extractions as
  additional patients or studies.
- Define test share only for cases with at least one classified finding in new
  primary analyses; retain v4 zero-assignment results as an explicitly named
  sensitivity. Describe missingness for every source and region.
- Source-specific paired rank change uses the same full/bedside ranking for each
  case and publication-cluster bootstrap intervals. Add reporting-group bootstrap
  sensitivity. These are information-restriction benchmarks, not diagnostic yield.
- Preserve ancestry, place of care and reporting affiliation as distinct variables.
  A population assigned from care country cannot establish an ancestry effect after
  conditioning on that country. Report lack of overlapping groups as non-identifiability.
- Gene-specific analyses are exploratory; set support thresholds before testing,
  report all eligible comparisons and adjust across the complete tested family.
  Matching a gene does not equate disease subtype or allelic mechanism.
- HPOA frequency categories are intervals, not invented exact proportions; missing
  denominators remain missing. An unmentioned phenotype remains unassessed/unknown.
- Freeze training/reference versus evaluation publications before frequency-model
  comparison. Check shared patients and families where identifiable, and exclude
  development examples from independent evaluation.

## State, 2026-10-29

- Local analysis repository started clean at 531f4bee7. Manuscript source is at
  f3d9861; an author-facing review memo is untracked and is being retained.
- No project process from the interrupted goal was found running locally or in
  Unimatrix's job queue.
- Recovered original HPO inputs from unimatrix01:/data/hohndor/pavs-popmodel/hpo/.
- Human review is outstanding; preparing reviewer materials is part of this work.

## Progress checkpoint

- Source-specific analysis completed: 18 paired source–region benchmarks with
  publication and reporting-group resampling, plus 180 regional tests across
  phenotype/ranking and total/author-adjusted families. Results are documented
  in `results/v5/RESULTS.md`; new estimates are not yet in the manuscript.
- HPOA format/provenance audit and case genetic-metadata coverage completed;
  cohort-level source verification and frequency modelling remain outstanding.
- Exploratory gene analysis completed under prespecified minimum support
  (10 cases and 3 publications in each group after matching). No test passes
  correction across the complete tested families. Care-conditioned gene effects
  lack supported overlap; variant-level follow-up remains outstanding.
- Paired-source sample: 110 publications selected before access checks, 59 paired
  inputs recovered. Ten PMC full-text fetch failures were all non-OA metadata
  returning HTTP 500; distinguish this from a successful source retrieval.
- Human reference packets prepared for Paul and Marwa, 20 publications and 40
  documents per reviewer, model predictions and pairing keys withheld.
- Original direct-HTTP pilot completed with transient gateway errors; an SSH
  tunnel through Unimatrix gives a successful small request in 1.42 seconds.
  Full paired extraction is now running through that tunnel, resumably; inspect
  the live process and log before any restart. Raw inputs, model outputs and
  review packets are under ignored `data/v5_private/`.
