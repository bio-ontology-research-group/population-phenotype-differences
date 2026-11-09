# Language-model judge prompts (extraction and mapping audits)

Verbatim prompts given to the judge agents described in `docs/extraction_audit.md`
and in the Supplementary Methods ("Extraction performance"). Each file starts with
the agent description, the session id, the UTC timestamp at which the agent was
started, and the model alias requested.

| File | Audit | Cases | Judged file(s) |
|---|---|---|---|
| `pmcp_batch1_2026-10-26.md` .. `pmcp_batch3_2026-10-26.md` | PMC-Patients extraction audit, three region-blind batches of 20 | 60 | `results/v4/audit_pmcp/batch{1,2,3}_judged.json` |
| `mapping_blind_judge_2026-10-26.md` | blind pairwise judgement of the 254 phrases where the lexical and RAG mappers disagree | 254 phrases | see RESULTS.md, "Extraction and normalisation audit" |
| `abstracts_batches1-2_2026-10-28.md`, `abstracts_batches3-4_2026-10-28.md` | PubMed-abstract extraction audit, four region-blind batches of 30 | 120 | `results/v4/audit_abs/batch{1,2,3,4}_judged.json` |

Judge model and settings
- Model: Claude Sonnet 5.5 (`claude-sonnet-5-5`, Anthropic), run as Claude Code
  subagents (Agent tool, model alias `sonnet`, general-purpose agent type) on
  2026-10-26 (PMC-Patients, mapping judge) and 2026-10-28 (abstracts). The
  session logs record the resolved model identifier for every agent turn.
- Sampling: the Agent tool's default settings; no temperature or system prompt
  was set beyond the prompt text in these files. The agents had read access to
  the batch files and the HPO `hp.obo` named in the prompt and wrote the
  `*_judged.json` files directly.
- Sample selection and batch assignment: `results/v4/audit_abs/regions.json`
  (abstracts: 20 per population, one patient per publication, seed 2026,
  shuffled into four batches of 30); PMC-Patients: 20 cases each from EASIA, EUR
  and MENA shuffled into three batches of 20. Region labels, author lists and
  affiliations were withheld from the judges.
- Adjudication: none; each phrase received one judgement. Scoring is described
  in `docs/extraction_audit.md`.

The batch input files held PMC-Patients summaries (CC BY-NC-SA) and abstract
text and are not redistributed; the paths in the prompts are those of the
session scratchpad in which the audit ran.
