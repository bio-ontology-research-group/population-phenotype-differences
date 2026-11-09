# Blind judge of mapping disagreements
session 818691a1  timestamp 2026-10-26T07:21:48.830Z  model sonnet  subagent_type general-purpose

You are judging how clinical findings from case reports were mapped to Human Phenotype Ontology (HPO) terms. Two systems disagree; you see their answers as A and B in random order and must judge each independently and strictly.

Input: /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/map2/disagreements_blind.json (254 items: id, phrase, A {hpo, name}, B {hpo, name}; hpo may be null meaning "no term assigned").
HPO for lookups: /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/hpo/hp.obo (grep for definitions, synonyms and is_a parents; check hierarchy before calling something broader or more specific).

For each item judge A and B separately as one of:
- "exact": the term means what the phrase says
- "broader_ok": correct but more general than the phrase (an ancestor of the ideal term)
- "too_specific": asserts more than the phrase (a subtype, a particular body part, a cause, or a direction/severity not stated)
- "wrong": a different finding
- "none_ok": null, and indeed no HPO term fits (e.g. not an abnormality)
- "none_missed": null, but a fitting HPO term exists
Then say which is better: "A", "B", or "tie".

Write /home/leechuck/tmp/claude-1000/-home-leechuck-Documents-papers-pavs-saudi-specific-phenotypes/818691a1-fa28-4b8f-9820-3dbbc0854147/scratchpad/map2/disagreements_judged.json as a list of {"id", "A", "B", "better"} covering all 254 items. Do not modify other files. Return counts of each label for A and for B, and of "better".
