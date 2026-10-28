#!/usr/bin/env python3
"""Remove verbatim source-text fields from language-model output files before distribution.

The origin and setting extractions (store_origin_extract.py, pmcp_extract.py origin, qwen_tasks.py setting) ask the
model to quote the sentence that supports each label. Those quotes are article text and are not redistributed; this
script drops every *evidence* field, the free-text notes and the raw model response, in place, and keeps the labels.

  tools/strip_evidence.py FILE.jsonl [...]
"""
import json, sys

DROP_EXACT = {"notes", "raw"}


def strip(obj):
    if isinstance(obj, dict):
        return {k: strip(v) for k, v in obj.items() if not (k in DROP_EXACT or "evidence" in k)}
    if isinstance(obj, list):
        return [strip(x) for x in obj]
    return obj


for path in sys.argv[1:]:
    rows = [strip(json.loads(l)) for l in open(path) if l.strip()]
    with open(path, "w") as f:
        for r in rows: f.write(json.dumps(r, ensure_ascii=False) + "\n")
    print(path, len(rows))
