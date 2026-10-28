#!/usr/bin/env python3
"""Validate phenopacket JSON files against the GA4GH Phenopacket v2 protobuf schema.

Strict parse after removing `subject.extensions` (PAVS-style additions with no v2 slot);
reports any other schema violation. Needs the `phenopackets` package.
Usage: validate_phenopackets.py DIR
"""
import collections, glob, json, os, sys
from google.protobuf import json_format
from phenopackets import Phenopacket

bad = collections.Counter(); n = 0; ex = []
for f in glob.glob(os.path.join(sys.argv[1], "*", "*.json")):
    d = json.load(open(f)); d["subject"].pop("extensions", None); n += 1
    try:
        json_format.Parse(json.dumps(d), Phenopacket())
    except json_format.ParseError as e:
        msg = str(e).split(":")[-1].strip()[:120]; bad[msg] += 1
        if len(ex) < 5: ex.append((f, str(e)[:300]))
print(f"{n} files, {sum(bad.values())} invalid")
for k, v in bad.most_common(10): print(v, k)
for e in ex: print(e)
