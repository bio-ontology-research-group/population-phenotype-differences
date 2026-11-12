#!/usr/bin/env python3
"""Deterministic case selection for tool_restriction: per region, cases sorted by sha1(id) hex digest.

  --per-region N   take the first N per region (pilot: 10)
  --cap C          stratified sample: all cases of regions with <= C cases, the first C (by sha1) of larger ones
  (neither)        all cases
"""
import argparse, collections, hashlib, json


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--profiles", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--per-region", type=int); ap.add_argument("--cap", type=int)
    a = ap.parse_args()
    by = collections.defaultdict(list)
    for l in open(a.profiles):
        c = json.loads(l); by[c["region"]].append(c["id"])
    out = []
    for reg, ids in sorted(by.items()):
        ids = sorted(ids, key=lambda i: hashlib.sha1(i.encode()).hexdigest())
        k = a.per_region if a.per_region else (a.cap if a.cap and len(ids) > a.cap else len(ids))
        out += ids[:k]
        print(reg, len(ids), "->", min(k, len(ids)))
    open(a.out, "w").write("\n".join(out) + "\n")
    print(len(out), "selected")


if __name__ == "__main__":
    main()
