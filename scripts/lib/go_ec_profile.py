#!/usr/bin/env python3
"""
go_ec_profile.py -- compare the GO and EC functional profiles of the reference
and the (projected) query annotation.

EC numbers give the enzymatic/metabolic repertoire; GO terms give the
process/function/component breakdown.  Both are taken straight from the
annotation, so no additional database is needed.

Outputs (in --outdir)
  go_term_comparison.tsv
  go_namespace_summary.tsv
  ec_number_comparison.tsv
  ec_class_summary.tsv
  go_ec_summary.json
"""
import argparse
import json
import os
import re
import sys
from collections import Counter, defaultdict

EC_NAMES = {"1": "Oxidoreductases", "2": "Transferases", "3": "Hydrolases",
            "4": "Lyases", "5": "Isomerases", "6": "Ligases",
            "7": "Translocases"}


def parse(path, want_go=True, want_ec=True):
    go, ec = Counter(), Counter()
    n = 0
    if not path or not os.path.exists(path):
        return go, ec, n
    for line in open(path, errors="replace"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        n += 1
        attrs = {}
        for kv in f[8].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                attrs.setdefault(k, v)
        if want_go:
            for t in re.split(r"[,\s]+", attrs.get("Ontology_term", "")):
                t = t.strip()
                if t.startswith("GO:"):
                    go[t] += 1
        if want_ec:
            for t in re.split(r"[,\s]+", attrs.get("EC_number", "")):
                t = t.strip()
                if re.match(r"^\d+\.", t):
                    ec[t] += 1
    return go, ec, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--qry-gff", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rgo, rec, rn = parse(a.ref_gff)
    qgo, qec, qn = parse(a.qry_gff)
    print(f"[go_ec] reference CDS={rn} GO={len(rgo)} EC={len(rec)}")
    print(f"[go_ec] query     CDS={qn} GO={len(qgo)} EC={len(qec)}")

    with open(os.path.join(a.outdir, "go_term_comparison.tsv"), "w") as fh:
        fh.write("GO_term\treference_count\tquery_count\tstatus\n")
        for t in sorted(set(rgo) | set(qgo)):
            r_, q_ = rgo.get(t, 0), qgo.get(t, 0)
            st = "shared" if r_ and q_ else ("reference_only" if r_ else "query_only")
            fh.write(f"{t}\t{r_}\t{q_}\t{st}\n")

    with open(os.path.join(a.outdir, "ec_number_comparison.tsv"), "w") as fh:
        fh.write("EC_number\tEC_class\treference_count\tquery_count\tstatus\n")
        for t in sorted(set(rec) | set(qec)):
            r_, q_ = rec.get(t, 0), qec.get(t, 0)
            st = "shared" if r_ and q_ else ("reference_only" if r_ else "query_only")
            fh.write(f"{t}\t{t.split('.')[0]}\t{r_}\t{q_}\t{st}\n")

    with open(os.path.join(a.outdir, "ec_class_summary.tsv"), "w") as fh:
        fh.write("EC_class\tdescription\treference_count\tquery_count\t"
                 "difference\n")
        for cls, name in EC_NAMES.items():
            r_ = sum(v for k, v in rec.items() if k.startswith(cls + "."))
            q_ = sum(v for k, v in qec.items() if k.startswith(cls + "."))
            fh.write(f"{cls}\t{name}\t{r_}\t{q_}\t{q_-r_:+d}\n")

    # GO terms carry their namespace in the reference GFF only implicitly, so a
    # lightweight breakdown by term-count is reported instead of inventing one.
    with open(os.path.join(a.outdir, "go_namespace_summary.tsv"), "w") as fh:
        fh.write("metric\treference\tquery\tdifference\n")
        fh.write(f"distinct_GO_terms\t{len(rgo)}\t{len(qgo)}\t{len(qgo)-len(rgo):+d}\n")
        fh.write(f"total_GO_annotations\t{sum(rgo.values())}\t"
                 f"{sum(qgo.values())}\t{sum(qgo.values())-sum(rgo.values()):+d}\n")
        fh.write(f"distinct_EC_numbers\t{len(rec)}\t{len(qec)}\t"
                 f"{len(qec)-len(rec):+d}\n")
        fh.write(f"total_EC_annotations\t{sum(rec.values())}\t"
                 f"{sum(qec.values())}\t{sum(qec.values())-sum(rec.values()):+d}\n")
        fh.write(f"CDS_annotated\t{rn}\t{qn}\t{qn-rn:+d}\n")

    summary = dict(
        reference_CDS=rn, query_CDS=qn,
        reference_go_terms=len(rgo), query_go_terms=len(qgo),
        reference_ec_numbers=len(rec), query_ec_numbers=len(qec),
        ec_class_totals={
            "reference": {c: sum(v for k, v in rec.items() if k.startswith(c + "."))
                          for c in EC_NAMES},
            "query": {c: sum(v for k, v in qec.items() if k.startswith(c + "."))
                      for c in EC_NAMES}},
        go_reference_only=sorted(set(rgo) - set(qgo))[:200],
        go_query_only=sorted(set(qgo) - set(rgo))[:200],
        ec_reference_only=sorted(set(rec) - set(qec))[:200],
        ec_query_only=sorted(set(qec) - set(rec))[:200],
    )
    with open(os.path.join(a.outdir, "go_ec_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[go_ec] GO reference-only={len(summary['go_reference_only'])}, "
          f"query-only={len(summary['go_query_only'])}; "
          f"EC reference-only={len(summary['ec_reference_only'])}, "
          f"query-only={len(summary['ec_query_only'])}")


if __name__ == "__main__":
    main()
