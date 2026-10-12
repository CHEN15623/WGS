#!/usr/bin/env python3
"""
cog_profile.py -- assign COG functional categories to the reference and query
proteomes and compare the two profiles.

The NCBI COG2020 release is used:
    cog-20.fa.gz     member protein sequences   (searched with DIAMOND)
    cog-20.cog.csv   protein -> COG assignment
    cog-20.def.tab   COG -> functional category letter + name
and the COG functional categories (J A K L B D Y V T M N Z W U O C G E F H I
P Q R S) are reported per genome, which is what makes the "metabolism" question
answerable without a proprietary pathway database.

Outputs (in --outdir)
  cog_assignments_reference.tsv / cog_assignments_query.tsv
  cog_category_comparison.tsv     per-category counts, reference vs query
  cog_id_comparison.tsv           per-COG-id counts, reference vs query
  cog_profile_summary.json
"""
import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

COG_CATEGORIES = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Post-translational modification, protein turnover, chaperones",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}

METABOLISM_LETTERS = set("CEFGHIPQ")


def canon(acc):
    """Canonical accession with the version suffix removed.

    The COG csv writes versions with a dot (NP_811616.1) while the COG FASTA
    mangles them to an underscore (NP_811616_1), so both separators are handled.
    """
    acc = acc.strip().upper()
    m = re.match(r"^(.*?)[._](\d+)$", acc)
    return m.group(1) if m else acc


def run(cmd, **kw):
    return subprocess.run(["bash", "-lc", cmd], capture_output=True,
                          text=True, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cog-dir", required=True,
                    help="directory holding cog-20.fa.gz / .cog.csv / .def.tab")
    ap.add_argument("--ref-faa", required=True)
    ap.add_argument("--qry-faa", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--identity", type=float, default=30.0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    fa_gz = os.path.join(a.cog_dir, "cog-20.fa.gz")
    csv_f = os.path.join(a.cog_dir, "cog-20.cog.csv")
    def_f = os.path.join(a.cog_dir, "cog-20.def.tab")
    for f in (fa_gz, csv_f, def_f):
        if not os.path.exists(f):
            print(f"[cog] missing {f}", file=sys.stderr)
            return

    if not shutil.which("diamond"):
        print("[cog] DIAMOND not available -- skipping COG profiling",
              file=sys.stderr)
        return

    # ---- COG -> category ---------------------------------------------------
    cog2cat, cog2name = {}, {}
    for line in open(def_f, errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            cog2cat[f[0]] = f[1]
            cog2name[f[0]] = f[2]

    # ---- protein accession -> COG -----------------------------------------
    acc2cog = {}
    with open(csv_f, errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split(",")
            if len(f) >= 7 and f[6].startswith("COG") and f[2]:
                acc2cog[canon(f[2])] = f[6]
    print(f"[cog] protein->COG mappings: {len(acc2cog):,}")

    # ---- DIAMOND database --------------------------------------------------
    # the COG database is a shared reference resource: build it in the COG
    # directory (work/cog), never inside the result tree
    fa = os.path.join(a.cog_dir, "cog-20.fa")
    db = os.path.join(a.cog_dir, "cog-20.dmnd")
    if not os.path.exists(fa):
        with gzip.open(fa_gz, "rb") as fi, open(fa, "wb") as fo:
            shutil.copyfileobj(fi, fo, length=1 << 22)
    if not os.path.exists(db):
        r = run(f"diamond makedb --in '{fa}' -d '{db}' --quiet")
        if r.returncode != 0:
            print(f"[cog] diamond makedb failed: {r.stderr[:300]}",
                  file=sys.stderr)
            return

    # ---- search both proteomes --------------------------------------------
    results = {}
    for label, faa in (("reference", a.ref_faa), ("query", a.qry_faa)):
        out = os.path.join(a.outdir, f"diamond_{label}.tsv")
        if not os.path.exists(out):
            # NOTE: option names differ between DIAMOND releases; the identity
            # threshold is applied in Python below instead of on the command line.
            r = run(f"diamond blastp -q '{faa}' -d '{db}' -o '{out}' "
                    f"--evalue {a.evalue} "
                    f"--max-target-seqs 1 --threads {a.threads} --quiet "
                    f"-f 6 qseqid sseqid pident length evalue bitscore")
            if r.returncode != 0:
                print(f"[cog] diamond blastp ({label}) failed: {r.stderr[:300]}",
                      file=sys.stderr)
                continue
        results[label] = out

    # ---- summarise ---------------------------------------------------------
    per_genome = {}
    for label, path in results.items():
        rows = []
        cat_counts = Counter()
        cog_counts = Counter()
        assigned = 0
        total = 0
        for line in open(path, errors="replace"):
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            total += 1
            qid, sid, pid, aln, ev, bits = f[:6]
            try:
                if float(pid) < a.identity:
                    continue
            except ValueError:
                pass
            cog = acc2cog.get(canon(sid), "")
            cats = cog2cat.get(cog, "")
            if cog:
                assigned += 1
                cog_counts[cog] += 1
                for L in (cats or "S"):
                    if L in COG_CATEGORIES:
                        cat_counts[L] += 1
            rows.append((qid, sid, pid, aln, ev, cog, cats))
        with open(os.path.join(a.outdir, f"cog_assignments_{label}.tsv"),
                  "w") as fh:
            fh.write("query_protein\tcog_member_protein\tpercent_identity\t"
                     "alignment_length\tevalue\tCOG_id\tCOG_categories\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        per_genome[label] = dict(proteins_with_hit=total, assigned=assigned,
                                 category_counts=dict(cat_counts),
                                 cog_counts=dict(cog_counts))
        print(f"[cog] {label}: {assigned}/{total} proteins assigned a COG")

    if len(per_genome) < 2:
        with open(os.path.join(a.outdir, "cog_profile_summary.json"), "w") as fh:
            json.dump(per_genome, fh, indent=2)
        return

    rc = per_genome["reference"]["category_counts"]
    qc = per_genome["query"]["category_counts"]
    with open(os.path.join(a.outdir, "cog_category_comparison.tsv"), "w") as fh:
        fh.write("COG_letter\tCOG_category\tmetabolism\t"
                 "reference_count\tquery_count\tdifference\n")
        for L, name in COG_CATEGORIES.items():
            r_, q_ = rc.get(L, 0), qc.get(L, 0)
            fh.write(f"{L}\t{name}\t{'yes' if L in METABOLISM_LETTERS else 'no'}\t"
                     f"{r_}\t{q_}\t{q_-r_:+d}\n")

    rcid = per_genome["reference"]["cog_counts"]
    qcid = per_genome["query"]["cog_counts"]
    with open(os.path.join(a.outdir, "cog_id_comparison.tsv"), "w") as fh:
        fh.write("COG_id\tCOG_category\tCOG_name\treference_count\t"
                 "query_count\tstatus\n")
        for cid in sorted(set(rcid) | set(qcid)):
            r_, q_ = rcid.get(cid, 0), qcid.get(cid, 0)
            st = "shared" if r_ and q_ else ("reference_only" if r_ else "query_only")
            fh.write(f"{cid}\t{cog2cat.get(cid,'')}\t{cog2name.get(cid,'')}\t"
                     f"{r_}\t{q_}\t{st}\n")

    summary = dict(
        reference=per_genome["reference"],
        query=per_genome["query"],
        metabolism_category_totals={
            "reference": sum(rc.get(L, 0) for L in METABOLISM_LETTERS),
            "query": sum(qc.get(L, 0) for L in METABOLISM_LETTERS)},
        cog_ids_reference_only=sorted(set(rcid) - set(qcid)),
        cog_ids_query_only=sorted(set(qcid) - set(rcid)),
    )
    summary["reference"].pop("cog_counts", None)
    summary["query"].pop("cog_counts", None)
    with open(os.path.join(a.outdir, "cog_profile_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print("[cog] metabolism-category assignments: "
          f"reference={summary['metabolism_category_totals']['reference']} "
          f"query={summary['metabolism_category_totals']['query']}")


if __name__ == "__main__":
    main()
