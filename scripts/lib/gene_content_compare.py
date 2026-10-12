#!/usr/bin/env python3
"""
gene_content_compare.py -- which genes are shared, absent or divergent between
                           the reference and the query assembly.

Method
  The asm5 BAM holds the query contigs aligned to the reference.  Walking every
  CIGAR once yields the set of reference positions that are covered and the set
  of query positions that are covered.  A gene whose own positions are largely
  NOT covered is either absent from the other genome or too divergent to align.

  Optionally the query proteins are searched against the reference proteins with
  DIAMOND/BLAST to attach a sequence-identity figure to every shared gene.

Outputs (in --outdir)
  reference_gene_presence.tsv   coverage of each reference CDS in the query
  query_gene_presence.tsv       coverage of each query CDS in the reference
  absent_from_query.tsv         reference genes missing in the query
  query_specific_genes.tsv      query genes missing in the reference
  gene_content_summary.json
  shared_gene_identity.tsv      (when diamond/blastp is available)
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict


def parse_gff_cds(path, key_attr="locus_tag"):
    genes = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        d = {}
        for kv in f[8].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                d.setdefault(k, v)
        locus = d.get(key_attr) or d.get("ID") or d.get("Parent") or ""
        genes.append(dict(contig=f[0], start=int(f[3]), end=int(f[4]),
                          strand=f[6], locus=locus,
                          gene=d.get("gene", d.get("Name", "")),
                          product=d.get("product", "")))
    return genes


def walk_bam(bam_path):
    """Return the set of covered REFERENCE positions, {seqname: set(0-based)}.

    Reference coordinates come straight from the CIGAR and are always correct.
    Query coordinates are deliberately NOT derived here: the asm5 alignment uses
    supplementary records with hard clips, so SEQ indices are not contig
    coordinates.  Query coverage is computed from the NUCmer blocks instead
    (see `block_coverage`).
    """
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    cov_ref = defaultdict(set)
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.cigartuples is None:
            continue
        rname = bam.get_reference_name(read.reference_id)
        rp = read.reference_start
        for op, ln in read.cigartuples:
            if op in (0, 7, 8):                      # M / = / X
                for k in range(ln):
                    cov_ref[rname].add(rp + k)
                rp += ln
            elif op in (2, 3):                       # D / N
                rp += ln
    bam.close()
    return cov_ref, {}


def block_coverage(coords_path):
    """Query positions covered by the 1-to-1 alignment, from NUCmer blocks."""
    cov = defaultdict(list)
    if not coords_path or not os.path.exists(coords_path):
        return cov
    for line in open(coords_path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            f = line.split()
        if len(f) < 13:
            continue
        try:
            s2, e2 = int(f[2]), int(f[3])
        except ValueError:
            continue
        qname = f[-1]
        lo, hi = min(s2, e2), max(s2, e2)
        cov[qname].append((lo - 1, hi))
    for k in cov:
        cov[k].sort()
    return cov


def coverage_fraction(gene, covered):
    """Fraction of a gene's positions that the alignment covers.

    `covered` is either a set of 0-based positions (reference side) or a list of
    (start, end) half-open intervals (query side)."""
    s = covered.get(gene["contig"])
    n = gene["end"] - gene["start"] + 1
    if not s or n <= 0:
        return 0.0
    gs, ge = gene["start"] - 1, gene["end"]
    if isinstance(s, list):
        hit = 0
        for a, b in s:
            ov = min(ge, b) - max(gs, a)
            if ov > 0:
                hit += ov
            if a >= ge:
                break
        hit = min(hit, n)
    else:
        hit = sum(1 for p in range(gs, ge) if p in s)
    return 100.0 * hit / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--qry-gff", required=True)
    ap.add_argument("--ref-faa")
    ap.add_argument("--qry-faa")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--bam", default=None,
                    help="asm5 BAM; defaults to the pipeline location")
    ap.add_argument("--absent-threshold", type=float, default=80.0)
    ap.add_argument("--coords", default=None,
                    help="show-coords -rclTHd of the 1-to-1 delta, used for the "
                         "query-side coverage")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    bam = a.bam
    if bam is None:
        # pipeline default
        bam = os.path.join(os.path.dirname(a.outdir), "..", "..", "work", "bam",
                           "CCUG4856_to_NCTC9343.asm5.bam")
        bam = os.path.normpath(bam)
    if not os.path.exists(bam):
        print(f"[genecontent] BAM not found ({bam}); cannot assess presence",
              file=sys.stderr)
        return

    ref_genes = parse_gff_cds(a.ref_gff)
    qry_genes = parse_gff_cds(a.qry_gff, key_attr="locus_tag")
    print(f"[genecontent] reference CDS: {len(ref_genes)}, query CDS: {len(qry_genes)}")

    cov_ref, _ = walk_bam(bam)
    # query coverage must come from genome-coordinate blocks, not from SEQ indices
    cov_qry = block_coverage(a.coords) if a.coords else {}
    print(f"[genecontent] covered reference positions: "
          f"{ {k: len(v) for k, v in cov_ref.items()} }")
    print(f"[genecontent] covered query positions    : "
          f"{ {k: len(v) for k, v in cov_qry.items()} }")

    # ---- reference genes in the query ------------------------------------
    ref_pres = []
    for g in ref_genes:
        g = dict(g); g["coverage_pct"] = round(coverage_fraction(g, cov_ref), 2)
        ref_pres.append(g)
    with open(f"{a.outdir}/reference_gene_presence.tsv", "w") as fh:
        fh.write("locus_tag\tgene\tproduct\tcontig\tstart\tend\tstrand\t"
                 "query_coverage_pct\tstatus\n")
        for g in ref_pres:
            st = ("present" if g["coverage_pct"] >= a.absent_threshold
                  else "partial" if g["coverage_pct"] > 20 else "absent")
            fh.write(f"{g['locus']}\t{g['gene']}\t{g['product']}\t{g['contig']}\t"
                     f"{g['start']}\t{g['end']}\t{g['strand']}\t"
                     f"{g['coverage_pct']}\t{st}\n")

    absent = [g for g in ref_pres if g["coverage_pct"] < a.absent_threshold]
    with open(f"{a.outdir}/absent_from_query.tsv", "w") as fh:
        fh.write("locus_tag\tgene\tproduct\tcontig\tstart\tend\t"
                 "query_coverage_pct\n")
        for g in sorted(absent, key=lambda x: x["coverage_pct"]):
            fh.write(f"{g['locus']}\t{g['gene']}\t{g['product']}\t{g['contig']}\t"
                     f"{g['start']}\t{g['end']}\t{g['coverage_pct']}\n")

    # ---- query genes in the reference ------------------------------------
    qry_pres = []
    for g in qry_genes:
        g = dict(g); g["coverage_pct"] = round(coverage_fraction(g, cov_qry), 2)
        qry_pres.append(g)
    with open(f"{a.outdir}/query_gene_presence.tsv", "w") as fh:
        fh.write("locus_tag\tgene\tproduct\tcontig\tstart\tend\tstrand\t"
                 "reference_coverage_pct\tstatus\n")
        for g in qry_pres:
            st = ("present" if g["coverage_pct"] >= a.absent_threshold
                  else "partial" if g["coverage_pct"] > 20 else "absent")
            fh.write(f"{g['locus']}\t{g['gene']}\t{g['product']}\t{g['contig']}\t"
                     f"{g['start']}\t{g['end']}\t{g['strand']}\t"
                     f"{g['coverage_pct']}\t{st}\n")

    q_specific = [g for g in qry_pres if g["coverage_pct"] < a.absent_threshold]
    with open(f"{a.outdir}/query_specific_genes.tsv", "w") as fh:
        fh.write("locus_tag\tgene\tproduct\tcontig\tstart\tend\t"
                 "reference_coverage_pct\n")
        for g in sorted(q_specific, key=lambda x: x["coverage_pct"]):
            fh.write(f"{g['locus']}\t{g['gene']}\t{g['product']}\t{g['contig']}\t"
                     f"{g['start']}\t{g['end']}\t{g['coverage_pct']}\n")

    # ---- optional protein-level identity ---------------------------------
    n_ident = 0
    if a.ref_faa and a.qry_faa and os.path.exists(a.ref_faa) and os.path.exists(a.qry_faa):
        diamond = subprocess.run(["bash", "-lc", "command -v diamond"],
                                 capture_output=True, text=True).stdout.strip()
        blastp = subprocess.run(["bash", "-lc", "command -v blastp"],
                                capture_output=True, text=True).stdout.strip()
        if diamond or blastp:
            with tempfile.TemporaryDirectory() as td:
                out_tsv = os.path.join(td, "hits.tsv")
                if diamond:
                    db = os.path.join(td, "refdb")
                    cmd = (f"diamond makedb --in '{a.ref_faa}' -d '{db}' >/dev/null 2>&1 && "
                           f"diamond blastp -q '{a.qry_faa}' -d '{db}' -o '{out_tsv}' "
                           f"-f 6 qseqid sseqid pident length qlen slen evalue bitscore "
                           f"--max-target-seqs 1 --threads {a.threads} --quiet")
                else:
                    db = os.path.join(td, "refdb")
                    cmd = (f"makeblastdb -in '{a.ref_faa}' -dbtype prot -out '{db}' >/dev/null 2>&1 && "
                           f"blastp -query '{a.qry_faa}' -db '{db}' -out '{out_tsv}' "
                           f"-outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' "
                           f"-max_target_seqs 1 -num_threads {a.threads}")
                rc = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
                if rc.returncode == 0 and os.path.exists(out_tsv):
                    with open(f"{a.outdir}/shared_gene_identity.tsv", "w") as fh:
                        fh.write("query_protein\treference_protein\tpident\t"
                                 "aln_length\tquery_len\tsubject_len\tevalue\t"
                                 "bitscore\tquery_cov_pct\n")
                        for line in open(out_tsv):
                            f = line.rstrip("\n").split("\t")
                            if len(f) < 8:
                                continue
                            ql = float(f[4]) if f[4] not in ("", ".") else 0
                            al = float(f[3])
                            cov = 100.0 * al / ql if ql else 0
                            fh.write("\t".join(f) + f"\t{cov:.1f}\n")
                            n_ident += 1
                    print(f"[genecontent] protein identity for {n_ident} query proteins")

    summary = dict(
        reference_cds=len(ref_genes), query_cds=len(qry_genes),
        reference_genes_absent_or_partial_in_query=len(absent),
        query_genes_absent_or_partial_in_reference=len(q_specific),
        absent_from_query=[dict(locus=g["locus"], product=g["product"],
                                coverage=g["coverage_pct"]) for g in absent],
        query_specific=[dict(locus=g["locus"], product=g["product"],
                             coverage=g["coverage_pct"]) for g in q_specific],
        protein_identity_rows=n_ident,
        absent_threshold_pct=a.absent_threshold,
    )
    with open(f"{a.outdir}/gene_content_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[genecontent] reference genes absent/partial in query: {len(absent)}; "
          f"query genes absent/partial in reference: {len(q_specific)}")


if __name__ == "__main__":
    main()
