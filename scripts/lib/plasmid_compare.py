#!/usr/bin/env python3
"""
plasmid_compare.py -- put the reference plasmid and the query plasmid side by
side: length, orientation, nucleotide differences and the genes those
differences fall in.

Both replicons are located automatically (a FASTA header containing "plasmid",
otherwise the shortest contig of each assembly).  Alignment is done with NUCmer
when available and falls back to a direct position-by-position comparison.
"""
import argparse
import json
import os
import subprocess
import sys
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotate_variants import read_fasta, parse_gff, revcomp  # noqa: E402


def pick_replicons(path):
    """Return (chromosome_name, plasmid_name, {name: sequence}, {name: header})."""
    seqs, heads, name, chunks = OrderedDict(), {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                header = line[1:].strip()
                name = header.split()[0]
                heads[name] = header
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks)

    plasmid = None
    for n, h in heads.items():
        if "plasmid" in h.lower():
            plasmid = n
            break
    chrom = max(seqs, key=lambda n: len(seqs[n]))
    if plasmid is None:
        plasmid = min(seqs, key=lambda n: len(seqs[n]))
    return chrom, plasmid, seqs, heads


def run(cmd):
    return subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-fna", required=True)
    ap.add_argument("--qry-fna", required=True)
    ap.add_argument("--ref-gff")
    ap.add_argument("--bam")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    rchrom, rplas, rseqs, rheads = pick_replicons(a.ref_fna)
    qchrom, qplas, qseqs, qheads = pick_replicons(a.qry_fna)

    rep = OrderedDict()
    rep["reference_chromosome"] = dict(
        name=rchrom, length=len(rseqs[rchrom]), header=rheads[rchrom])
    rep["query_chromosome"] = dict(
        name=qchrom, length=len(qseqs[qchrom]), header=qheads[qchrom])
    rep["reference_plasmid"] = dict(
        name=rplas, length=len(rseqs[rplas]), header=rheads[rplas])
    rep["query_plasmid"] = dict(
        name=qplas, length=len(qseqs[qplas]), header=qheads[qplas])

    rp = rseqs[rplas].upper()
    qp_raw = qseqs[qplas].upper()

    # ---- circularity-aware comparison -------------------------------------
    # Plasmids are circular, so two identical plasmids can be deposited with
    # different origins and/or in opposite orientations.  A naive
    # position-by-position comparison of such a pair reports ~25 % identity
    # (i.e. chance) even though the molecules are identical.  The comparison
    # therefore searches for the rotation offset that best aligns the two
    # sequences, in both orientations.
    def best_circular_identity(a, b):
        """Return (identity_percent, offset, orientation) for circular a vs b."""
        L = len(a)
        if L != len(b):
            return None, None, None
        k = 31
        best = (0.0, 0, "forward")
        for orient, bb in (("forward", b), ("reverse_complement", revcomp(b))):
            # locate a unique k-mer of `a` inside `bb` to fix the offset
            offsets = []
            for start in range(0, min(L - k, 200000), max(1, L // 200)):
                probe = a[start:start + k]
                if probe.count("N"):
                    continue
                pos = bb.find(probe)
                if pos >= 0 and bb.find(probe, pos + 1) < 0:   # unique
                    offsets.append((pos - start) % L)
                if len(offsets) >= 5:
                    break
            cand = offsets or [0]
            for off in cand:
                rolled = bb[off:] + bb[:off]
                same = sum(1 for x, y in zip(a, rolled) if x == y)
                ident = 100.0 * same / L
                if ident > best[0]:
                    best = (ident, off, orient)
        return best

    ident, offset, orientation = best_circular_identity(rp, qp_raw)
    if ident is None:                      # different lengths: no circular match
        orientation = "indeterminate"
        offset = None
    # store the query in the reference's frame so the sequences are directly
    # comparable and the written FASTA is meaningful
    if orientation == "reverse_complement":
        qp = revcomp(qp_raw)
    elif orientation == "forward":
        qp = qp_raw
    else:
        qp = qp_raw
    if offset:
        qp = qp[offset:] + qp[:offset]

    rep["orientation"] = orientation
    rep["circular_rotation_offset_bp"] = offset
    rep["length_equal"] = len(rseqs[rplas]) == len(qseqs[qplas])
    rep["identity_circular_pct"] = None if ident is None else round(ident, 3)
    rep["identity_direct_pct"] = None if ident is None else round(ident, 3)
    rep["circular_note"] = (
        "plasmids are circular: the two molecules are compared after allowing "
        "for an arbitrary rotation offset and both orientations, so a reported "
        "identity of ~100 % means the molecules are the same sequence even if "
        "they were deposited from different origins")

    # ---- write the plasmid pair for downstream tools ----------------------
    rf = os.path.join(a.outdir, "reference_plasmid.fna")
    qf = os.path.join(a.outdir, "query_plasmid.fna")
    with open(rf, "w") as fh:
        fh.write(f">{rplas}\n")
        for i in range(0, len(rp), 70):
            fh.write(rp[i:i + 70] + "\n")
    # the query plasmid is written in the reference's frame (orientation and
    # rotation applied) so the two files are directly comparable
    with open(qf, "w") as fh:
        fh.write(f">{qplas}_in_reference_frame\n")
        for i in range(0, len(qp), 70):
            fh.write(qp[i:i + 70] + "\n")

    # ---- differences ------------------------------------------------------
    diffs = []
    nucmer = run("command -v nucmer").stdout.strip()
    show_snps = run("command -v show-snps").stdout.strip()
    if nucmer and show_snps:
        prefix = os.path.join(a.outdir, "plasmid_nucmer")
        run(f"nucmer --maxmatch -c 50 -l 20 -p '{prefix}' '{rf}' '{qf}'")
        delta = f"{prefix}.delta"
        if os.path.exists(delta):
            out = run(f"show-snps -rlTHC -I 70 '{delta}'")
            with open(os.path.join(a.outdir, "plasmid_nucmer_snps.tsv"), "w") as fh:
                fh.write("ref_pos\tref_base\tqry_base\tqry_pos\tbuff\tdist\t"
                         "lenR\tlenQ\tfrm\ttags\n")
                fh.write(out.stdout)
            for line in out.stdout.splitlines():
                f = line.split("\t")
                if len(f) >= 3:
                    diffs.append(dict(pos=int(f[0]), ref=f[1], alt=f[2],
                                      qry_pos=f[3] if len(f) > 3 else "",
                                      buff=f[4] if len(f) > 4 else "",
                                      dist=f[5] if len(f) > 5 else ""))
    if not diffs and len(rp) == len(qp):
        for i, (x, y) in enumerate(zip(rp, qp), 1):
            if x != y:
                diffs.append(dict(pos=i, ref=x, alt=y, qry_pos=i, buff="", dist=""))

    # ---- gene context -----------------------------------------------------
    genes = []
    if a.ref_gff and os.path.exists(a.ref_gff):
        for g in parse_gff(a.ref_gff):
            if g["contig"] == rplas and g["type"] == "CDS":
                genes.append(g)

    def gene_at(pos):
        for g in genes:
            if g["start"] <= pos <= g["end"]:
                return g
        return None

    with open(os.path.join(a.outdir, "plasmid_differences.tsv"), "w") as fh:
        fh.write("position\tref_base\tqry_base\tquery_position\tbuff\tdist\t"
                 "region\tlocus_tag\tgene\tproduct\n")
        for d in diffs:
            g = gene_at(d["pos"])
            fh.write(f"{d['pos']}\t{d['ref']}\t{d['alt']}\t{d['qry_pos']}\t"
                     f"{d['buff']}\t{d['dist']}\t"
                     f"{'CDS' if g else 'intergenic'}\t"
                     f"{g['locus'] if g else ''}\t{g['gene'] if g else ''}\t"
                     f"{g['product'] if g else ''}\n")

    # ---- gene content -----------------------------------------------------
    with open(os.path.join(a.outdir, "plasmid_genes.tsv"), "w") as fh:
        fh.write("locus_tag\tgene\tproduct\tstart\tend\tstrand\n")
        for g in sorted(genes, key=lambda x: x["start"]):
            fh.write(f"{g['locus']}\t{g['gene']}\t{g['product']}\t{g['start']}\t"
                     f"{g['end']}\t{g['strand']}\n")

    rep["nucleotide_differences"] = len(diffs)
    rep["reference_plasmid_genes"] = len(genes)
    rep["differences_in_genes"] = sum(1 for d in diffs if gene_at(d["pos"]))
    with open(os.path.join(a.outdir, "plasmid_comparison.json"), "w") as fh:
        json.dump(rep, fh, indent=2)

    print(f"[plasmid] ref {rplas} {len(rp):,} bp vs qry {qplas} {len(qp):,} bp; "
          f"orientation={orientation}; identity="
          f"{rep['identity_direct_pct']}%; differences={len(diffs)}")


if __name__ == "__main__":
    main()
