#!/usr/bin/env python3
"""
assembly_stats.py -- reference-free assembly quality metrics for a FASTA file.

Reports contig count, total length, N50/L50/N75/L75/N90, GC%, ambiguous bases,
longest/shortest contig, and per-contig detail.  Written as a dependency-free
fallback/complement to QUAST so that QC is always produced.

Usage: assembly_stats.py <genome.fna> <label> <outdir>
"""
import sys
import os
import json
from collections import OrderedDict


def read_fasta(path):
    """Yield (name, sequence) -- tolerates gzip transparently."""
    opener = None
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    else:
        opener = open
    name, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def nx(lengths_sorted_desc, total, frac):
    """Nx: length of the contig at which the cumulative length reaches frac*total."""
    target = total * frac
    cum = 0
    for i, ln in enumerate(lengths_sorted_desc, start=1):
        cum += ln
        if cum >= target:
            return ln, i
    return 0, 0


AMBIG = set("NRYKMSWBDHVnrykmswbdhv")


def stats_for(seqs):
    lengths = [len(s) for _, s in seqs]
    lengths.sort(reverse=True)
    total = sum(lengths)
    n_contigs = len(lengths)
    gc = amb = 0
    at = 0
    for _, s in seqs:
        u = s.upper()
        gc += u.count("G") + u.count("C")
        amb += sum(u.count(c) for c in AMBIG)
        at += u.count("A") + u.count("T")
    unambiguous = total - amb
    out = OrderedDict()
    out["num_contigs"] = n_contigs
    out["total_length_bp"] = total
    out["largest_contig_bp"] = lengths[0] if lengths else 0
    out["smallest_contig_bp"] = lengths[-1] if lengths else 0
    out["mean_contig_bp"] = round(total / n_contigs, 1) if n_contigs else 0
    for label, frac in (("N50", 0.50), ("N75", 0.75), ("N90", 0.90)):
        val, idx = nx(lengths, total, frac)
        out[f"{label}_bp"] = val
        out[f"L{label[1:]}"] = idx
    out["GC_percent"] = round(100.0 * gc / unambiguous, 2) if unambiguous else 0.0
    out["AT_percent"] = round(100.0 * at / unambiguous, 2) if unambiguous else 0.0
    out["ambiguous_bases"] = amb
    out["ambiguous_percent"] = round(100.0 * amb / total, 4) if total else 0.0
    out["contigs_gt_1kb"] = sum(1 for x in lengths if x > 1000)
    out["contigs_gt_10kb"] = sum(1 for x in lengths if x > 10000)
    # assembly "continuity" descriptors
    out["assembly_level_guess"] = (
        "Complete Genome" if (n_contigs <= 5 and out["contigs_gt_10kb"] == n_contigs)
        else "Draft"
    )
    out["_lengths"] = lengths
    return out


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    fasta, label, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)

    seqs = list(read_fasta(fasta))
    st = stats_for(seqs)
    lengths = st.pop("_lengths")

    # ---- summary TSV ------------------------------------------------------
    tsv = os.path.join(outdir, f"{label}_assembly_stats.tsv")
    with open(tsv, "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in st.items():
            fh.write(f"{k}\t{v}\n")

    # ---- per-contig TSV ---------------------------------------------------
    per = os.path.join(outdir, f"{label}_per_contig.tsv")
    with open(per, "w") as fh:
        fh.write("contig\tlength_bp\tGC_percent\tambiguous\n")
        for name, s in sorted(seqs, key=lambda x: -len(x[1])):
            u = s.upper()
            g = 100.0 * (u.count("G") + u.count("C")) / len(s) if s else 0
            a = sum(u.count(c) for c in AMBIG)
            fh.write(f"{name}\t{len(s)}\t{g:.2f}\t{a}\n")

    # ---- cumulative-length curve (for the Nx plot) -------------------------
    curve = os.path.join(outdir, f"{label}_cumulative_length.tsv")
    total = sum(lengths)
    cum = 0
    with open(curve, "w") as fh:
        fh.write("contig_rank\tlength_bp\tcumulative_bp\tcumulative_percent\n")
        for i, ln in enumerate(lengths, 1):
            cum += ln
            fh.write(f"{i}\t{ln}\t{cum}\t{100.0*cum/total:.4f}\n")

    # ---- JSON (for the plotting / report layer) ---------------------------
    with open(os.path.join(outdir, f"{label}_assembly_stats.json"), "w") as fh:
        json.dump(st, fh, indent=2)

    print(f"[assembly_stats] {label}: {st['num_contigs']} contigs, "
          f"{st['total_length_bp']:,} bp, N50={st['N50_bp']:,}, "
          f"GC={st['GC_percent']}%")


if __name__ == "__main__":
    main()
