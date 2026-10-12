#!/usr/bin/env python3
"""
structural_variants.py -- large-scale differences between the reference and the
query genome, derived from the 1-to-1 NUCmer alignment blocks.

It distinguishes three very different things that a naive position-by-position
comparison would conflate:

  * ORIGIN SHIFT (re-linearisation)  -- the deposited query sequence simply
    starts at a different position of the same circular chromosome.  This is a
    file/coordinate convention, NOT a biological difference.  Reporting it
    explicitly prevents the ~100 % of "positions" it touches from being mistaken
    for variation.

  * INVERSIONS -- a block that aligns to the reference in the reverse
    orientation and is bounded by forward-aligned flanks.  Genuine structural
    variation.

  * TRANSLOCATIONS -- a block whose reference and query neighbours disagree.

Outputs (in --outdir)
  alignment_blocks.tsv        every 1-to-1 block with orientation
  structural_variants.tsv     the rearrangements that were called
  origin_shift.tsv            the re-linearisation offset per replicon
  structural_summary.json
"""
import argparse
import json
import os
import statistics
import sys
from collections import OrderedDict, Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotate_variants import parse_gff  # noqa: E402


def parse_coords(path):
    """show-coords -rclTHd output -> list of blocks."""
    blocks = []
    if not os.path.exists(path):
        return blocks
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            f = line.split()
        if len(f) < 13:
            continue
        try:
            s1, e1, s2, e2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
            len1, len2 = int(f[4]), int(f[5])
            pid = float(f[6])
            lenR, lenQ = int(f[7]), int(f[8])
        except ValueError:
            continue
        # with -d the last four fields are: ref_dir, qry_dir, ref_tag, qry_tag
        try:
            qdir = int(f[-3])
        except ValueError:
            qdir = 1
        blocks.append(dict(s1=s1, e1=e1, s2=s2, e2=e2, len1=len1, len2=len2,
                           pid=pid, lenR=lenR, lenQ=lenQ,
                           strand="+" if qdir >= 0 else "-",
                           ref_tag=f[-2], qry_tag=f[-1]))
    return blocks


def read_fasta_lengths(path):
    lens, name = OrderedDict(), None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    lens[name] = n
                name = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name:
        lens[name] = n
    return lens


def genes_at(genes, contig, start, end, window=5000):
    pad_s, pad_e = start - window, end + window
    out = []
    for g in genes:
        if g["contig"] != contig:
            continue
        if g["end"] >= pad_s and g["start"] <= pad_e:
            out.append(g)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coords", required=True, help="show-coords -rclTHd output (1-to-1)")
    ap.add_argument("--coords-all", default=None,
                    help="show-coords over ALL alignments, used to test whether a "
                         "reverse block is really an inversion or just an inverted repeat")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--qry", required=True)
    ap.add_argument("--gff")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-inversion", type=int, default=10000,
                    help="minimum size for a confident inversion call")
    ap.add_argument("--min-block", type=int, default=500)
    ap.add_argument("--offset-tolerance", type=int, default=1000,
                    help="bp tolerance when clustering origin-shift offsets")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    blocks = [b for b in parse_coords(a.coords) if b["len1"] >= a.min_block]
    all_blocks = parse_coords(a.coords_all) if a.coords_all else []
    ref_len = read_fasta_lengths(a.ref)
    qry_len = read_fasta_lengths(a.qry)
    genes = parse_gff(a.gff) if a.gff and os.path.exists(a.gff) else []

    # ---- write all blocks -------------------------------------------------
    with open(f"{a.outdir}/alignment_blocks.tsv", "w") as fh:
        fh.write("ref_contig\tref_start\tref_end\tqry_contig\tqry_start\t"
                 "qry_end\tstrand\tlength_bp\tidentity_pct\n")
        for b in sorted(blocks, key=lambda x: (x["ref_tag"], x["s1"])):
            fh.write(f"{b['ref_tag']}\t{b['s1']}\t{b['e1']}\t{b['qry_tag']}\t"
                     f"{b['s2']}\t{b['e2']}\t{b['strand']}\t{b['len1']}\t"
                     f"{b['pid']:.2f}\n")

    # ---- origin shift -----------------------------------------------------
    # On a re-linearised (rotated) genome every long forward block shares one
    # offset between reference and query coordinates:
    #     qry_pos = (ref_pos - offset) mod replicon_length
    # Block boundaries differ by a few bp, so offsets are clustered with a
    # tolerance rather than compared exactly.
    origin_rows = []
    offsets_by_ref = {}
    for rc, rl in ref_len.items():
        fwd = [b for b in blocks
               if b["ref_tag"] == rc and b["strand"] == "+" and b["len1"] >= 20000]
        if not fwd:
            continue
        raw = [((b["s2"] - b["s1"]) % rl) for b in fwd]
        clusters = []                       # list of [centre, total_bp, members]
        for off, b in zip(raw, fwd):
            placed = False
            for c in clusters:
                d = min(abs(off - c[0]), rl - abs(off - c[0]))
                if d <= a.offset_tolerance:
                    c[1] += b["len1"]; c[2].append(off)
                    # keep a running representative value
                    c[0] = int(statistics.median(c[2]))
                    placed = True
                    break
            if not placed:
                clusters.append([off, b["len1"], [off]])
        clusters.sort(key=lambda c: -c[1])
        best = clusters[0]
        covered = best[1] / rl * 100.0
        # a rotated genome has ONE cluster covering essentially everything;
        # an unrotated genome has offset 0 covering everything.
        shifted = covered > 50 and not (best[0] <= a.offset_tolerance or
                                        rl - best[0] <= a.offset_tolerance)
        offsets_by_ref[rc] = best[0] if shifted else 0
        origin_rows.append(dict(
            ref_contig=rc, ref_length=rl, n_forward_blocks=len(fwd),
            n_offset_clusters=len(clusters),
            dominant_offset_bp=best[0],
            n_blocks_at_dominant_offset=len(best[2]),
            pct_reference_bp_at_dominant_offset=round(covered, 2),
            origin_shifted=bool(shifted),
            interpretation=(
                f"query sequence is re-linearised (rotated): reference position 1 "
                f"corresponds to query position {best[0] or rl:,} of "
                f"{qry_len.get(fwd[0]['qry_tag'], 0):,} bp.  This is a deposited-"
                f"sequence coordinate convention, NOT a biological difference."
                if shifted else
                "no systematic origin shift detected")))

    with open(f"{a.outdir}/origin_shift.tsv", "w") as fh:
        fh.write("ref_contig\tref_length\tn_forward_blocks\tn_offset_clusters\t"
                 "dominant_offset_bp\tn_blocks_at_dominant_offset\t"
                 "pct_reference_bp_at_dominant_offset\torigin_shifted\t"
                 "interpretation\n")
        for r in origin_rows:
            fh.write("\t".join(str(r[k]) for k in
                               ["ref_contig", "ref_length", "n_forward_blocks",
                                "n_offset_clusters", "dominant_offset_bp",
                                "n_blocks_at_dominant_offset",
                                "pct_reference_bp_at_dominant_offset",
                                "origin_shifted", "interpretation"]) + "\n")

    # ---- inversions / reverse blocks --------------------------------------
    def has_forward_copy(ref_tag, s1, e1, qry_tag, min_frac=0.8):
        """Is there a long *forward* alignment of this reference segment?

        If yes the reverse block is an inverted repeat rather than a true
        inversion (the segment is present in both orientations)."""
        span = e1 - s1 + 1
        for b in all_blocks:
            if b["ref_tag"] != ref_tag or b["strand"] != "+":
                continue
            ov = min(b["e1"], e1) - max(b["s1"], s1) + 1
            if ov > 0 and ov / span >= min_frac and b["pid"] >= 90:
                return True
        return False

    sv = []
    by_ref = defaultdict(list)
    for b in blocks:
        by_ref[b["ref_tag"]].append(b)

    for rc, bl in by_ref.items():
        bl.sort(key=lambda x: x["s1"])
        rl = ref_len.get(rc, 0)
        for b in bl:
            if b["strand"] != "-":
                continue
            frac = b["len1"] / rl if rl else 0
            left = [x for x in bl if x["strand"] == "+" and x["e1"] <= b["s1"] + 1000]
            right = [x for x in bl if x["strand"] == "+" and x["s1"] >= b["e1"] - 1000]
            flanked = bool(left) and bool(right)
            repeat = has_forward_copy(b["ref_tag"], b["s1"], b["e1"], b["qry_tag"])
            if frac > 0.8:
                typ = "whole_replicon_reverse_orientation"
            elif b["len1"] >= a.min_inversion and flanked and not repeat:
                typ = "inversion"
            elif b["len1"] >= a.min_inversion and repeat:
                typ = "inverted_repeat"
            else:
                typ = "small_reverse_block"
            sv.append(dict(
                type=typ, ref_contig=rc, ref_start=b["s1"], ref_end=b["e1"],
                qry_contig=b["qry_tag"], qry_start=b["s2"], qry_end=b["e2"],
                length_bp=b["len1"], identity_pct=round(b["pid"], 2),
                left_flank_closed=bool(left), right_flank_closed=bool(right),
                forward_copy_elsewhere=bool(repeat),
                genes_at_breakpoints=""))

    # ---- annotate breakpoints with genes ---------------------------------
    for s in sv:
        hits = genes_at(genes, s["ref_contig"], s["ref_start"], s["ref_end"])
        s["genes_at_breakpoints"] = ";".join(
            f"{g['locus']}({g['gene'] or g['product'][:30]})" for g in hits[:12])

    with open(f"{a.outdir}/structural_variants.tsv", "w") as fh:
        cols = ["type", "ref_contig", "ref_start", "ref_end", "qry_contig",
                "qry_start", "qry_end", "length_bp", "identity_pct",
                "left_flank_closed", "right_flank_closed",
                "forward_copy_elsewhere", "genes_at_breakpoints"]
        fh.write("\t".join(cols) + "\n")
        for s in sorted(sv, key=lambda x: -x["length_bp"]):
            fh.write("\t".join(str(s.get(c, "")) for c in cols) + "\n")

    inversions = [s for s in sv if s["type"] == "inversion"]
    summary = dict(
        n_alignment_blocks=len(blocks),
        total_aligned_reference_bp=sum(b["len1"] for b in blocks),
        n_inversions=len(inversions),
        n_inverted_repeats=len([s for s in sv if s["type"] == "inverted_repeat"]),
        n_small_reverse_blocks=len([s for s in sv if s["type"] == "small_reverse_block"]),
        largest_inversion_bp=max([s["length_bp"] for s in inversions], default=0),
        origin_shift=[r for r in origin_rows if r["origin_shifted"]],
        structural_variants=sv,
    )
    with open(f"{a.outdir}/structural_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[sv] {len(blocks)} alignment blocks; "
          f"{summary['n_inversions']} confirmed inversion(s), "
          f"largest {summary['largest_inversion_bp']:,} bp; "
          f"{len(summary['origin_shift'])} replicon(s) with an origin shift")


if __name__ == "__main__":
    main()
