#!/usr/bin/env python3
"""
transfer_annotation.py -- derive a full query-genome annotation by projecting
every reference gene through the assembly-to-reference alignment.

Why this is legitimate here
  The two genomes are 99.99 % identical along a single collinear alignment (one
  88 kb inversion aside), so the reference's curated RefSeq annotation --
  products, GO terms, EC numbers, locus tags, pseudogene calls -- transfers gene
  by gene.  The transfer goes through the alignment, so every projected feature
  carries the query's own coordinates *and* the query's own sequence: genes that
  a variant touches are automatically written out with the corrected allele.

Implementation notes
  * The query allele SEQUENCE is read from each alignment record's own SEQ.  The
    asm5 alignment contains supplementary records with hard clips, so SEQ indices
    are NOT contig coordinates; using the record's own bases sidesteps that
    entirely and handles reverse-strand alignments correctly.
  * The query COORDINATES written into the GFF come from the 1-to-1 NUCmer
    blocks, which are expressed in true genome coordinates.

Outputs (in --outdir)
  query.gff3              projected annotation in query coordinates
  query.faa               projected query proteome
  projection_table.tsv    per gene: coordinates, protein identity, ORF status
  projection_summary.json
"""
import argparse
import json
import os
import sys
from collections import OrderedDict, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotate_variants import (read_fasta, revcomp, translate,  # noqa: E402
                               parse_gff, project_region)


def parse_blocks(coords_path):
    """1-to-1 NUCmer blocks, strand-aware, in true genome coordinates."""
    blocks = []
    if not coords_path or not os.path.exists(coords_path):
        return blocks
    for line in open(coords_path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            f = line.split()
        if len(f) < 13:
            continue
        try:
            s1, e1, s2, e2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
            qdir = int(f[-3])
        except ValueError:
            continue
        blocks.append(dict(ref=f[-2], rs=s1, re=e1, qry=f[-1], qs=s2, qe=e2,
                           strand=1 if qdir >= 0 else -1))
    return blocks


class CoordMap:
    """Map reference positions to query positions using the 1-to-1 blocks."""

    def __init__(self, blocks):
        self.by_ref = defaultdict(list)
        for b in blocks:
            self.by_ref[b["ref"]].append(b)
        for k in self.by_ref:
            self.by_ref[k].sort(key=lambda b: b["rs"])

    def lookup(self, contig, pos):
        for b in self.by_ref.get(contig, []):
            if b["rs"] <= pos <= b["re"]:
                if b["strand"] > 0:
                    q = b["qs"] + (pos - b["rs"])
                else:
                    q = b["qe"] - (pos - b["rs"])
                return b["qry"], q, b["strand"]
        return None, None, 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--ref-fna", required=True)
    ap.add_argument("--qry-fna", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--coords", default=None,
                    help="show-coords -rclTHd of the 1-to-1 delta")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-projection", type=float, default=0.9)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    import pysam
    bam = pysam.AlignmentFile(a.bam, "rb")
    ref = read_fasta(a.ref_fna)
    blocks = parse_blocks(a.coords)
    cmap = CoordMap(blocks)
    if not blocks:
        print("[transfer] WARNING: no coordinate blocks supplied; the query GFF "
              "will carry no coordinates", file=sys.stderr)

    feats = parse_gff(a.ref_gff)
    cds = [f for f in feats if f["type"] == "CDS"]
    others = [f for f in feats if f["type"] != "CDS"]

    prot_f = os.path.join(a.outdir, "query.faa")
    gff_f = os.path.join(a.outdir, "query.gff3")
    tab_f = os.path.join(a.outdir, "projection_table.tsv")
    stats = defaultdict(int)

    with open(prot_f, "w") as fp, open(gff_f, "w") as fg, open(tab_f, "w") as ft:
        fg.write("##gff-version 3\n")
        fg.write("# annotation transferred from the reference by alignment "
                 "projection\n")
        ft.write("locus_tag\tproduct\tref_contig\tref_start\tref_end\tstrand\t"
                 "qry_contig\tqry_start\tqry_end\tprojected_fraction\t"
                 "ref_protein_len\tqry_protein_len\tprotein_identical\t"
                 "query_internal_stops\tref_internal_stops\torf_status\t"
                 "effect_note\n")

        for g in cds:
            contig, s0, e0, strand = g["contig"], g["start"], g["end"], g["strand"]
            phase = g.get("phase", 0) or 0
            length = e0 - s0 + 1

            ref_gen = ref[contig][s0 - 1:e0]
            ref_tx = revcomp(ref_gen) if strand == "-" else ref_gen
            ref_tx = ref_tx[phase:] if phase else ref_tx
            ref_prot = translate(ref_tx)

            # ---- query allele, from the alignment records themselves -------
            ref2q, ins, covered = project_region(bam, contig, s0 - 1, e0)
            parts = []
            for p in range(s0 - 1, e0):
                parts.append(ref2q.get(p, "N"))
                if p in ins:
                    parts.extend(ins[p])
            qgen = "".join(parts)
            qtx = revcomp(qgen) if strand == "-" else qgen
            qtx = qtx[phase:] if phase else qtx
            qprot = translate(qtx)

            ncov = len([p for p in range(s0 - 1, e0) if p in covered])
            frac = ncov / length if length else 0.0

            # ---- query coordinates from the 1-to-1 blocks ------------------
            qname = qstart = qend = ""
            if blocks:
                qn1, qa, _ = cmap.lookup(contig, s0)
                qn2, qb, _ = cmap.lookup(contig, e0)
                if qn1 and qn2 and qn1 == qn2:
                    qname = qn1
                    qstart, qend = min(qa, qb), max(qa, qb)

            ref_c = ref_prot[:-1] if ref_prot.endswith("*") else ref_prot
            qry_c = qprot[:-1] if qprot.endswith("*") else qprot
            identical = (ref_c == qry_c and ref_c != "")
            q_stops = qry_c.count("*")
            r_stops = ref_c.count("*")

            if frac < a.min_projection:
                note, orf = "partial_projection", ""
                stats["partial"] += 1
            elif identical:
                note = "identical_to_reference"
                orf = "intact_in_both" if r_stops == 0 else "disrupted_in_both"
                stats["identical"] += 1
            else:
                note = "protein_differs"
                if r_stops == 0 and q_stops == 0:
                    orf = "intact_in_both"
                elif r_stops > 0 and q_stops == 0:
                    orf = "query_ORF_restored"
                elif r_stops == 0 and q_stops > 0:
                    orf = "query_ORF_disrupted"
                else:
                    orf = "disrupted_in_both"
                stats["differs"] += 1
                stats[f"orf_{orf}"] += 1

            ft.write(f"{g['locus']}\t{g['product']}\t{contig}\t{s0}\t{e0}\t"
                     f"{strand}\t{qname}\t{qstart}\t{qend}\t{frac:.3f}\t"
                     f"{len(ref_c)}\t{len(qry_c)}\t"
                     f"{'yes' if identical else 'no'}\t{q_stops}\t{r_stops}\t"
                     f"{orf}\t{note}\n")

            if frac >= a.min_projection:
                fp.write(f">{g['locus']}|{g['gene']}|{g['product']}\n{qprot}\n")
                if qname and qstart != "":
                    attrs = [f"ID=cds-{g['locus']}"]
                    for k, key in (("locus", "locus_tag"), ("gene", "gene"),
                                   ("product", "product"),
                                   ("protein_id", "protein_id"),
                                   ("ec", "EC_number")):
                        if g.get(k):
                            attrs.append(f"{key}={g[k]}")
                    if g.get("go"):
                        attrs.append(f"Ontology_term={g['go']}")
                    if g.get("pseudo"):
                        attrs.append("pseudo=true")
                    attrs.append("transferred_from_reference=true")
                    fg.write(f"{qname}\tWGS_transfer\tCDS\t{qstart}\t{qend}\t.\t"
                             f"{strand}\t0\t{';'.join(attrs)}\n")
                    stats["written"] += 1

        for g in others:
            qn, qa, _ = cmap.lookup(g["contig"], g["start"])
            qn2, qb, _ = cmap.lookup(g["contig"], g["end"])
            if not (qn and qn2 and qn == qn2):
                continue
            fg.write(f"{qn}\tWGS_transfer\t{g['type']}\t{min(qa,qb)}\t"
                     f"{max(qa,qb)}\t.\t{g['strand']}\t0\t"
                     f"ID={g['type']}-{g['locus']};locus_tag={g['locus']};"
                     f"product={g['product']}\n")
            stats[f"nc_{g['type']}"] += 1

    bam.close()

    summary = dict(
        reference_CDS=len(cds),
        projected_and_written=stats.get("written", 0),
        identical_protein=stats.get("identical", 0),
        protein_differs=stats.get("differs", 0),
        query_ORF_restored=stats.get("orf_query_ORF_restored", 0),
        query_ORF_disrupted=stats.get("orf_query_ORF_disrupted", 0),
        sequence_differs_but_intact_in_both=stats.get("orf_intact_in_both", 0),
        disrupted_in_both=stats.get("orf_disrupted_in_both", 0),
        partial_projection=stats.get("partial", 0),
        non_coding_transferred=sum(v for k, v in stats.items()
                                   if k.startswith("nc_")),
        min_projection_fraction=a.min_projection,
    )
    with open(os.path.join(a.outdir, "projection_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print("[transfer] " + ", ".join(f"{k}={v}" for k, v in summary.items()))


if __name__ == "__main__":
    main()
