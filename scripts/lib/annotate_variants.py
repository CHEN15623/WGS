#!/usr/bin/env python3
"""
annotate_variants.py -- functional consequence annotation of the variant set.

Two independent layers are produced:

  A. VARIANT-CENTRIC (reference coordinates)
     Every variant is located relative to the reference annotation: which gene,
     which codon, which amino-acid change, which effect class, conservation
     (BLOSUM62) of the substitution, and distance to the nearest gene for
     intergenic calls.

  B. GENE-CENTRIC (allele-aware, alignment-projected)
     For every gene touched by a variant the *actual query allele* is recovered
     by projecting the reference gene interval through the minimap2 asm5 CIGAR
     onto the query contig.  Both alleles are translated, so the true protein
     consequence (including frameshifts and downstream premature stops) is
     reported rather than merely predicted from the local codon.

Outputs (in --outdir)
  variant_annotation.tsv     one row per variant, full consequence annotation
  gene_impact.tsv            one row per affected gene, reference vs query protein
  affected_genes.faa         query (mutant) protein sequences
  reference_genes.faa        matching reference protein sequences
  annotation_summary.json    headline counts
"""
import argparse
import json
import os
import sys
from collections import OrderedDict, Counter, defaultdict

# --------------------------------------------------------------------------- #
#  sequence utilities
# --------------------------------------------------------------------------- #
COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv",
                     "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(s):
    return s.translate(COMP)[::-1]


CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(seq):
    seq = seq.upper()
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X")
                   for i in range(0, len(seq) - len(seq) % 3, 3))


# ---- BLOSUM62 (only the score lookup is needed) ---------------------------
_B62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62 = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
_B62_MAT = {}
_rows = [r for r in _B62.strip().split("\n") if r.strip()]
for _i, _r in enumerate(_rows):
    for _j, _v in enumerate(_r.split()):
        _B62_MAT[(_B62_ORDER[_i], _B62_ORDER[_j])] = int(_v)


def blosum62(a, b):
    if a in "*-X" or b in "*-X":
        return ""
    return _B62_MAT.get((a, b), "")


# --------------------------------------------------------------------------- #
#  parsing
# --------------------------------------------------------------------------- #

def read_fasta(path):
    seqs, name, chunks = OrderedDict(), None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]; chunks = []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks)
    return seqs


def parse_gff(path):
    """Return (features, genes).

    features : list of dicts for CDS / tRNA / rRNA / tmRNA / ncRNA
    genes    : {locus_tag: dict} for CDS, with the transcription-ordered
               genomic position map needed for codon arithmetic.
    """
    feats = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        ftype = f[2]
        if ftype not in ("CDS", "tRNA", "rRNA", "tmRNA", "ncRNA",
                         "RNase_P_RNA", "SRP_RNA"):
            continue
        d = {}
        for kv in f[8].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                d.setdefault(k, v)
        pseudo = (d.get("pseudo", "").lower() == "true"
                  or ftype == "pseudogene"
                  or "frameshift" in d.get("Note", "").lower()
                  or "internal stop" in d.get("Note", "").lower())
        feats.append(dict(
            contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
            phase=int(f[7]) if f[7] in "012" else 0, type=ftype,
            locus=d.get("locus_tag", d.get("ID", "")),
            gene=d.get("gene", d.get("Name", "")),
            product=d.get("product", ""),
            protein_id=d.get("protein_id", ""),
            pseudo=pseudo,
            note=d.get("Note", ""),
            ec=d.get("EC_number", ""),
            go=d.get("Ontology_term", ""),
        ))
    return feats


# --------------------------------------------------------------------------- #
#  CIGAR-based projection of reference intervals onto the query assembly
# --------------------------------------------------------------------------- #

def project_region(bam, contig, start0, end0):
    """Project reference interval [start0,end0) (0-based, half-open) to query.

    Returns (ref2qry_base, insertions, ref_covered_by_alignment).
      ref2qry_base  : {ref_pos0: query base}
      insertions    : {anchor_ref_pos0: [extra query bases]}  (insertion follows
                      the reference base at anchor_ref_pos0)
    """
    ref2q = {}
    ins = defaultdict(list)
    covered = set()
    try:
        reads = list(bam.fetch(contig, start0, end0))
    except Exception:
        return ref2q, ins, covered
    for read in reads:
        if read.is_unmapped or read.cigartuples is None:
            continue
        qseq = read.query_sequence
        if qseq is None:
            continue
        rp = read.reference_start
        qp = read.query_alignment_start
        for op, ln in read.cigartuples:
            if op in (0, 7, 8):                       # M / = / X
                lo = max(start0, rp)
                hi = min(end0, rp + ln)
                for p in range(lo, hi):
                    ref2q[p] = qseq[qp + (p - rp)]
                    covered.add(p)
                rp += ln; qp += ln
            elif op == 1:                             # I : query has extra bases
                anchor = rp - 1
                if start0 - 1 <= anchor < end0:
                    ins[anchor].extend(qseq[qp:qp + ln])
                qp += ln
            elif op == 2:                             # D : reference bases absent
                for p in range(max(start0, rp), min(end0, rp + ln)):
                    covered.add(p)
                rp += ln
            elif op == 3:                             # N
                rp += ln
            # op 4 (S, soft clip) must NOT advance qp: query_alignment_start
            # already excludes the leading soft-clipped bases, so advancing here
            # would count them twice and read the wrong query bases.
            # 5=H, 6=P -> nothing
    return ref2q, ins, covered


def project_allele(bam, contig, start0, end0, strand):
    """Return the query allele sequence for a reference interval, transcription
    oriented (i.e. as it appears in the protein)."""
    ref2q, ins, covered = project_region(bam, contig, start0, end0)
    parts = []
    for p in range(start0, end0):
        parts.append(ref2q.get(p, "N"))
        if p in ins:
            parts.extend(ins[p])
    seq = "".join(parts)
    if strand == "-":
        seq = revcomp(seq)
    return seq, covered


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", required=True,
                    help="variants_highconf.tsv (or variants_merged.tsv)")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--bam", help="asm5 BAM used to recover the query allele")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--updown", type=int, default=300,
                    help="bp window for upstream/downstream classification")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    genome = read_fasta(a.ref)
    feats = parse_gff(a.gff)

    # ---- gene lookup structures -------------------------------------------
    cds = [f for f in feats if f["type"] == "CDS"]
    noncds = [f for f in feats if f["type"] != "CDS"]
    by_contig = defaultdict(list)
    for f in cds:
        by_contig[f["contig"]].append(f)
    for c in by_contig:
        by_contig[c].sort(key=lambda x: x["start"])

    import bisect
    starts = {c: [f["start"] for f in v] for c, v in by_contig.items()}

    def find_overlapping(contig, pos):
        v = by_contig.get(contig, [])
        if not v:
            return []
        i = bisect.bisect_right(starts[contig], pos)
        hits = []
        for j in (i - 1, i, i - 2):
            if 0 <= j < len(v) and v[j]["start"] <= pos <= v[j]["end"]:
                hits.append(v[j])
        return hits

    def nearest_gene(contig, pos):
        v = by_contig.get(contig, [])
        best, bestd = None, None
        for f in v:
            if f["end"] < pos:
                d = pos - f["end"]
                side = "upstream" if f["strand"] == "+" else "downstream"
            elif f["start"] > pos:
                d = f["start"] - pos
                side = "downstream" if f["strand"] == "+" else "upstream"
            else:
                return f, 0, "inside"
            if bestd is None or d < bestd:
                best, bestd, bestside = f, d, side
        return (best, bestd, bestside) if best else (None, None, None)

    # ---- optional BAM for allele projection -------------------------------
    bam = None
    if a.bam and os.path.exists(a.bam):
        try:
            import pysam
            bam = pysam.AlignmentFile(a.bam, "rb")
        except Exception as exc:
            print(f"[annot] pysam/BAM unavailable ({exc}); "
                  f"gene-centric protein projection skipped", file=sys.stderr)

    # ---- load variants ----------------------------------------------------
    rows = []
    with open(a.variants) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            f += [""] * (len(hdr) - len(f))
            rows.append(dict(zip(hdr, f)))

    # ---- codon arithmetic helper -----------------------------------------
    def reference_cds(gene):
        """The reference CDS as it is translated: transcription oriented, so a
        minus-strand gene is reverse COMPLEMENTED (not merely reversed), and
        phase-trimmed."""
        s, e = gene["start"], gene["end"]
        raw = genome[gene["contig"]][s - 1:e]
        seq = revcomp(raw) if gene["strand"] == "-" else raw
        ph = gene.get("phase", 0) or 0
        return seq[ph:] if ph else seq

    def cds_offset(gene, pos):
        """0-based offset of a genomic position within the oriented CDS, or None."""
        s, e, strand = gene["start"], gene["end"], gene["strand"]
        if not (s <= pos <= e):
            return None
        off = (pos - s) if strand == "+" else (e - pos)
        off -= gene.get("phase", 0) or 0
        return off if off >= 0 else None

    annotations = []
    affected = OrderedDict()

    for r in rows:
        contig, pos = r["contig"], int(r["pos"])
        ref_a, alt_a, vtype = r["ref"], r["alt"], r["type"]
        rec = dict(r)
        rec.update(effect="", impact="", gene="", locus_tag="", product="",
                   strand="", codon_ref="", codon_alt="", aa_ref="", aa_alt="",
                   aa_pos="", blosum62="", distance_to_gene="",
                   transcript_effect="", protein_change="", cds_identity="",
                   reference_gene_status="")

        hits = find_overlapping(contig, pos)
        region = "intergenic"
        if hits:
            g = hits[0]
            region = "CDS"
            r_cds_seq = reference_cds(g)
            idx = cds_offset(g, pos)
            rec.update(gene=g["gene"], locus_tag=g["locus"],
                       product=g["product"], strand=g["strand"],
                       reference_gene_status=("pseudogene" if g.get("pseudo")
                                              else "protein_coding"))
            if g.get("pseudo"):
                # The reference CDS is itself disrupted (frameshift / internal
                # stop).  Translating it produces a meaningless protein, so the
                # reference-vs-query protein identity must not be interpreted as
                # a measure of divergence.
                rec["cds_identity_note"] = ("reference gene is annotated as a "
                                            "pseudogene; protein-level identity "
                                            "is not meaningful")
            if idx is not None:
                cidx = idx // 3
                frame_off = idx % 3
                rec["aa_pos"] = cidx + 1
                cod_start = cidx * 3
                if cod_start + 3 <= len(r_cds_seq):
                    ref_codon = r_cds_seq[cod_start:cod_start + 3]
                    rec["codon_ref"] = ref_codon
                    if vtype == "SNP":
                        alt_codon = list(ref_codon)
                        alt_codon[frame_off] = alt_a
                        alt_codon = "".join(alt_codon)
                        rec["codon_alt"] = alt_codon
                        ra, aa = translate(ref_codon), translate(alt_codon)
                        rec["aa_ref"], rec["aa_alt"] = ra, aa
                        rec["blosum62"] = blosum62(ra, aa)
                        # effect
                        if ra == aa:
                            eff = "synonymous_variant"
                        elif aa == "*":
                            eff = "stop_gained"
                        elif ra == "*":
                            eff = "stop_lost"
                        elif cidx == 0 and aa != "M":
                            eff = "start_lost"
                        else:
                            eff = "missense_variant"
                        # consequence of the first codon only; refined below
                        rec["effect"] = eff
                        rec["protein_change"] = f"p.{ra}{cidx+1}{aa}"
                    else:
                        rec["effect"] = ("frameshift_variant"
                                         if len(alt_a) % 3 else "inframe_variant")
            # ------- refine with the projected query allele -----------------
            if bam is not None:
                qseq, _cov = project_allele(bam, contig, g["start"] - 1,
                                            g["end"], g["strand"])
                rseq = reference_cds(g)
                rprot, qprot = translate(rseq), translate(qseq)
                # trim trailing stop for comparison
                if rprot.endswith("*"):
                    rprot_c = rprot[:-1]
                else:
                    rprot_c = rprot
                if qprot.endswith("*"):
                    qprot_c = qprot[:-1]
                else:
                    qprot_c = qprot
                if rprot_c and qprot_c and not g.get("pseudo"):
                    same = sum(1 for x, y in zip(rprot_c, qprot_c) if x == y)
                    ident = 100.0 * same / max(len(rprot_c), len(qprot_c))
                    rec["cds_identity"] = f"{ident:.1f}"
                    # first difference
                    fd = ""
                    for i, (x, y) in enumerate(zip(rprot_c, qprot_c)):
                        if x != y:
                            fd = f"p.{x}{i+1}{y}"
                            break
                    if not fd and len(rprot_c) != len(qprot_c):
                        fd = f"p.{rprot_c[min(len(rprot_c),len(qprot_c))-1]}" \
                             f"{min(len(rprot_c),len(qprot_c))}fs" \
                             if len(qprot_c) != len(rprot_c) else ""
                    rec["protein_change"] = fd or rec["protein_change"] or "p.="
                    # a premature stop inside the query protein => truncation
                    qstop = qprot_c.find("*")
                    if qstop != -1 and qstop < len(rprot_c) - 1:
                        rec["transcript_effect"] = f"premature_stop_at_aa{qstop+1}"
                        if rec["effect"] in ("missense_variant",
                                             "synonymous_variant", ""):
                            rec["effect"] = ("frameshift_variant"
                                             if vtype != "SNP"
                                             else "stop_gained")
                    elif len(qprot_c) != len(rprot_c):
                        rec["transcript_effect"] = (
                            f"protein_length_{len(rprot_c)}_to_{len(qprot_c)}")
            affected.setdefault(g["locus"] or g["gene"] or f"{contig}:{g['start']}",
                                dict(gene=g, variants=[]))
            affected[g["locus"] or g["gene"] or f"{contig}:{g['start']}"]["variants"].append(rec)
        elif noncds and any(f["contig"] == contig and f["start"] <= pos <= f["end"]
                            for f in noncds):
            f = [f for f in noncds
                 if f["contig"] == contig and f["start"] <= pos <= f["end"]][0]
            region = f["type"]
            rec["effect"] = f"{f['type']}_variant"
            rec.update(gene=f["gene"], locus_tag=f["locus"], product=f["product"])
        else:
            g, d, side = nearest_gene(contig, pos)
            if g is not None and d is not None and d <= a.updown:
                region = "intergenic"
                rec["effect"] = f"{side}_gene_variant"
                rec["distance_to_gene"] = d
                rec.update(gene=g["gene"], locus_tag=g["locus"], product=g["product"])
            else:
                rec["effect"] = "intergenic_variant"
                if g is not None:
                    rec["distance_to_gene"] = d
                    rec.update(gene=g["gene"], locus_tag=g["locus"],
                               product=g["product"])
        rec["region"] = region
        annotations.append(rec)

    # ---- impact ranking ---------------------------------------------------
    ORDER = {"stop_gained": 0, "frameshift_variant": 0, "start_lost": 1,
             "stop_lost": 1, "inframe_variant": 2, "missense_variant": 3,
             "synonymous_variant": 4, "tRNA_variant": 2, "rRNA_variant": 2,
             "upstream_gene_variant": 5, "downstream_gene_variant": 5,
             "intergenic_variant": 6}
    for rec in annotations:
        e = rec.get("effect", "")
        rec["impact"] = ("HIGH" if ORDER.get(e, 9) <= 1 else
                         "MODERATE" if ORDER.get(e, 9) <= 3 else
                         "LOW" if ORDER.get(e, 9) <= 5 else "MODIFIER")

    # ---- write variant-centric table -------------------------------------
    cols = ["contig", "pos", "ref", "alt", "type", "confidence",
            "n_independent_methods", "n_callers", "callers", "flags", "effect",
            "impact", "region", "locus_tag", "gene", "product",
            "reference_gene_status", "strand", "aa_pos", "codon_ref",
            "codon_alt", "aa_ref", "aa_alt", "blosum62", "protein_change",
            "transcript_effect", "cds_identity", "distance_to_gene",
            "obs_depth", "obs_alt_frac", "asm_mapq", "asm_depth", "sim_depth",
            "mummer_buff", "mummer_dist"]
    var_f = os.path.join(a.outdir, "variant_annotation.tsv")
    with open(var_f, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for rec in sorted(annotations, key=lambda x: (x["contig"], int(x["pos"]))):
            fh.write("\t".join(str(rec.get(c, "")) for c in cols) + "\n")

    # ---- write gene-centric table + FASTA --------------------------------
    gene_f = os.path.join(a.outdir, "gene_impact.tsv")
    gcols = ["locus_tag", "gene", "product", "contig", "start", "end", "strand",
             "reference_gene_status", "reference_note", "n_variants",
             "variant_types", "worst_effect", "worst_impact", "protein_change",
             "cds_identity", "transcript_effect", "ref_protein_len",
             "qry_protein_len", "ref_internal_stops", "qry_internal_stops",
             "query_orf_status", "variant_positions"]
    qaa_f = os.path.join(a.outdir, "affected_genes.faa")
    raa_f = os.path.join(a.outdir, "reference_genes.faa")
    nq = nr = 0
    with open(gene_f, "w") as fh, open(qaa_f, "w") as fq, open(raa_f, "w") as fr:
        fh.write("\t".join(gcols) + "\n")
        for key, info in affected.items():
            g = info["gene"]
            vs = info["variants"]
            worst = min(vs, key=lambda v: ORDER.get(v.get("effect", ""), 9))
            types = Counter(v["type"] for v in vs)
            prot_change = next((v["protein_change"] for v in vs
                                if v.get("protein_change")), "")
            ident = next((v["cds_identity"] for v in vs if v.get("cds_identity")), "")
            tex = next((v["transcript_effect"] for v in vs
                        if v.get("transcript_effect")), "")
            # --- ORF integrity of both alleles -----------------------------
            # The reference annotation marks some CDSs as pseudogenes with an
            # internal frameshift, so simply translating them yields a
            # meaningless protein.  Counting internal stop codons in each allele
            # makes it explicit whether the query ORF is intact, restored, or
            # further disrupted.
            ref_len = qry_len = ref_stops = qry_stops = ""
            orf_status = ""
            if bam is not None:
                try:
                    rseq_g = reference_cds(g)
                    qseq_g, _ = project_allele(bam, g["contig"], g["start"] - 1,
                                               g["end"], g["strand"])
                    rp = translate(rseq_g)
                    qp = translate(qseq_g)
                    rc_ = rp[:-1] if rp.endswith("*") else rp
                    qc_ = qp[:-1] if qp.endswith("*") else qp
                    ref_len, qry_len = len(rc_), len(qc_)
                    ref_stops = rc_.count("*")
                    qry_stops = qc_.count("*")
                    if ref_stops == 0 and qry_stops == 0:
                        orf_status = "intact_in_both"
                    elif ref_stops > 0 and qry_stops == 0:
                        orf_status = "query_ORF_restored"
                    elif ref_stops == 0 and qry_stops > 0:
                        orf_status = "query_ORF_disrupted"
                    else:
                        orf_status = "disrupted_in_both"
                except Exception:
                    pass

            fh.write("\t".join(str(x) for x in [
                g["locus"], g["gene"], g["product"], g["contig"], g["start"],
                g["end"], g["strand"],
                "pseudogene" if g.get("pseudo") else "protein_coding",
                g.get("note", ""),
                len(vs),
                ";".join(f"{k}:{v}" for k, v in types.items()),
                worst.get("effect", ""), worst.get("impact", ""),
                prot_change, ident, tex,
                ref_len, qry_len, ref_stops, qry_stops, orf_status,
                ";".join(f"{v['contig']}:{v['pos']}{v['ref']}>{v['alt']}"
                         for v in sorted(vs, key=lambda x: int(x["pos"])))]) + "\n")
            if bam is not None:
                try:
                    qseq, _ = project_allele(bam, g["contig"], g["start"] - 1,
                                             g["end"], g["strand"])
                    rseq = reference_cds(g)
                    if qseq.strip("N"):
                        fq.write(f">{g['locus']}|{g['gene']}|{g['product']}|"
                                 f"QUERY\n{translate(qseq)}\n"); nq += 1
                        fr.write(f">{g['locus']}|{g['gene']}|{g['product']}|"
                                 f"REFERENCE\n{translate(rseq)}\n"); nr += 1
                except Exception:
                    pass

    # ---- summary ----------------------------------------------------------
    summary = dict(
        total_variants=len(annotations),
        effects=dict(Counter(r.get("effect", "") for r in annotations)),
        impacts=dict(Counter(r.get("impact", "") for r in annotations)),
        affected_genes=len(affected),
        affected_gene_list=[dict(locus=i["gene"]["locus"],
                                 gene=i["gene"]["gene"],
                                 product=i["gene"]["product"],
                                 n_variants=len(i["variants"]))
                            for i in affected.values()],
        high_impact_genes=[dict(locus=i["gene"]["locus"],
                                gene=i["gene"]["gene"],
                                product=i["gene"]["product"],
                                effect=min(i["variants"],
                                           key=lambda v: ORDER.get(v.get("effect", ""), 9)
                                           ).get("effect", ""),
                                protein_change=next((v["protein_change"]
                                                     for v in i["variants"]
                                                     if v.get("protein_change")), ""))
                           for i in affected.values()
                           if min(i["variants"],
                                  key=lambda v: ORDER.get(v.get("effect", ""), 9)
                                  ).get("impact") in ("HIGH", "MODERATE")],
    )
    with open(os.path.join(a.outdir, "annotation_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    if bam is not None:
        bam.close()

    print(f"[annot] {len(annotations)} variants annotated; "
          f"{len(affected)} genes affected")
    print(f"[annot] effects: {dict(summary['effects'])}")
    print(f"[annot] impacts: {dict(summary['impacts'])}")
    print(f"[annot] wrote {var_f}, {gene_f}, {qaa_f}({nq} seqs)")


if __name__ == "__main__":
    main()
