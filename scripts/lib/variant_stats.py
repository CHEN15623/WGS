#!/usr/bin/env python3
"""
variant_stats.py -- descriptive statistics for the reconciled variant set.

Produces (in --outdir):
  variant_type_counts.tsv        SNP / INS / DEL / MNV totals
  snp_substitution_spectrum.tsv  12-channel and 6-class substitution spectrum
  ts_tv.tsv                      transition/transversion ratio
  per_contig_variants.tsv        counts and density per replicon
  variant_context.tsv            genic vs intergenic, strand, gene overlap
  window_density.tsv             10 kb windowed density (for plotting)
  indel_length_distribution.tsv
  variant_stats_summary.json
"""
import argparse
import json
import os
from collections import OrderedDict, Counter, defaultdict

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
PURINES = set("AG")
PYRIMIDINES = set("CT")


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


def read_gff_genes(path):
    """Return {contig: [(start0, end, strand, locus, gene, product)]} for CDS."""
    genes = defaultdict(list)
    locus = gene = product = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] != "CDS":
                continue
            attr = f[8]
            d = {}
            for kv in attr.split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    d[k] = v
            locus = d.get("locus_tag", d.get("ID", ""))
            gene = d.get("gene", d.get("Name", ""))
            product = d.get("product", d.get("Parent", ""))
            genes[f[0]].append((int(f[3]) - 1, int(f[4]), f[6], locus, gene, product))
    for c in genes:
        genes[c].sort()
    return genes


def find_gene(genes, contig, pos1):
    """pos1 is 1-based; return the CDS dict overlapping it (or None)."""
    pos0 = pos1 - 1
    for g in genes.get(contig, []):
        if g[0] <= pos0 < g[1]:
            return g
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged", required=True)
    ap.add_argument("--confident", default=None,
                    help="verified confident variant table (optional); when given a "
                         "confident-only statistics block is also written")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--window", type=int, default=10000)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    ref = read_fasta(a.ref)
    genes = read_gff_genes(a.gff)

    rows = []
    with open(a.merged) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                f += [""] * (len(hdr) - len(f))
            rows.append(dict(zip(hdr, f)))

    # ---- type counts ------------------------------------------------------
    types = Counter(r["type"] for r in rows)
    with open(f"{a.outdir}/variant_type_counts.tsv", "w") as fh:
        fh.write("variant_type\tcount\tpercent\n")
        for t, n in types.most_common():
            fh.write(f"{t}\t{n}\t{100.0*n/len(rows):.2f}\n" if rows else f"{t}\t{n}\t0\n")

    # ---- substitution spectrum + ts/tv -----------------------------------
    spec12 = Counter()
    spec6 = Counter()
    ts = tv = 0
    for r in rows:
        if r["type"] != "SNP":
            continue
        rb, ab = r["ref"].upper(), r["alt"].upper()
        if rb not in "ACGT" or ab not in "ACGT":
            continue
        spec12[f"{rb}>{ab}"] += 1
        pair = (rb, ab)
        if pair in TRANSITIONS:
            ts += 1
            spec6["transition"] += 1
            if rb in PURINES:
                spec6["purine_transition"] += 1
            else:
                spec6["pyrimidine_transition"] += 1
        else:
            tv += 1
            spec6["transversion"] += 1
            spec6["purine_transversion" if rb in PURINES
                  else "pyrimidine_transversion"] += 1

    with open(f"{a.outdir}/snp_substitution_spectrum.tsv", "w") as fh:
        fh.write("substitution\tcount\n")
        order = ["A>C", "A>G", "A>T", "C>A", "C>G", "C>T",
                 "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"]
        for k in order:
            fh.write(f"{k}\t{spec12.get(k, 0)}\n")
        fh.write("\n")
        fh.write("class\tcount\n")
        for k, v in spec6.most_common():
            fh.write(f"{k}\t{v}\n")

    tstv = (ts / tv) if tv else float("nan")
    with open(f"{a.outdir}/ts_tv.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"transitions\t{ts}\n")
        fh.write(f"transversions\t{tv}\n")
        fh.write(f"ts_tv_ratio\t{tstv:.3f}\n")

    # ---- per-contig -------------------------------------------------------
    count_by_contig = Counter(r["contig"] for r in rows)
    with open(f"{a.outdir}/per_contig_variants.tsv", "w") as fh:
        fh.write("contig\tlength_bp\tn_variants\tn_SNP\tn_INS\tn_DEL\t"
                 "variants_per_kb\n")
        for c, seq in ref.items():
            n = count_by_contig.get(c, 0)
            sub = Counter(r["type"] for r in rows if r["contig"] == c)
            fh.write(f"{c}\t{len(seq)}\t{n}\t{sub.get('SNP',0)}\t"
                     f"{sub.get('INS',0)}\t{sub.get('DEL',0)}\t"
                     f"{1000.0*n/len(seq):.4f}\n")

    # ---- context: genic vs intergenic, strand, homopolymer ----------------
    ctx = Counter()
    with open(f"{a.outdir}/variant_context.tsv", "w") as fh:
        fh.write("contig\tpos\tref\talt\ttype\tconfidence\tregion\tstrand\t"
                 "locus_tag\tgene\tproduct\thomopolymer_len\n")
        for r in rows:
            c, p = r["contig"], int(r["pos"])
            g = find_gene(genes, c, p)
            if g:
                region = "CDS"
                strand, locus, gene, product = g[2], g[3], g[4], g[5]
            else:
                region, strand, locus, gene, product = "intergenic", ".", "", "", ""
            ctx[region] += 1
            if strand == "+":
                ctx["CDS_forward"] += 1
            elif strand == "-":
                ctx["CDS_reverse"] += 1
            # homopolymer length immediately 5' of an indel (assembly-artefact proxy)
            homo = 0
            seq = ref.get(c, "")
            if r["type"] in ("INS", "DEL") and seq:
                base = seq[p - 2] if p - 2 >= 0 else ""
                if base:
                    i = p - 2
                    while i >= 0 and seq[i] == base:
                        homo += 1; i -= 1
                    j = p - 1
                    while j < len(seq) and seq[j] == base:
                        homo += 1; j += 1
            fh.write(f"{c}\t{p}\t{r['ref']}\t{r['alt']}\t{r['type']}\t"
                     f"{r['confidence']}\t{region}\t{strand}\t{locus}\t{gene}\t"
                     f"{product}\t{homo}\n")

    # ---- windowed density -------------------------------------------------
    win = a.window
    with open(f"{a.outdir}/window_density.tsv", "w") as fh:
        fh.write("contig\tstart\tend\tmid\tn_variants\tn_SNP\tn_INDEL\n")
        for c, seq in ref.items():
            nwin = (len(seq) + win - 1) // win
            vc = Counter(); sc = Counter(); ic = Counter()
            for r in rows:
                if r["contig"] != c:
                    continue
                w = (int(r["pos"]) - 1) // win
                vc[w] += 1
                if r["type"] == "SNP":
                    sc[w] += 1
                else:
                    ic[w] += 1
            for w in range(nwin):
                s, e = w * win + 1, min((w + 1) * win, len(seq))
                fh.write(f"{c}\t{s}\t{e}\t{(s+e)//2}\t{vc.get(w,0)}\t"
                         f"{sc.get(w,0)}\t{ic.get(w,0)}\n")

    # ---- indel length distribution ---------------------------------------
    lens = Counter()
    for r in rows:
        if r["type"] == "INS":
            lens[len(r["alt"])] += 1
        elif r["type"] == "DEL":
            lens[-(len(r["ref"]))] += 1
    with open(f"{a.outdir}/indel_length_distribution.tsv", "w") as fh:
        fh.write("indel_length_bp\tcount\n")
        for k in sorted(lens):
            fh.write(f"{k}\t{lens[k]}\n")

    # ---- summary ----------------------------------------------------------
    summary = dict(
        total_variants=len(rows),
        types=dict(types),
        transitions=ts, transversions=tv,
        ts_tv_ratio=None if tv == 0 else round(ts / tv, 3),
        genic=ctx.get("CDS", 0), intergenic=ctx.get("intergenic", 0),
        cds_forward=ctx.get("CDS_forward", 0),
        cds_reverse=ctx.get("CDS_reverse", 0),
        per_contig=dict(count_by_contig),
        confidence=dict(Counter(r["confidence"] for r in rows)),
        window_size=win,
    )
    with open(f"{a.outdir}/variant_stats_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    # ---- confident-only statistics ---------------------------------------
    if a.confident and os.path.exists(a.confident):
        crows = []
        with open(a.confident) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                f += [""] * (len(hdr) - len(f))
                crows.append(dict(zip(hdr, f)))
        ctypes = Counter(r["type"] for r in crows)
        with open(f"{a.outdir}/confident_variant_stats.tsv", "w") as fh:
            fh.write("metric\tvalue\n")
            fh.write(f"confident_variant_events\t{len(crows)}\n")
            fh.write(f"confident_SNPs\t{ctypes.get('SNP',0)}\n")
            fh.write(f"confident_insertions\t{ctypes.get('INS',0)}\n")
            fh.write(f"confident_deletions\t{ctypes.get('DEL',0)}\n")
            fh.write(f"confident_MNV\t{ctypes.get('MNV',0)}\n")
            per = Counter(r["contig"] for r in crows)
            for k, v in per.items():
                fh.write(f"confident_on_{k}\t{v}\n")
        with open(f"{a.outdir}/confident_events.tsv", "w") as fh:
            cols = [c for c in ["contig", "pos", "ref", "alt", "type",
                                "n_independent_methods", "callers",
                                "obs_depth", "obs_alt_frac", "event_span"]
                    if c in (crows[0] if crows else {}) or True]
            fh.write("\t".join(cols) + "\n")
            for r in crows:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
        print(f"[vstats] confident events: {len(crows)} {dict(ctypes)}")

    print(f"[vstats] {len(rows)} variants | types={dict(types)} | "
          f"ts/tv={summary['ts_tv_ratio']} | genic={summary['genic']} "
          f"intergenic={summary['intergenic']}")


if __name__ == "__main__":
    main()
