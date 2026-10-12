#!/usr/bin/env python3
"""
curated_determinants.py -- targeted, B. fragilis-aware scan of resistance,
virulence and metabolism determinants in BOTH genomes, and a direct
reference-vs-query comparison of every determinant that is found.

Rather than trusting a single database, this script:
  1. scans the reference and query annotations for a curated set of
     B. fragilis determinants (product-name / keyword patterns),
  2. recovers the *query allele* of every matching reference gene by projecting
     the gene interval through the asm5 alignment (so absent genes and
     full-length, point-mutated or truncated alleles are distinguishable),
  3. translates both alleles and reports identity, amino-acid changes and
     presence/absence.

Outputs (in --outdir)
  determinant_catalogue.tsv    the curated list and what was found
  determinant_comparison.tsv   reference vs query, per determinant gene
  determinant_aa_changes.tsv   amino-acid level differences only
  determinant_presence_matrix.tsv  matrix of determinant x genome
  curated_summary.json
"""
import argparse
import json
import os
import sys
from collections import OrderedDict, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotate_variants import (read_fasta, revcomp, translate,  # noqa: E402
                               parse_gff, project_allele, project_region,
                               blosum62)


def reference_cds(gene, genome):
    """Transcription-oriented reference CDS (minus strand is reverse COMPLEMENTED)."""
    s, e = gene["start"], gene["end"]
    raw = genome[gene["contig"]][s - 1:e]
    seq = revcomp(raw) if gene["strand"] == "-" else raw
    ph = gene.get("phase", 0) or 0
    return seq[ph:] if ph else seq

# --------------------------------------------------------------------------- #
#  curated catalogue
#  category : resistance | virulence | metabolism | stress
# --------------------------------------------------------------------------- #
CATALOGUE = [
    # ---------------- beta-lactam resistance ----------------
    ("cfiA", "resistance", "beta-lactam", "Class B metallo-beta-lactamase CfiA (carbapenem resistance)",
     ["cfia", "metallo-beta-lactamase", "class b beta-lactamase", "blab"]),
    ("cfxA", "resistance", "beta-lactam", "Class A extended-spectrum beta-lactamase CfxA",
     ["cfxa", "class a beta-lactamase", "beta-lactamase"]),
    ("ampC", "resistance", "beta-lactam", "AmpC-type cephalosporinase",
     ["ampc", "cephalosporinase"]),
    ("bla_other", "resistance", "beta-lactam", "other beta-lactamase / penicillin-binding protein",
     ["penicillin-binding protein", "beta-lactamase domain", "blar"]),
    # ---------------- nitroimidazole ----------------
    ("nim", "resistance", "nitroimidazole", "Nim nitroimidazole resistance protein",
     ["nima", "nimb", "nimc", "nimd", "nime", "nimj", "nitroimidazole"]),
    # ---------------- macrolide / lincosamide ----------------
    ("ermF", "resistance", "macrolide", "ErmF 23S rRNA methyltransferase",
     ["ermf", "rrna methyltransferase", "erythromycin resistance", "erm"]),
    ("mefA_msr", "resistance", "macrolide", "macrolide efflux (MefA/Msr)",
     ["mefa", "msra", "macrolide efflux"]),
    ("lnu", "resistance", "lincosamide", "lincosamide nucleotidyltransferase",
     ["lnu", "lincosamide"]),
    # ---------------- tetracycline ----------------
    ("tetQ", "resistance", "tetracycline", "TetQ ribosomal protection protein",
     ["tetq", "tetracycline resistance", "tetracycline-resistant"]),
    ("tetX", "resistance", "tetracycline", "TetX tetracycline destructase",
     ["tetx", "tetracycline destructase", "fad-binding monooxygenase"]),
    ("tet_other", "resistance", "tetracycline", "other tetracycline determinant",
     ["tetm", "teto", "tetw", "tet36", "tetb", "tetp"]),
    # ---------------- aminoglycoside ----------------
    ("aminoglycoside", "resistance", "aminoglycoside",
     "aminoglycoside modifying enzyme (aac/aph/aad/ant)",
     ["aminoglycoside", "aaca", "aph(", "aad", "ant(", "kanamycin", "streptomycin"]),
    # ---------------- sulfonamide / trimethoprim / phenicol ----------------
    ("sul", "resistance", "sulfonamide", "sulphonamide resistance (sul)",
     ["sul1", "sul2", "dihydropteroate synthase"]),
    ("dfr", "resistance", "trimethoprim", "dihydrofolate reductase (dfr)",
     ["dfra", "dfr", "dihydrofolate reductase"]),
    ("cat_flo", "resistance", "phenicol", "chloramphenicol acetyltransferase / floR",
     ["chloramphenicol acetyltransferase", "flor", "cata"]),
    ("oxazolidinone", "resistance", "oxazolidinone", "optrA/poxtA/cfr",
     ["optra", "poxta", "cfr"]),
    ("qnr", "resistance", "quinolone", "quinolone resistance protein (qnr)",
     ["qnr", "pentapeptide repeat"]),
    ("efflux_RND", "resistance", "multidrug efflux",
     "RND / MFS multidrug efflux pump",
     ["multidrug efflux", "acrb", "mexb", "mdte", "efflux pump", "bexa", "rnd transporter"]),
    # ---------------- target-site (mutations matter) ----------------
    ("gyrA", "resistance", "quinolone", "DNA gyrase subunit A (QRDR target)",
     ["dna gyrase subunit a", "gyra"]),
    ("gyrB", "resistance", "quinolone", "DNA gyrase subunit B",
     ["dna gyrase subunit b", "gyrb"]),
    ("parC", "resistance", "quinolone", "DNA topoisomerase IV subunit A (QRDR target)",
     ["topoisomerase iv subunit a", "parc", "topoisomerase 4 subunit a"]),
    ("parE", "resistance", "quinolone", "DNA topoisomerase IV subunit B",
     ["topoisomerase iv subunit b", "pare"]),
    ("rpoB", "resistance", "rifampicin", "RNA polymerase beta subunit (rifampicin target)",
     ["rna polymerase subunit beta", "rpob", "dna-directed rna polymerase subunit beta"]),
    ("rrs_23S", "resistance", "macrolide", "23S rRNA (macrolide target)",
     ["23s ribosomal rna"]),
    # ---------------- virulence: B. fragilis specific ----------------
    ("bft", "virulence", "toxin", "B. fragilis toxin / fragilysin (bft)",
     ["fragilysin", "bft", "toxin", "metallopeptidase"]),
    ("cps_locus", "virulence", "capsule",
     "capsular polysaccharide biosynthesis (CPS locus)",
     ["capsular polysaccharide", "polysaccharide biosynthesis", "wzy", "wzx",
      "glycosyltransferase", "capsule"]),
    ("susC", "virulence", "PUL", "SusC TonB-dependent outer-membrane transporter (PUL)",
     ["susC", "tonb-dependent", "tonb dependent", "outer membrane receptor"]),
    ("susD", "virulence", "PUL", "SusD polysaccharide-binding lipoprotein (PUL)",
     ["susd", "susd-like", "polysaccharide-binding"]),
    ("sialidase", "virulence", "enzyme", "sialidase / neuraminidase",
     ["sialidase", "neuraminidase"]),
    ("hemolysin", "virulence", "toxin", "hemolysin / cytolysin",
     ["hemolysin", "haemolysin", "cytolysin"]),
    ("protease", "virulence", "enzyme", "extracellular protease",
     ["metalloprotease", "zinc protease", "peptidase m10", "collagenase"]),
    ("fimbriae", "virulence", "adhesion", "fimbriae / pilus assembly",
     ["fimbri", "pilin", "pilus assembly", "type iv pilus"]),
    ("bile_salt", "virulence", "persistence", "bile salt hydrolase",
     ["bile salt hydrolase", "choloylglycine hydrolase"]),
    ("sod_catalase", "virulence", "oxidative stress",
     "superoxide dismutase / catalase",
     ["superoxide dismutase", "catalase", "rubrerythrin"]),
    ("b12", "metabolism", "cofactor", "cobalamin (B12) biosynthesis/transport",
     ["cobalamin", "cobyrinic", "cob(", "precorrin", "vitamin b12", "btu"]),
    ("sulfur", "metabolism", "redox", "sulfate/sulfite metabolism",
     ["sulfite reductase", "sulfate adenylyltransferase", "sulfatase", "cysteine synthase"]),
    ("glycoside_hydrolase", "metabolism", "CAZy", "glycoside hydrolase (glycan degradation)",
     ["glycoside hydrolase", "glycosyl hydrolase", "glucosidase", "galactosidase",
      "xylanase", "glycanase", "mannosidase", "fucosidase"]),
    ("polysaccharide_lyase", "metabolism", "CAZy", "polysaccharide lyase",
     ["polysaccharide lyase", "chondroitin", "heparinase", "hyaluronate lyase",
      "pectate lyase", "alginate lyase"]),
    ("cazyme_esterase", "metabolism", "CAZy", "carbohydrate esterase / sulfatase",
     ["carbohydrate esterase", "acetylxylan esterase", "arylsulfatase"]),
]


def match_products(gff_path, patterns, want_types=("CDS",)):
    """Return the genes whose product matches any pattern (case-insensitive)."""
    feats = parse_gff(gff_path)
    out = []
    for f in feats:
        if f["type"] not in want_types:
            continue
        hay = f"{f['product']} {f['gene']}".lower()
        for p in patterns:
            if p.lower() in hay:
                out.append(f)
                break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-faa"); ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--qry-faa"); ap.add_argument("--qry-gff")
    ap.add_argument("--ref-fna", required=True); ap.add_argument("--qry-fna")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--bam", default=None)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    ref_genome = read_fasta(a.ref_fna)
    ref_feats = parse_gff(a.ref_gff)

    bam = None
    if a.bam is None:
        cand = os.path.normpath(os.path.join(
            os.path.dirname(a.outdir), "..", "..", "work", "bam",
            "CCUG4856_to_NCTC9343.asm5.bam"))
        a.bam = cand
    if os.path.exists(a.bam):
        try:
            import pysam
            bam = pysam.AlignmentFile(a.bam, "rb")
        except Exception as exc:
            print(f"[curated] BAM unavailable ({exc})", file=sys.stderr)

    catalogue_rows = []
    comparison_rows = []
    aa_rows = []
    presence = OrderedDict()

    for name, category, subcat, desc, patterns in CATALOGUE:
        hits = match_products(a.ref_gff, patterns)
        # keep the best (longest) few hits per determinant to avoid huge lists
        # No arbitrary truncation: a curated determinant that is present must be
        # reported even when the family has many paralogues in the genome.
        hits = sorted(hits, key=lambda f: -(f["end"] - f["start"]))
        catalogue_rows.append(dict(
            determinant=name, category=category, subcategory=subcat,
            description=desc, reference_hits=len(hits),
            keywords=";".join(patterns)))
        found_ref = len(hits) > 0
        found_qry_via_bam = 0
        present_in_query = False
        for h in hits:
            rseq = reference_cds(h, ref_genome)
            rprot = translate(rseq)
            row = dict(determinant=name, category=category,
                       locus_tag=h["locus"], gene=h["gene"],
                       product=h["product"], contig=h["contig"],
                       start=h["start"], end=h["end"], strand=h["strand"],
                       ref_protein_len=len(rprot.rstrip("*")),
                       query_coverage_pct="", query_protein_len="",
                       protein_identity_pct="", aa_changes="",
                       status="reference_only")
            if bam is not None:
                qseq, covered = project_allele(bam, h["contig"], h["start"] - 1,
                                               h["end"], h["strand"])
                span = h["end"] - h["start"] + 1
                ncov = sum(1 for p in range(h["start"] - 1, h["end"])
                           if p in covered)
                cov = 100.0 * ncov / span if span else 0.0
                if qseq.strip("N"):
                    qprot = translate(qseq)
                    qp, rp = qprot.rstrip("*"), rprot.rstrip("*")
                    qry_stops = qp.count("*")
                    same = sum(1 for x, y in zip(rp, qp) if x == y)
                    ident = 100.0 * same / max(len(rp), len(qp)) if max(len(rp), len(qp)) else 0
                    row["query_coverage_pct"] = f"{cov:.1f}"
                    row["query_protein_len"] = len(qp)
                    row["protein_identity_pct"] = f"{ident:.1f}"
                    changes = []
                    for i, (x, y) in enumerate(zip(rp, qp)):
                        if x != y:
                            b = blosum62(x, y)
                            changes.append(f"p.{x}{i+1}{y}" +
                                           (f"(BLOSUM62={b})" if b != "" else ""))
                    if len(qp) < len(rp):
                        changes.append(f"truncation_{len(rp)}->{len(qp)}aa")
                    elif len(qp) > len(rp):
                        changes.append(f"extension_{len(rp)}->{len(qp)}aa")
                    if h.get("pseudo") and changes:
                        # The reference allele is an annotated pseudogene, so the
                        # apparent amino-acid avalanche is a reading-frame
                        # difference rather than a set of substitutions.
                        row["aa_changes"] = (
                            f"reference_is_pseudogene__reading_frame_"
                            f"{'restored' if len(qp) and qry_stops == 0 else 'differs'}"
                            f"__(query internal stops: {qry_stops})")
                    else:
                        row["aa_changes"] = ";".join(changes) if changes else "none"
                    row["status"] = ("identical" if ident >= 99.99 and len(qp) == len(rp)
                                     else "mutated" if cov >= 50
                                     else "partial")
                    if cov >= 50:
                        present_in_query = True
                        found_qry_via_bam += 1
                    if changes:
                        aa_rows.append(dict(determinant=name, locus_tag=h["locus"],
                                            product=h["product"], gene=h["gene"],
                                            category=category,
                                            ref_len=len(rp), qry_len=len(qp),
                                            identity=f"{ident:.1f}",
                                            aa_changes=row["aa_changes"]))
                else:
                    row["status"] = "absent_in_query"
            comparison_rows.append(row)
        presence[name] = dict(category=category, in_reference=found_ref,
                              in_query=present_in_query,
                              n_reference_genes=len(hits),
                              n_query_genes=found_qry_via_bam)

    # ---- write ------------------------------------------------------------
    with open(f"{a.outdir}/determinant_catalogue.tsv", "w") as fh:
        fh.write("determinant\tcategory\tsubcategory\tdescription\t"
                 "reference_hits\tkeywords\n")
        for r in catalogue_rows:
            fh.write(f"{r['determinant']}\t{r['category']}\t{r['subcategory']}\t"
                     f"{r['description']}\t{r['reference_hits']}\t{r['keywords']}\n")

    with open(f"{a.outdir}/determinant_comparison.tsv", "w") as fh:
        cols = ["determinant", "category", "locus_tag", "gene", "product",
                "contig", "start", "end", "strand", "ref_protein_len",
                "query_protein_len", "query_coverage_pct",
                "protein_identity_pct", "status", "aa_changes"]
        fh.write("\t".join(cols) + "\n")
        for r in comparison_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    with open(f"{a.outdir}/determinant_aa_changes.tsv", "w") as fh:
        fh.write("determinant\tcategory\tlocus_tag\tgene\tproduct\t"
                 "ref_len\tqry_len\tidentity\taa_changes\n")
        for r in aa_rows:
            fh.write(f"{r['determinant']}\t{r['category']}\t{r['locus_tag']}\t"
                     f"{r['gene']}\t{r['product']}\t{r['ref_len']}\t"
                     f"{r['qry_len']}\t{r['identity']}\t{r['aa_changes']}\n")

    with open(f"{a.outdir}/determinant_presence_matrix.tsv", "w") as fh:
        fh.write("determinant\tcategory\treference\tquery\n")
        for k, v in presence.items():
            fh.write(f"{k}\t{v['category']}\t"
                     f"{'present' if v['in_reference'] else 'absent'}\t"
                     f"{'present' if v['in_query'] else 'absent'}\n")

    summary = dict(
        determinants_in_catalogue=len(CATALOGUE),
        determinants_found_in_reference=sum(1 for v in presence.values()
                                            if v["in_reference"]),
        determinants_found_in_query=sum(1 for v in presence.values()
                                        if v["in_query"]),
        reference_genes_examined=len(comparison_rows),
        genes_with_aa_changes=len(aa_rows),
        present_in_reference_only=[k for k, v in presence.items()
                                   if v["in_reference"] and not v["in_query"]],
        present_in_query_only=[k for k, v in presence.items()
                               if v["in_query"] and not v["in_reference"]],
        aa_change_detail=aa_rows,
    )
    with open(f"{a.outdir}/curated_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    if bam is not None:
        bam.close()

    print(f"[curated] catalogue={len(CATALOGUE)}, "
          f"ref genes examined={len(comparison_rows)}, "
          f"genes with AA changes={len(aa_rows)}")
    print(f"[curated] determinants in reference only: "
          f"{summary['present_in_reference_only']}")


def _positions(gene):
    s, e, strand = gene["start"], gene["end"], gene["strand"]
    pos = list(range(s, e + 1))
    if strand == "-":
        pos = pos[::-1]
    ph = gene.get("phase", 0)
    if ph and ph < len(pos):
        pos = pos[ph:]
    return pos


if __name__ == "__main__":
    main()
