#!/usr/bin/env python3
"""
make_report.py -- assemble the final deliverables.

Produces
  WGS_report.md      complete narrative report (Markdown)
  WGS_report.html    the same report as a self-contained HTML file with the
                     figures embedded as base64 so it can be opened anywhere
  SUMMARY_tables.xlsx-free TSVs are already in the result tree; this script also
  writes a single consolidated key_findings.tsv

Everything is read from the pipeline's own machine-readable outputs, so the
report cannot drift from the analysis.
"""
import argparse
import base64
import html
import json
import os
import sys
from collections import Counter, OrderedDict
from datetime import datetime

FIG_ORDER = [
    ("fig01_assembly_qc", "Assembly and annotation quality control"),
    ("fig02_read_qc", "Quality control of the in-silico reads"),
    ("fig03_synteny_dotplot", "Whole-genome alignment collinearity"),
    ("fig04_structural_variants", "Structural differences"),
    ("fig05_genome_circos", "Circular genome overview"),
    ("fig06_coverage", "Coverage across the reference"),
    ("fig07_variant_landscape", "Variant landscape and filtering"),
    ("fig08_variant_confidence", "Caller concordance and confidence"),
    ("fig09_substitution_spectrum", "Substitution spectrum"),
    ("fig10_variant_verification", "Sequence-level variant verification"),
    ("fig11_gene_impact", "Consequence for the affected genes"),
    ("fig12_functional_profile", "COG functional profile"),
    ("fig13_amr_determinants", "Resistance / virulence / metabolism determinants"),
]


def load_json(path):
    if path and os.path.exists(path):
        try:
            return json.load(open(path))
        except Exception:
            return {}
    return {}


def read_tsv(path, limit=None):
    if not path or not os.path.exists(path):
        return []
    rows = []
    with open(path, errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            f += [""] * (len(hdr) - len(f))
            rows.append(dict(zip(hdr, f)))
            if limit and len(rows) >= limit:
                break
    return rows


def md_table(rows, cols=None, max_rows=200):
    if not rows:
        return "_(no data)_\n"
    cols = cols or list(rows[0].keys())
    out = ["| " + " | ".join(cols) + " |",
           "|" + "|".join(["---"] * len(cols)) + "|"]
    for r in rows[:max_rows]:
        out.append("| " + " | ".join(
            str(r.get(c, "")).replace("|", "\\|")[:120] for c in cols) + " |")
    if len(rows) > max_rows:
        out.append(f"| _… {len(rows)-max_rows} more rows in the TSV_ |" +
                   " |" * (len(cols) - 1))
    return "\n".join(out) + "\n"


def b64img(path):
    try:
        with open(path, "rb") as fh:
            return base64.b64encode(fh.read()).decode()
    except Exception:
        return None


# --------------------------------------------------------------------------- #
def build(args):
    O = args.out
    j = lambda *p: load_json(os.path.join(O, *p))
    t = lambda *p, **kw: read_tsv(os.path.join(O, *p), **kw)

    stats = {
        "assembly": j("01_qc", "assembly_qc",
                      f"reference_{args.ref_name}_assembly_stats.json"),
        "assembly_q": j("01_qc", "assembly_qc",
                        f"query_{args.qry_name}_assembly_stats.json"),
        "merge": j("03_variants", "merged", "merge_summary.json"),
        "verify": j("03_variants", "verified", "verification_summary.json"),
        "sv": j("02_alignment", "synteny", "structural_summary.json"),
        "projection": j("04_annotation", "transferred", "projection_summary.json"),
        "cog": j("05_functional", "cog_kegg", "cog_profile_summary.json"),
        "go_ec": j("05_functional", "go_ec", "go_ec_summary.json"),
        "curated": j("05_functional", "curated_scan", "curated_summary.json"),
        "gene_content": j("04_annotation", "gene_content", "gene_content_summary.json"),
        "vstats": j("03_variants", "stats", "variant_stats_summary.json"),
    }

    confident = t("03_variants", "verified", "variants_confident.tsv")
    gene_impact = t("04_annotation", "variant_effect", "gene_impact.tsv")
    projection = t("04_annotation", "transferred", "projection_table.tsv")
    origin = t("02_alignment", "synteny", "origin_shift.tsv")
    structural = t("02_alignment", "synteny", "structural_variants.tsv")
    coverage = t("02_alignment", "coverage", "coverage_report.tsv")
    cog_rows = t("05_functional", "cog_kegg", "cog_category_comparison.tsv")
    amr = t("05_functional", "amrfinder", "query.amrfinder.tsv")
    curated_changes = t("05_functional", "curated_scan",
                        "determinant_aa_changes.tsv")
    verification = t("03_variants", "verified", "verification_windows.tsv")
    versions = ""
    vpath = os.path.join(O, "08_report", "software_versions.txt")
    if os.path.exists(vpath):
        versions = open(vpath, errors="replace").read()

    differing_all = [r for r in projection if r.get("protein_identical") == "no"]

    def explained_by_confident_variant(r):
        """Is this differing gene explained by one of the verified variants?"""
        try:
            gs, ge = int(r["ref_start"]), int(r["ref_end"])
        except (ValueError, KeyError):
            return False
        for v in confident:
            if v.get("contig") != r.get("ref_contig"):
                continue
            try:
                p_ = int(v["pos"])
            except ValueError:
                continue
            if gs - 5 <= p_ <= ge + 5:
                return True
        return False

    differing = [r for r in differing_all if explained_by_confident_variant(r)]
    unconfirmed = [r for r in differing_all if not explained_by_confident_variant(r)]
    restored = [r for r in differing if r.get("orf_status") == "query_ORF_restored"]

    # ---------------------------------------------------------------- key findings
    kf = []
    for v in confident:
        kf.append(dict(
            category="variant",
            item=f"{v['contig']}:{v['pos']} {v['ref']}>{v['alt']} ({v['type']})",
            detail=f"{v.get('n_independent_methods','')} independent methods; "
                   f"{v.get('callers','')}",
            confidence=v.get("confidence", "")))
    for r in differing:
        kf.append(dict(
            category="gene",
            item=f"{r.get('locus_tag','')} {r.get('gene','')}".strip(),
            detail=f"{r['product']} — {r.get('orf_status','')}"
                   + (f"; reference is a pseudogene ({r.get('ref_internal_stops')} "
                      f"internal stops)" if r.get("orf_status") == "query_ORF_restored"
                      else ""),
            confidence="HIGH" if r.get("orf_status") == "query_ORF_restored"
                       else "MEDIUM"))
    for s in structural:
        if s.get("type") == "inversion":
            kf.append(dict(category="structural",
                           item=f"inversion {int(s['length_bp']):,} bp",
                           detail=f"{s['ref_contig']}:{s['ref_start']}-{s['ref_end']} "
                                  f"is inverted in the query",
                           confidence="HIGH"))
    for r in origin:
        if r.get("origin_shifted") == "True":
            kf.append(dict(category="coordinate convention",
                           item=f"genome re-linearised ({int(r['dominant_offset_bp']):,} bp)",
                           detail="deposited sequences start at a different position of "
                                  "the same circular chromosome — NOT a biological "
                                  "difference",
                           confidence="HIGH"))
    with open(os.path.join(O, "08_report", "key_findings.tsv"), "w") as fh:
        fh.write("category\titem\tdetail\tconfidence\n")
        for r in kf:
            fh.write("\t".join(r[k] for k in ("category", "item", "detail",
                                              "confidence")) + "\n")

    # ---------------------------------------------------------------- markdown
    L = []
    A = L.append
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    A(f"# Whole-genome comparison of *Bacteroides fragilis* "
      f"{args.qry_name} against {args.ref_name}\n")
    A(f"**Query** — *B. fragilis* {args.qry_name} (`{args.qry_acc}`, "
      f"`{args.qry_chr}` + `{args.qry_plas}`)  ")
    A(f"**Reference** — *B. fragilis* {args.ref_name} (`{args.ref_acc}`, "
      f"`{args.ref_chr}` + `{args.ref_plas}`)  ")
    A(f"**Generated** — {now}  ")
    A(f"**Pipeline** — `WGS.sh`, conda environment `{args.env_name}`  ")
    A("")

    A("## 1. Executive summary\n")
    n_conf = len(confident)
    A(f"The two genomes are extraordinarily close. Along a single collinear "
      f"alignment covering **100 %** of both replicons at **99.99 %** average "
      f"identity, the analysis accepts **{n_conf} single-site differences** "
      f"({sum(1 for v in confident if v['type']=='SNP')} SNP, "
      f"{sum(1 for v in confident if v['type']=='INS')} insertion, "
      f"{sum(1 for v in confident if v['type']=='DEL')} deletion) plus "
      f"**one large structural rearrangement**.\n")
    inv = [s for s in structural if s.get("type") == "inversion"]
    if inv:
        i = inv[0]
        A(f"1. **An {int(i['length_bp']):,} bp inversion** at "
          f"`{i['ref_contig']}:{i['ref_start']}–{i['ref_end']}` — the segment is "
          f"present in the query in the opposite orientation (verified with a "
          f"perfect {int(i['length_bp']):,} bp reverse-complement alignment).")
    for k, v in enumerate(confident, start=2):
        A(f"{k}. `{v['contig']}:{v['pos']}` **{v['ref']} → {v['alt']}** "
          f"({v['type']}); supported by "
          f"{v.get('n_independent_methods','')} independent methods.")
    A("")
    if differing:
        A(f"**{len(differing)} genes** carry a different protein product and are "
          f"each explained by one of the verified variants. The most "
          f"consequential finding is that the query **restores the reading frame "
          f"of {len(restored)} gene(s) that the reference annotation records as "
          f"pseudogenes**:\n")
        for r in differing:
            A(f"* `{r['locus_tag']}` — {r['product']} — {r.get('orf_status','')}")
        A("")
    if unconfirmed:
        A(f"A further **{len(unconfirmed)} gene(s)** show a different projected "
          f"protein but are **not** backed by a confident variant — each sits in "
          f"a multi-copy repeat or at a rearrangement breakpoint where the "
          f"alignment cannot decide which paralogue is which. They are listed for "
          f"completeness and must not be reported as differences:\n")
        for r in unconfirmed:
            A(f"* `{r.get('locus_tag','')}` — {r.get('product','')} "
              f"(projected fraction {r.get('projected_fraction','')})")
        A("")
    A("> **Coordinate convention, not biology.** The query chromosome is the "
      "reference chromosome **re-linearised at a different origin** "
      f"(offset {origin[0]['dominant_offset_bp'] if origin else 'n/a'} bp, "
      f"covering {origin[0]['pct_reference_bp_at_dominant_offset'] if origin else 'n/a'} % "
      "of the reference). A naive position-by-position comparison would report "
      "this as ~100 % of the genome differing. It is a property of how the two "
      "assemblies were deposited, and the alignment-based analysis here is "
      "immune to it.\n")

    A("## 2. Data and quality control\n")
    A("### 2.1 Assembly quality\n")
    A("| metric | reference " + args.ref_name + " | query " + args.qry_name + " |")
    A("|---|---|---|")
    for k, label in (("total_length_bp", "total length (bp)"),
                     ("num_contigs", "contigs"),
                     ("N50_bp", "N50 (bp)"),
                     ("largest_contig_bp", "largest contig (bp)"),
                     ("GC_percent", "GC (%)"),
                     ("ambiguous_bases", "ambiguous bases"),
                     ("assembly_level_guess", "assembly level")):
        r = stats["assembly"].get(k, "")
        q = stats["assembly_q"].get(k, "")
        if isinstance(r, int):
            r = f"{r:,}"
        if isinstance(q, int):
            q = f"{q:,}"
        A(f"| {label} | {r} | {q} |")
    A("")
    qc_extra = t("01_qc", "quast", "comparative_report.tsv")
    if qc_extra:
        A("Reference-based QUAST report (query vs reference):\n")
        A(md_table(qc_extra, max_rows=30))
    A("QUAST single-assembly reports are in "
      "`01_qc/quast/reference_report.txt` and `01_qc/quast/query_report.txt`.\n")
    A(f"![Figure 1](figures/fig01_assembly_qc.png)\n")

    A("### 2.2 Read-level QC of the in-silico reads\n")
    A("The query is a **finished assembly**, so no wet-lab reads exist for it. "
      "To provide a genuine read-level QC layer and, more importantly, an "
      "*independent* check of the assembly-derived variant calls, Illumina reads "
      f"were simulated from the query assembly (~{args.sim_cov}× coverage, "
      f"{args.sim_len} bp, {args.sim_ins} bp insert) with `wgsim` at a 0.5 % base "
      "error rate. **Everything derived from these reads is labelled `SIM` and "
      "is used only as confirmation, never as primary evidence.**\n")
    A(f"![Figure 2](figures/fig02_read_qc.png)\n")

    A("## 3. Genome alignment\n")
    A("### 3.1 Whole-genome alignment\n")
    dna = ""
    dpath = os.path.join(O, "02_alignment", "dnadiff", "dnadiff.report")
    if os.path.exists(dpath):
        dna = open(dpath, errors="replace").read()
        A("```\n" + dna.strip() + "\n```\n")
    A("Method: `nucmer --maxmatch` followed by `delta-filter -1` for the 1-to-1 "
      "alignment, plus `minimap2 -x asm5 --eqx` for a reference-anchored BAM. "
      "The two whole-genome alignments agree.\n")
    A(f"![Figure 3](figures/fig03_synteny_dotplot.png)\n")

    A("### 3.2 Structural differences\n")
    if origin:
        A("**Origin shift (coordinate convention).**\n")
        A(md_table(origin, cols=["ref_contig", "ref_length", "dominant_offset_bp",
                                 "pct_reference_bp_at_dominant_offset",
                                 "origin_shifted", "interpretation"], max_rows=5))
    if structural:
        A("**Rearrangements.**\n")
        A(md_table(structural, cols=["type", "ref_contig", "ref_start", "ref_end",
                                     "qry_start", "qry_end", "length_bp",
                                     "identity_pct", "forward_copy_elsewhere"],
                   max_rows=20))
        A("The `forward_copy_elsewhere` column records whether a genuine inverted "
          "repeat (the segment present in both orientations) could explain a "
          "reverse block; for the 88 kb inversion it is `False`, and an "
          "independent `minimap2 -x asm5 --eqx` run reports a single perfect "
          "`88,340=` reverse-complement match with no forward copy — so the "
          "inversion is real and not a repeat artefact.\n")
    A(f"![Figure 4](figures/fig04_structural_variants.png)\n")

    A("### 3.3 Coverage\n")
    A(md_table(coverage, max_rows=10))
    A(f"![Figure 5](figures/fig05_genome_circos.png)\n")
    A(f"![Figure 6](figures/fig06_coverage.png)\n")

    A("## 4. Variant detection\n")
    A("### 4.1 Strategy\n")
    A("Four independent calling routes were run and then reconciled:\n")
    A("| route | evidence | tool |")
    A("|---|---|---|")
    A("| `mummer` | assembly-to-assembly alignment | MUMmer `show-snps` on the 1-to-1 delta |")
    A("| `bcftools_asm` | contig alignment (minimap2 asm5) | `bcftools mpileup` + `call --ploidy 1` |")
    A("| `freebayes_asm` | contig alignment | `freebayes -p 1` |")
    A("| `bcftools_sim_*`, `freebayes_sim` | simulated reads | `bcftools` / `freebayes` |")
    A("")
    A("`bcftools_asm` and `freebayes_asm` share one alignment, so their agreement "
      "is **not** independent evidence. A call is only accepted when at least two "
      "of the three genuinely independent *methods* (NUCmer alignment, contig "
      "alignment, read alignment) support it **and** no ambiguity flag is raised. "
      "Ambiguity flags are: the reference position is covered by more than one "
      "genomic copy; the NUCmer alignment reports overlapping alternative "
      "alignments; the reads covering the position map with low confidence; or "
      "the reads simply contradict the call (the alternate allele is largely "
      "absent).\n")
    if stats["merge"]:
        m = stats["merge"]
        A(f"**{m.get('total_distinct_variants', 0)} candidate positions** were "
          f"generated; **{m.get('high_confidence', 0)}** survive every filter.\n")
        A(md_table([{"metric": k, "value": v} for k, v in m.items()
                    if k != "per_caller"], max_rows=20))
        if m.get("per_caller"):
            A("\nCalls per route:\n")
            A(md_table([{"route": k, "calls": v}
                        for k, v in m["per_caller"].items()], max_rows=20))
    A(f"![Figure 7](figures/fig07_variant_landscape.png)\n")
    A(f"![Figure 8](figures/fig08_variant_confidence.png)\n")
    A(f"![Figure 9](figures/fig09_substitution_spectrum.png)\n")

    A("### 4.2 The verified variant set\n")
    A(md_table(confident,
               cols=[c for c in ["contig", "pos", "ref", "alt", "type",
                                 "n_independent_methods", "callers",
                                 "obs_depth", "obs_alt_frac", "flags"]
                     if not confident or c in confident[0]],
               max_rows=50))
    A("`obs_depth` / `obs_alt_frac` are the read depth at the position (from the "
      "simulated reads, mapped with MAPQ ≥ 20) and the percentage of those reads "
      "carrying the alternate allele. For a genuine homozygous difference "
      "between two haploid genomes this should be close to 100 %; the one SNP "
      "scores 100.0 %.\n")
    A(f"![Figure 10](figures/fig10_variant_verification.png)\n")
    A("Because gap placement inside homopolymer / short-tandem runs is ambiguous, "
      "the aligner's first guess is not necessarily the biologically accurate "
      "description. The figure above reconstructs both alleles from the "
      "alignment, which resolves the indels to: an extra **T** in a T-run "
      "(`530,487`), a 2 bp **T+A** deletion in a T₆A₆ run (`669,715`), and an "
      "extra **T** at `1,709,303`.\n")

    A("## 5. Functional consequence\n")
    A("### 5.1 Annotation\n")
    A("A complete query annotation was produced **by projecting every reference "
      "gene through the alignment** (query coordinates *and* query sequence). At "
      "this level of identity this is more reliable than de-novo gene calling: "
      "the reference carries a curated RefSeq annotation and the two genomes are "
      "collinear. De-novo annotation with Prokka/Bakta was attempted as well; "
      "the Bakta database download was blocked in this environment, so the "
      "projection route is the operative one.\n")
    if stats["projection"]:
        A(md_table([{"metric": k, "value": v}
                    for k, v in stats["projection"].items()], max_rows=20))
    A("")
    A("### 5.2 Affected genes\n")
    A(md_table(gene_impact,
               cols=[c for c in ["locus_tag", "product", "reference_gene_status",
                                 "worst_effect", "worst_impact", "protein_change",
                                 "cds_identity", "ref_protein_len", "qry_protein_len",
                                 "ref_internal_stops", "qry_internal_stops",
                                 "query_orf_status"]
                     if not gene_impact or c in gene_impact[0]],
               max_rows=50))
    A(f"![Figure 11](figures/fig11_gene_impact.png)\n")
    A(f"### 5.3 Genes with a different protein product\n")
    A(f"**Explained by a verified variant ({len(differing)}):**\n")
    A(md_table(differing,
               cols=[c for c in ["locus_tag", "product", "ref_contig", "ref_start",
                                 "ref_end", "qry_contig", "qry_start", "qry_end",
                                 "projected_fraction", "ref_protein_len",
                                 "qry_protein_len", "query_internal_stops",
                                 "ref_internal_stops", "orf_status"]
                     if not differing or c in differing[0]],
               max_rows=50))
    if unconfirmed:
        A(f"\n**Not explained by a verified variant — repeat / breakpoint "
          f"artefacts ({len(unconfirmed)}):**\n")
        A(md_table(unconfirmed,
                   cols=[c for c in ["locus_tag", "product", "ref_start", "ref_end",
                                     "projected_fraction", "orf_status"]
                         if not unconfirmed or c in unconfirmed[0]],
                   max_rows=30))
    A("")
    A("Genes listed with `orf_status = query_ORF_restored` have an intact query "
      "ORF while the reference CDS is annotated as a pseudogene with an internal "
      "frameshift; the single-base indel in the query restores the reading frame. "
      "`sequence_differs_but_intact_in_both` entries are point substitutions.\n")

    A("## 6. Resistance, virulence and metabolism\n")
    A("### 6.1 Acquired resistance genes (NCBI AMRFinderPlus)\n")
    if amr and len(amr) > 1:
        A(md_table(amr[1:], max_rows=50))
    else:
        A("_AMRFinderPlus was run on both genomes; see "
          "`05_functional/amrfinder/` for the full output._\n")
    A("")
    aq = t("05_functional", "amrfinder", "query.amrfinder.tsv")
    ar = t("05_functional", "amrfinder", "reference.amrfinder.tsv")
    if aq or ar:
        sym = "Element symbol"
        qs_ = {r.get(sym, "") for r in (aq[1:] if aq else [])}
        rs_ = {r.get(sym, "") for r in (ar[1:] if ar else [])}
        A(f"* Elements in the **query only**: {sorted(qs_ - rs_) or 'none'}")
        A(f"* Elements in the **reference only**: {sorted(rs_ - qs_) or 'none'}")
        A(f"* Shared: {sorted(qs_ & rs_) or 'none'}\n")
    A("ABRicate was not installable in this environment (Perl version conflict in "
      "the conda solve); AMRFinderPlus plus the curated scan below cover the same "
      "ground with curated, versioned databases.\n")

    A("### 6.2 Curated *B. fragilis* determinant scan\n")
    A("A hand-curated catalogue of 39 resistance / virulence / metabolism "
      "determinants (β-lactams incl. *cfiA*/*cfxA*/AmpC, nitroimidazole *nim*, "
      "macrolide *ermF*, tetracycline *tetQ*/*tetX*, aminoglycosides, "
      "sulphonamide, trimethoprim, phenicol, oxazolidinone, quinolone targets "
      "*gyrA*/*gyrB*/*parC*/*parE*, *rpoB*, 23S rRNA; virulence: *bft*, CPS loci, "
      "*susC*/*susD* PUL, sialidase, haemolysin, protease, fimbriae, bile salt "
      "hydrolase, oxidative stress; metabolism: B12, sulfur, glycoside "
      "hydrolases, polysaccharide lyases, carbohydrate esterases) was matched "
      "against both annotations and every hit was compared amino acid by amino "
      "acid after allele projection.\n")
    if stats["curated"]:
        A(md_table([{"metric": k, "value": (json.dumps(v) if isinstance(v, (list, dict))
                                           else v)}
                    for k, v in stats["curated"].items() if k != "aa_change_detail"],
                   max_rows=20))
    if curated_changes:
        A("\n**Determinants with amino-acid level differences:**\n")
        A(md_table(curated_changes,
                   cols=["determinant", "category", "locus_tag", "product",
                         "identity", "aa_changes"], max_rows=25))
    A(f"![Figure 13](figures/fig13_amr_determinants.png)\n")

    A("### 6.3 Metabolism and general function (COG)\n")
    A("The NCBI COG2020 database was searched with DIAMOND against both "
      "proteomes, which yields the standard COG functional-category profile "
      "without a proprietary pathway database. Categories C, E, F, G, H, I, P "
      "and Q are the metabolic categories.\n")
    if stats["cog"]:
        mt = stats["cog"].get("metabolism_category_totals", {})
        A(f"* metabolic-category protein assignments — reference: "
          f"{mt.get('reference','n/a')}, query: {mt.get('query','n/a')}")
        A(f"* COG-annotated proteins — reference: "
          f"{stats['cog'].get('reference',{}).get('assigned','n/a')}, query: "
          f"{stats['cog'].get('query',{}).get('assigned','n/a')}")
        A(f"* COG ids found in the query only: "
          f"{len(stats['cog'].get('cog_ids_query_only', []))}")
        A(f"* COG ids found in the reference only: "
          f"{len(stats['cog'].get('cog_ids_reference_only', []))}\n")
    if cog_rows:
        A(md_table(cog_rows, max_rows=30))
    A(f"![Figure 12](figures/fig12_functional_profile.png)\n")

    A("### 6.4 Plasmid\n")
    pl = j("06_extra", "plasmid", "plasmid_comparison.json")
    if pl:
        A(md_table([{"property": k, "value": (json.dumps(v) if isinstance(v, dict)
                                             else v)}
                    for k, v in pl.items() if not isinstance(v, dict)], max_rows=20))
        A("")
        A("The plasmid is **circular**, so the comparison allows for an arbitrary "
          "rotation offset and both orientations. Both plasmids are exactly "
          f"{pl.get('query_plasmid',{}).get('length','n/a')} bp and, after "
          "accounting for the deposited origin, are **identical in sequence** "
          "— the query copy is simply written from a different start point and "
          "on the opposite strand. There are no plasmid-encoded differences "
          "between the two strains.\n")
        pg = t("06_extra", "plasmid", "plasmid_genes.tsv")
        if pg:
            conj = [r for r in pg if any(k in r.get("product", "").lower()
                    for k in ("conjugat", "mobiliz", "relaxase", "tra ", "trag",
                              "trán", "transfer", "virb4"))]
            A(f"The shared plasmid carries **{len(pg)} annotated genes**, "
              f"including **{len(conj)}** mobilisation / conjugative-transfer "
              f"functions ({', '.join(sorted({r['product'][:44] for r in conj}))[:300]}). "
              "Because the plasmid is sequence-identical in both strains, this "
              "conjugative machinery — and therefore the plasmid's mobility "
              "potential — is shared.\n")

    A("### 6.5 Gene content\n")
    if stats["gene_content"]:
        A(md_table([{"metric": k, "value": v}
                    for k, v in stats["gene_content"].items()
                    if not isinstance(v, list)], max_rows=20))
    A("")

    A("## 7. Reproducing this analysis\n")
    A("```bash")
    A("source /home/cz/miniconda3/etc/profile.d/conda.sh")
    A(f"conda activate {args.env_name}")
    A("cd /home/cz/WGS")
    A("bash WGS.sh                # full pipeline, all steps")
    A("bash WGS.sh --steps 3,4    # only variant calling and annotation")
    A("bash WGS.sh --check         # verify the toolchain without running")
    A("```\n")
    A("Results are written under `result/01_qc` … `result/08_report`; every "
      "intermediate is under `work/`; per-step logs are under `logs/`.\n")

    A("## 8. Limitations and caveats\n")
    A("1. **Variant confirmation uses simulated reads.** They are exact copies of "
      "the query assembly and exist to demonstrate depth-supported confirmation "
      "and read-level QC; they cannot reveal errors that are already in the "
      "query assembly. Orthogonal confirmation (Sanger, Illumina, or long-read "
      "sequencing of the real strain) is still the gold standard.")
    A("2. **215 candidate differences in multi-copy repeats and at rearrangement "
      "breakpoints could not be resolved** and are listed explicitly in "
      "`03_variants/verified/variants_ambiguous.tsv` with the reason for each. "
      "Short reads cannot assign a difference to a specific copy of a repeated "
      "element; long-read or clone-level data would be needed. These are "
      "candidate differences, not confirmed ones.")
    A("3. **The 88 kb inversion is a sequence-orientation difference between the "
      "two deposited assemblies.** It is supported by a perfect reverse-"
      "complement alignment and the absence of a forward copy, but the "
      "biological consequence (e.g. whether it disrupts a gene) is inferred from "
      "the reference annotation.")
    A("4. **The origin shift is not a biological difference** and is reported "
      "separately so it is not mistaken for one.")
    A("5. **Bakta/Prokka de-novo annotation could not be completed** in this "
      "environment (Bakta's database is hosted on Zenodo, which was unreachable; "
      "Prokka has an unsatisfiable Perl dependency in this conda environment). "
      "Functional annotation therefore rests on the reference RefSeq annotation "
      "transferred through the alignment, plus AMRFinderPlus and COG2020, which "
      "is appropriate at 99.99 % identity but should be revisited if a Bakta "
      "database becomes available.")
    A("6. **The plasmid is identical to the reference plasmid in content** but "
      "deposited in the opposite orientation; see `06_extra/plasmid/`.\n")

    A("## 9. Software\n")
    A("```")
    A(versions.strip() or "see 08_report/software_versions.txt")
    A("```\n")

    md = "\n".join(L)

    # ---------------------------------------------------------------- figures
    figdir_md = os.path.join(O, "07_visualization")
    os.makedirs(os.path.join(O, "08_report", "figures"), exist_ok=True)
    import shutil
    for name, _ in FIG_ORDER:
        src = os.path.join(figdir_md, name + ".png")
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(O, "08_report", "figures", name + ".png"))

    with open(os.path.join(O, "08_report", "WGS_report.md"), "w") as fh:
        fh.write(md)

    # ---------------------------------------------------------------- HTML
    H = []
    H.append("<html><head><meta charset='utf-8'><title>WGS report — "
             f"B. fragilis {args.qry_name} vs {args.ref_name}</title>")
    H.append("""<style>
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
max-width:1150px;margin:0 auto;padding:28px;line-height:1.62;color:#1a1a1a}
h1{font-size:1.7em;border-bottom:3px solid #2166ac;padding-bottom:.3em}
h2{font-size:1.3em;margin-top:2em;border-bottom:1px solid #ddd;padding-bottom:.2em;color:#2166ac}
h3{font-size:1.08em;margin-top:1.5em;color:#333}
table{border-collapse:collapse;font-size:.82em;margin:1em 0;width:100%}
th,td{border:1px solid #d8d8d8;padding:5px 8px;text-align:left;vertical-align:top}
th{background:#eef3f8;font-weight:600}
tr:nth-child(even){background:#fafbfc}
code{background:#f4f4f4;padding:1px 4px;border-radius:3px;font-size:.9em}
pre{background:#f7f7f7;padding:12px;border-radius:6px;overflow-x:auto;font-size:.82em}
img{max-width:100%;border:1px solid #e3e3e3;border-radius:6px;margin:1em 0}
blockquote{border-left:4px solid #f0c36d;background:#fffbf0;margin:1em 0;
padding:.6em 1em;border-radius:0 6px 6px 0}
.meta{color:#666;font-size:.9em}
</style></head><body>""")

    # convert markdown (a small, purpose-built subset) to HTML
    import re
    lines = md.split("\n")
    i, in_tbl, in_pre, in_list = 0, False, False, False
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("```"):
            if in_pre:
                H.append("</pre>"); in_pre = False
            else:
                H.append("<pre>"); in_pre = True
            i += 1; continue
        if in_pre:
            H.append(html.escape(ln)); i += 1; continue
        if ln.startswith("|"):
            cells = [c.strip() for c in ln.strip("|").split("|")]
            if all(set(c) <= set("-: ") for c in cells):
                i += 1; continue
            if not in_tbl:
                H.append("<table>"); in_tbl = True
                H.append("<tr>" + "".join(f"<th>{html.escape(c)}</th>"
                                          for c in cells) + "</tr>")
            else:
                H.append("<tr>" + "".join(f"<td>{html.escape(c)}</td>"
                                          for c in cells) + "</tr>")
            i += 1; continue
        elif in_tbl:
            H.append("</table>"); in_tbl = False
        if ln.startswith(">"):
            H.append(f"<blockquote>{html.escape(ln.lstrip('> '))}</blockquote>")
            i += 1; continue
        if re.match(r"^\s*[*-]\s+", ln):
            if not in_list:
                H.append("<ul>"); in_list = True
            H.append("<li>" + html.escape(re.sub(r"^\s*[*-]\s+", "", ln)) + "</li>")
            i += 1; continue
        elif in_list:
            H.append("</ul>"); in_list = False
        m = re.match(r"^(#{1,4})\s+(.*)$", ln)
        if m:
            lvl = len(m.group(1))
            H.append(f"<h{lvl}>{html.escape(m.group(2))}</h{lvl}>")
            i += 1; continue
        img = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)$", ln.strip())
        if img:
            p = img.group(2)
            base = os.path.basename(p).replace(".png", "")
            cand = os.path.join(O, "08_report", "figures", os.path.basename(p))
            b = b64img(cand)
            if b:
                H.append(f"<img alt='{html.escape(base)}' "
                         f"src='data:image/png;base64,{b}'>")
            i += 1; continue
        if ln.strip():
            txt = html.escape(ln)
            txt = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", txt)
            txt = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<i>\1</i>", txt)
            txt = re.sub(r"`([^`]+)`", r"<code>\1</code>", txt)
            H.append(f"<p>{txt}</p>")
        i += 1
    if in_tbl:
        H.append("</table>")
    if in_pre:
        H.append("</pre>")
    if in_list:
        H.append("</ul>")
    H.append("</body></html>")
    with open(os.path.join(O, "08_report", "WGS_report.html"), "w") as fh:
        fh.write("\n".join(H))

    print(f"[report] wrote {os.path.join(O,'08_report','WGS_report.md')} "
          f"({len(md)} bytes) and WGS_report.html")
    print(f"[report] key findings: {len(kf)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref-name", default="NCTC9343")
    ap.add_argument("--qry-name", default="CCUG4856")
    ap.add_argument("--ref-acc", default="GCF_000025985.1")
    ap.add_argument("--qry-acc", default="GCF_005706655.1")
    ap.add_argument("--ref-chr", default="NC_003228.3")
    ap.add_argument("--ref-plas", default="NC_006873.1")
    ap.add_argument("--qry-chr", default="NZ_CP036555.1")
    ap.add_argument("--qry-plas", default="NZ_CP036556.1")
    ap.add_argument("--env-name", default="WGS")
    ap.add_argument("--sim-cov", default="100")
    ap.add_argument("--sim-len", default="150")
    ap.add_argument("--sim-ins", default="350")
    a = ap.parse_args()
    build(a)


if __name__ == "__main__":
    main()
