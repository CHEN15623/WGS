#!/usr/bin/env python3
"""
functional_compare.py -- consolidate the functional-annotation evidence into
reference-vs-query comparisons along the three requested axes.

  RESISTANCE  : ABRicate (CARD / ResFinder / NCBI / ARG-ANNOT / MEGARES /
                PlasmidFinder / ISfinder / Victors) and NCBI AMRFinderPlus
  VIRULENCE   : ABRicate VFDB / Victors
  METABOLISM  : Bakta COG / KEGG / GO / EC profiles
  plus a COG functional-category breakdown (needs cog-20.def.tab; downloaded
  automatically when the network allows).

Outputs (in --outdir)
  amr_gene_comparison.tsv      every AMR/virulence gene call, ref vs query
  amr_differences.tsv          genes present in only one genome
  cog_category_comparison.tsv  COG functional category counts, ref vs query
  kegg_ko_comparison.tsv       KEGG KO presence differences
  go_term_comparison.tsv       GO term presence differences
  ec_class_comparison.tsv      EC class counts
  functional_summary.json
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
import urllib.request
from collections import OrderedDict, Counter, defaultdict

COG_URL = ("https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab")
COG_LETTERS = OrderedDict([
    ("J", "Translation, ribosomal structure and biogenesis"),
    ("A", "RNA processing and modification"),
    ("K", "Transcription"),
    ("L", "Replication, recombination and repair"),
    ("B", "Chromatin structure and dynamics"),
    ("D", "Cell cycle control, cell division, chromosome partitioning"),
    ("Y", "Nuclear structure"),
    ("V", "Defense mechanisms"),
    ("T", "Signal transduction mechanisms"),
    ("M", "Cell wall/membrane/envelope biogenesis"),
    ("N", "Cell motility"),
    ("Z", "Cytoskeleton"),
    ("W", "Extracellular structures"),
    ("U", "Intracellular trafficking, secretion, and vesicular transport"),
    ("O", "Post-translational modification, protein turnover, chaperones"),
    ("C", "Energy production and conversion"),
    ("G", "Carbohydrate transport and metabolism"),
    ("E", "Amino acid transport and metabolism"),
    ("F", "Nucleotide transport and metabolism"),
    ("H", "Coenzyme transport and metabolism"),
    ("I", "Lipid transport and metabolism"),
    ("P", "Inorganic ion transport and metabolism"),
    ("Q", "Secondary metabolites biosynthesis, transport and catabolism"),
    ("R", "General function prediction only"),
    ("S", "Function unknown"),
])


# --------------------------------------------------------------------------- #
def parse_bakta_tsv(path):
    """Bakta .tsv -> feature list.  All functional evidence is taken from the
    DbXrefs field so the parser is insensitive to column re-ordering."""
    feats = []
    if not path or not os.path.exists(path):
        return feats
    header = None
    for line in open(path, errors="replace"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if header is None:
            if f and f[0] == "ID":
                header = f
            continue
        if header is None:
            continue
        f += [""] * (len(header) - len(f))
        d = dict(zip(header, f))
        refs = defaultdict(list)
        raw = d.get("DbXrefs", "") or ""
        for item in re.split(r"[;,]", raw):
            item = item.strip()
            if ":" in item:
                pre, val = item.split(":", 1)
                refs[pre.strip()].append(val.strip())
        feats.append(dict(
            id=d.get("ID", ""), contig=d.get("Contig", ""),
            type=d.get("Type", ""), start=d.get("Start", ""),
            stop=d.get("Stop", ""), strand=d.get("Strand", ""),
            locus=d.get("Locus Tag", ""), gene=d.get("Gene", ""),
            product=d.get("Product", ""), refs=refs,
            cog=refs.get("COG", []), kegg=refs.get("KEGG", []),
            go=refs.get("GO", []), ec=refs.get("EC", []),
            pfam=refs.get("Pfam", []), uniref=refs.get("UniRef", [])))
    return feats


def parse_gff_dbxrefs(path):
    """Fallback: pull COG/KEGG/GO/EC out of a NCBI or Prokka GFF."""
    feats = []
    if not path or not os.path.exists(path):
        return feats
    for line in open(path, errors="replace"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = {}
        for kv in f[8].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                attrs.setdefault(k, v)
        refs = defaultdict(list)
        for raw in (attrs.get("Dbxref", "") + "," + attrs.get("Ontology_term", "")).split(","):
            raw = raw.strip()
            if ":" in raw:
                pre, val = raw.split(":", 1)
                refs[pre.strip()].append(val.strip())
            elif raw.startswith("GO"):
                refs["GO"].append(raw)
        ecs = [x for x in attrs.get("EC_number", "").split(",") if x]
        feats.append(dict(id=attrs.get("ID", ""), contig=f[0], type="CDS",
                          start=f[3], stop=f[4], strand=f[6],
                          locus=attrs.get("locus_tag", ""),
                          gene=attrs.get("gene", ""),
                          product=attrs.get("product", ""), refs=refs,
                          cog=refs.get("COG", []), kegg=refs.get("KEGG", []),
                          go=refs.get("GO", []), ec=ecs or refs.get("EC", []),
                          pfam=refs.get("Pfam", []), uniref=refs.get("UniRef", [])))
    return feats


def load_cog_definitions(cache_dir):
    """COG id -> (category letters, name).  Cached locally."""
    os.makedirs(cache_dir, exist_ok=True)
    local = os.path.join(cache_dir, "cog-20.def.tab")
    if not os.path.exists(local) or os.path.getsize(local) < 1000:
        try:
            urllib.request.urlretrieve(COG_URL, local)
        except Exception as exc:
            print(f"[func] COG definition download failed ({exc})", file=sys.stderr)
            return {}
    mapping = {}
    for line in open(local, errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3:
            mapping[f[0]] = (f[1], f[2])
    return mapping


def parse_abricate(path):
    """abricate .tab -> {(gene): row}"""
    res = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return res
    header = None
    for line in open(path, errors="replace"):
        if line.startswith("#"):
            header = line.lstrip("#").rstrip("\n").split("\t")
            continue
        f = line.rstrip("\n").split("\t")
        if header and len(f) < len(header):
            f += [""] * (len(header) - len(f))
        d = dict(zip(header or [], f))
        gene = d.get("GENE", d.get("#FILE", ""))
        if not gene:
            continue
        key = gene
        res[key] = dict(gene=gene, product=d.get("PRODUCT", ""),
                        db=d.get("DATABASE", ""),
                        identity=d.get("IDENTITY", d.get("%IDENTITY", "")),
                        coverage=d.get("COVERAGE", d.get("%COVERAGE", "")),
                        contig=d.get("SEQUENCE", ""),
                        start=d.get("START", ""), end=d.get("END", ""),
                        accession=d.get("ACCESSION", ""),
                        resistance=d.get("RESISTANCE", ""))
    return res


def parse_amrfinder(path):
    res = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return res
    header = None
    for line in open(path, errors="replace"):
        f = line.rstrip("\n").split("\t")
        if header is None:
            header = f
            continue
        if len(f) < len(header):
            f += [""] * (len(header) - len(f))
        d = dict(zip(header, f))
        name = d.get("Element symbol", d.get("Gene symbol", ""))
        if not name:
            continue
        res.setdefault(name, []).append(d)
    return res


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--qry-bakta"); ap.add_argument("--qry-bakta-tsv")
    ap.add_argument("--ref-bakta-tsv")
    ap.add_argument("--qry-prokka")
    ap.add_argument("--abricate-dir"); ap.add_argument("--amrfinder-dir")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cache", default=None)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    cache = a.cache or os.path.join(os.path.dirname(a.outdir), "..", "..",
                                    "work", "cog_cache")

    # ---- functional annotation profiles ----------------------------------
    ref_feats = parse_bakta_tsv(a.ref_bakta_tsv) or parse_gff_dbxrefs(a.ref_gff)
    qry_feats = parse_bakta_tsv(a.qry_bakta_tsv) or parse_gff_dbxrefs(a.qry_prokka)
    print(f"[func] reference features with functional evidence: {len(ref_feats)}")
    print(f"[func] query features with functional evidence    : {len(qry_feats)}")

    def collect(feats, attr):
        c = Counter()
        for f in feats:
            for v in f.get(attr, []) or []:
                c[v.split(".")[0] if attr == "ec" else v] += 1
        return c

    # ---- COG categories ---------------------------------------------------
    cog_def = load_cog_definitions(cache)
    cog_ref, cog_qry = Counter(), Counter()
    for feats, cnt in ((ref_feats, cog_ref), (qry_feats, cog_qry)):
        for f in feats:
            for cid in f.get("cog", []) or []:
                info = cog_def.get(cid)
                letters = info[0] if info else ""
                for L in (letters or "S"):
                    if L in COG_LETTERS:
                        cnt[L] += 1
    with open(f"{a.outdir}/cog_category_comparison.tsv", "w") as fh:
        fh.write("COG_letter\tCOG_category\treference_count\tquery_count\t"
                 "difference\n")
        for L, name in COG_LETTERS.items():
            r, q = cog_ref.get(L, 0), cog_qry.get(L, 0)
            fh.write(f"{L}\t{name}\t{r}\t{q}\t{q-r:+d}\n")

    # ---- COG identity differences ----------------------------------------
    set_ref = {(f["locus"] or f["id"], c) for f in ref_feats for c in f.get("cog", []) or []}
    with open(f"{a.outdir}/cog_id_comparison.tsv", "w") as fh:
        fh.write("COG_id\tCOG_category\treference_count\tquery_count\tdifference\n")
        ids = sorted(set(list(collect(ref_feats, "cog"))) | set(list(collect(qry_feats, "cog"))))
        rc, qc = collect(ref_feats, "cog"), collect(qry_feats, "cog")
        for cid in ids:
            letters = cog_def.get(cid, ("", ""))[0]
            fh.write(f"{cid}\t{letters}\t{rc.get(cid,0)}\t{qc.get(cid,0)}\t"
                     f"{qc.get(cid,0)-rc.get(cid,0):+d}\n")

    # ---- KEGG KO ----------------------------------------------------------
    kr, kq = collect(ref_feats, "kegg"), collect(qry_feats, "kegg")
    with open(f"{a.outdir}/kegg_ko_comparison.tsv", "w") as fh:
        fh.write("KO\treference_count\tquery_count\tdifference\tstatus\n")
        for ko in sorted(set(kr) | set(kq)):
            r, q = kr.get(ko, 0), kq.get(ko, 0)
            st = "shared" if r and q else ("reference_only" if r else "query_only")
            fh.write(f"{ko}\t{r}\t{q}\t{q-r:+d}\t{st}\n")

    # ---- GO ---------------------------------------------------------------
    gr, gq = collect(ref_feats, "go"), collect(qry_feats, "go")
    with open(f"{a.outdir}/go_term_comparison.tsv", "w") as fh:
        fh.write("GO_term\treference_count\tquery_count\tstatus\n")
        for go in sorted(set(gr) | set(gq)):
            r, q = gr.get(go, 0), gq.get(go, 0)
            st = "shared" if r and q else ("reference_only" if r else "query_only")
            fh.write(f"{go}\t{r}\t{q}\t{st}\n")

    # ---- EC classes -------------------------------------------------------
    er, eq = collect(ref_feats, "ec"), collect(qry_feats, "ec")
    with open(f"{a.outdir}/ec_class_comparison.tsv", "w") as fh:
        fh.write("EC_class\tdescription\treference_count\tquery_count\tdifference\n")
        EC_NAMES = {"1": "Oxidoreductases", "2": "Transferases",
                    "3": "Hydrolases", "4": "Lyases", "5": "Isomerases",
                    "6": "Ligases", "7": "Translocases"}
        for cls, nm in EC_NAMES.items():
            r = sum(v for k, v in er.items() if str(k).startswith(cls))
            q = sum(v for k, v in eq.items() if str(k).startswith(cls))
            fh.write(f"{cls}\t{nm}\t{r}\t{q}\t{q-r:+d}\n")

    # ---- ABRicate comparison ---------------------------------------------
    amr_rows = []
    if a.abricate_dir and os.path.isdir(a.abricate_dir):
        for ref_f in sorted(glob.glob(os.path.join(a.abricate_dir, "reference.*.tab"))):
            db = os.path.basename(ref_f).split(".", 2)[1]
            qry_f = os.path.join(a.abricate_dir, f"query.{db}.tab")
            R = parse_abricate(ref_f)
            Q = parse_abricate(qry_f)
            for gene in sorted(set(R) | set(Q)):
                amr_rows.append(dict(
                    database=db, gene=gene,
                    product=(Q.get(gene) or R.get(gene) or {}).get("product", ""),
                    resistance=(Q.get(gene) or R.get(gene) or {}).get("resistance", ""),
                    in_reference="yes" if gene in R else "no",
                    in_query="yes" if gene in Q else "no",
                    ref_identity=(R.get(gene) or {}).get("identity", ""),
                    qry_identity=(Q.get(gene) or {}).get("identity", ""),
                    ref_coverage=(R.get(gene) or {}).get("coverage", ""),
                    qry_coverage=(Q.get(gene) or {}).get("coverage", ""),
                    status=("shared" if gene in R and gene in Q else
                            "reference_only" if gene in R else "query_only")))
    with open(f"{a.outdir}/amr_gene_comparison.tsv", "w") as fh:
        cols = ["database", "gene", "product", "resistance", "in_reference",
                "in_query", "status", "ref_identity", "qry_identity",
                "ref_coverage", "qry_coverage"]
        fh.write("\t".join(cols) + "\n")
        for r in amr_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    diffs = [r for r in amr_rows if r["status"] != "shared"]
    with open(f"{a.outdir}/amr_differences.tsv", "w") as fh:
        cols = ["database", "gene", "product", "status", "in_reference",
                "in_query", "resistance"]
        fh.write("\t".join(cols) + "\n")
        for r in diffs:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # ---- AMRFinderPlus comparison ----------------------------------------
    amrf_rows = []
    if a.amrfinder_dir and os.path.isdir(a.amrfinder_dir):
        R = parse_amrfinder(os.path.join(a.amrfinder_dir, "reference.amrfinder.tsv"))
        Q = parse_amrfinder(os.path.join(a.amrfinder_dir, "query.amrfinder.tsv"))
        for gene in sorted(set(R) | set(Q)):
            r0 = (R.get(gene) or [{}])[0]
            q0 = (Q.get(gene) or [{}])[0]
            amrf_rows.append(dict(
                element=gene,
                element_type=(q0 or r0).get("Element type", ""),
                element_subtype=(q0 or r0).get("Element subtype", ""),
                class_=(q0 or r0).get("Class", ""),
                subclass=(q0 or r0).get("Subclass", ""),
                method=(q0 or r0).get("Method", ""),
                in_reference="yes" if gene in R else "no",
                in_query="yes" if gene in Q else "no",
                status=("shared" if gene in R and gene in Q else
                        "reference_only" if gene in R else "query_only"),
                ref_pident=(r0 or {}).get("% Identity to reference sequence", ""),
                qry_pident=(q0 or {}).get("% Identity to reference sequence", "")))
    with open(f"{a.outdir}/amrfinder_comparison.tsv", "w") as fh:
        cols = ["element", "element_type", "element_subtype", "class_",
                "subclass", "method", "in_reference", "in_query", "status",
                "ref_pident", "qry_pident"]
        fh.write("\t".join(cols) + "\n")
        for r in amrf_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # ---- summary ----------------------------------------------------------
    summary = dict(
        reference_features=len(ref_feats), query_features=len(qry_feats),
        cog_categories=dict(cog_ref), cog_categories_query=dict(cog_qry),
        n_cog_ids_reference=len(collect(ref_feats, "cog")),
        n_cog_ids_query=len(collect(qry_feats, "cog")),
        n_kegg_ko_reference=len(kr), n_kegg_ko_query=len(kq),
        kegg_reference_only=sorted(set(kr) - set(kq))[:100],
        kegg_query_only=sorted(set(kq) - set(kr))[:100],
        n_go_reference=len(gr), n_go_query=len(gq),
        n_ec_reference=len(er), n_ec_query=len(eq),
        abricate_calls=len(amr_rows),
        abricate_differences=len(diffs),
        amrfinder_calls=len(amrf_rows),
        amrfinder_differences=len([r for r in amrf_rows if r["status"] != "shared"]),
        cog_definitions_loaded=len(cog_def),
    )
    with open(f"{a.outdir}/functional_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[func] COG ids ref={summary['n_cog_ids_reference']} "
          f"qry={summary['n_cog_ids_query']}; KEGG KO ref={len(kr)} qry={len(kq)}")
    print(f"[func] ABRicate calls={len(amr_rows)} differences={len(diffs)}; "
          f"AMRFinderPlus calls={len(amrf_rows)}")


if __name__ == "__main__":
    main()
