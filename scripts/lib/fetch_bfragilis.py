#!/usr/bin/env python3
"""
fetch_bfragilis.py -- download finished B. fragilis RefSeq assemblies to place
the query strain in a phylogenetic context.

Reads the per-species NCBI assembly_summary, keeps the newest version of each
assembly that is at least "Contig" level (preferring Complete Genome and
representative/reference assemblies), and fetches `*_genomic.fna.gz` over HTTPS.
"""
import argparse
import gzip
import os
import re
import sys
import time
import urllib.request

SUMMARY_URLS = [
    "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/"
    "Bacteroides_fragilis/assembly_summary.txt",
    "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/"
    "Bacteroides_fragilis/assembly_summary.txt",
]
FALLBACK_ASSEMBLY = [
    # a hand-picked, well-known panel used if the FTP index is unreachable
    ("GCF_000025985.1", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/025/985/"
     "GCF_000025985.1_ASM2598v1/GCF_000025985.1_ASM2598v1_genomic.fna.gz"),
    ("GCF_005706655.1", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/005/706/655/"
     "GCF_005706655.1_ASM570665v1/GCF_005706655.1_ASM570665v1_genomic.fna.gz"),
]


def fetch_summary():
    for url in SUMMARY_URLS:
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as exc:
            print(f"[fetch] {url} failed: {exc}", file=sys.stderr)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max", type=int, default=40)
    ap.add_argument("--taxid", default="817")
    ap.add_argument("--min-level", default="Contig")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    text = fetch_summary()
    rows = []
    if text:
        header = None
        for line in text.splitlines():
            if line.startswith("# assembly_accession"):
                header = line.lstrip("# ").split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if header is None or len(f) < len(header):
                continue
            d = dict(zip(header, f))
            if d.get("taxid") != a.taxid:
                continue
            if d.get("assembly_status") not in ("current", ""):
                continue
            if d.get("version_status") != "latest":
                continue
            rows.append(d)

    # rank: complete genome first, then reference/representative, then size
    LEVEL = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    CAT = {"reference genome": 0, "representative genome": 1, "na": 2}

    def key(d):
        return (LEVEL.get(d.get("assembly_level", "Contig"), 9),
                CAT.get(d.get("refseq_category", "na"), 9),
                -int(d.get("genome_size", 0) or 0))

    rows.sort(key=key)

    # keep one assembly per BioSample/BioProject to reduce redundancy
    seen = set()
    chosen = []
    for d in rows:
        k = d.get("biosample", "")  or d.get("bioproject", "")
        if k in seen:
            continue
        seen.add(k)
        chosen.append(d)
        if len(chosen) >= a.max:
            break

    if not chosen:
        print("[fetch] assembly index unavailable -- using the fallback panel",
              file=sys.stderr)
        for acc, url in FALLBACK_ASSEMBLY:
            out = os.path.join(a.outdir, f"{acc}.fna")
            if os.path.exists(out):
                continue
            try:
                urllib.request.urlretrieve(url, out + ".gz")
                with gzip.open(out + ".gz", "rb") as fi, open(out, "wb") as fo:
                    fo.write(fi.read())
                os.remove(out + ".gz")
                print(f"[fetch] {acc}")
            except Exception as exc:
                print(f"[fetch] {acc} failed: {exc}", file=sys.stderr)
        return

    manifest = os.path.join(a.outdir, "panel_manifest.tsv")
    with open(manifest, "w") as mf:
        mf.write("assembly_accession\torganism\tstrain\tspecies_taxid\t"
                 "assembly_level\trefseq_category\tgenome_size\tgc_percent\t"
                 "ftp_path\tlocal_file\n")
        for d in chosen:
            acc = d["assembly_accession"]
            ftp = d.get("ftp_path", "")
            if not ftp or ftp == "na":
                continue
            base = ftp.rstrip("/").split("/")[-1]
            url = f"{ftp.replace('ftp://', 'https://')}/{base}_genomic.fna.gz"
            out = os.path.join(a.outdir, f"{acc}.fna")
            if not os.path.exists(out) or os.path.getsize(out) < 1000:
                try:
                    tmp = out + ".gz"
                    urllib.request.urlretrieve(url, tmp)
                    with gzip.open(tmp, "rb") as fi, open(out, "wb") as fo:
                        fo.write(fi.read())
                    os.remove(tmp)
                    print(f"[fetch] {acc}  {d.get('organism_name','')}")
                except Exception as exc:
                    print(f"[fetch] {acc} failed: {exc}", file=sys.stderr)
                    if os.path.exists(out) and os.path.getsize(out) < 1000:
                        os.remove(out)
                    continue
            mf.write("\t".join([
                acc, d.get("organism_name", ""), d.get("infraspecific_name", ""),
                d.get("taxid", ""), d.get("assembly_level", ""),
                d.get("refseq_category", ""), d.get("genome_size", ""),
                d.get("gc_percent", ""), ftp, out]) + "\n")
    n = len([f for f in os.listdir(a.outdir) if f.endswith(".fna")])
    print(f"[fetch] {n} genomes available in {a.outdir}")


if __name__ == "__main__":
    main()
