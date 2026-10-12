#!/usr/bin/env python3
"""
pangenome_stats.py -- summarise a Roary pan-genome and, most importantly,
report which pan-genome genes distinguish the query strain from the reference.

Outputs (in --outdir)
  pangenome_summary.tsv        core / soft-core / shell / cloud counts
  strain_specific_genes.tsv    genes present in only one of ref / query
  query_unique_genes.tsv       genes found in the query but not the reference
  reference_unique_genes.tsv   genes found in the reference but not the query
  accessory_profile.tsv        per-gene presence across the panel
"""
import argparse
import csv
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def parse_summary(path):
    stats = {}
    if not path or not os.path.exists(path):
        return stats
    for line in open(path):
        if ":" in line:
            k, v = line.split(":", 1)
            stats[k.strip()] = v.strip()
    return stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roary-summary")
    ap.add_argument("--presence-absence")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ref-strain", default="reference_NCTC9343")
    ap.add_argument("--qry-strain", default="query_CCUG4856")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    stats = parse_summary(a.roary_summary)
    rows = []
    with open(f"{a.outdir}/pangenome_summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in stats.items():
            fh.write(f"{k}\t{v}\n")
            rows.append((k, v))

    if not a.presence_absence or not os.path.exists(a.presence_absence):
        with open(f"{a.outdir}/pangenome_stats.json", "w") as fh:
            json.dump({"summary": stats}, fh, indent=2)
        print("[pangenome] wrote summary only (no presence/absence matrix)")
        return

    with open(a.presence_absence, newline="", errors="replace") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        # find the first column that looks like a genome name
        try:
            gstart = next(i for i, h in enumerate(header)
                          if h and (h.lower().startswith("gcf_")
                                    or h.lower().startswith("query")
                                    or h.lower().startswith("reference")))
        except StopIteration:
            gstart = 14
        genomes = header[gstart:]
        body = [r for r in rd if r and len(r) > gstart]

    core = shell = cloud = 0
    qonly, ronly, single = [], [], []
    with open(f"{a.outdir}/accessory_profile.tsv", "w") as fh:
        fh.write("gene\tannotation\t" + "\t".join(genomes) + "\tn_genomes\n")
        for r in body:
            gene, annot = r[0], (r[1] if len(r) > 1 else "")
            present = [g for g, v in zip(genomes, r[gstart:gstart + len(genomes)])
                       if v.strip()]
            n = len(present)
            fh.write(f"{gene}\t{annot}\t" +
                     "\t".join("1" if v.strip() else "0"
                               for v in r[gstart:gstart + len(genomes)]) +
                     f"\t{n}\n")
            if n == len(genomes):
                core += 1
            elif n > len(genomes) * 0.95:
                shell += 1
            else:
                cloud += 1
            if n == 1:
                single.append((gene, annot, present[0]))
            if any(g == a.qry_strain for g in present) and \
               not any(g == a.ref_strain for g in present):
                qonly.append((gene, annot, n))
            if any(g == a.ref_strain for g in present) and \
               not any(g == a.qry_strain for g in present):
                ronly.append((gene, annot, n))

    with open(f"{a.outdir}/strain_specific_genes.tsv", "w") as fh:
        fh.write("gene\tannotation\tdirection\tn_genomes_in_panel\n")
        for g, an, n in qonly:
            fh.write(f"{g}\t{an}\tquery_only\t{n}\n")
        for g, an, n in ronly:
            fh.write(f"{g}\t{an}\treference_only\t{n}\n")

    for name, data in (("query_unique_genes.tsv", qonly),
                       ("reference_unique_genes.tsv", ronly)):
        with open(os.path.join(a.outdir, name), "w") as fh:
            fh.write("gene\tannotation\tn_genomes_in_panel\n")
            for g, an, n in data:
                fh.write(f"{g}\t{an}\t{n}\n")

    # core/accessory curve data (rarefaction-free, based on the panel)
    summary = dict(
        roary=stats, n_genomes=len(genomes), n_genes_total=len(body),
        core_genes=core, soft_core_genes=shell, cloud_genes=cloud,
        query_only_genes=len(qonly), reference_only_genes=len(ronly),
        single_genome_genes=len(single),
        genomes=genomes,
    )
    with open(f"{a.outdir}/pangenome_stats.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[pangenome] {len(genomes)} genomes, {len(body)} genes, "
          f"core={core}, cloud={cloud}, query-only={len(qonly)}, "
          f"reference-only={len(ronly)}")


if __name__ == "__main__":
    main()
