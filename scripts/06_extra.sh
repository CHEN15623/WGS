#!/usr/bin/env bash
# =============================================================================
#  06_extra.sh -- ADDITIONAL ANALYSES (optional, degrade gracefully)
# =============================================================================
#    (a) plasmid comparison   : the two 36.5 kb plasmids head to head
#    (b) prophage detection   : PhiSpy on both genomes
#    (c) CRISPR arrays        : minced + cas gene scan
#    (d) mobile elements      : IS / transposase landscape
#    (e) phylogenomics        : download B. fragilis RefSeq genomes, ANI,
#                               core-genome SNP tree (IQ-TREE)
#    (f) pan-genome           : Roary core/accessory/unique gene statistics
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

EX="$OUT/06_extra"
mkdir -p "$EX/plasmid" "$EX/prophage" "$EX/crispr" "$EX/mobile" \
         "$EX/phylogeny" "$EX/pangenome" "$WORK/phylo"

# =============================================================================
msg "06-EX (a): plasmid comparison (reference plasmid vs query plasmid)"
# =============================================================================
python3 "$SCRIPTS/lib/plasmid_compare.py" \
  --ref-fna "$REF_FA" --qry-fna "$QRY_FA" \
  --ref-gff "$REF_GFF" \
  --bam "$WORK/bam/${QRY_NAME}_to_${REF_NAME}.asm5.bam" \
  --outdir "$EX/plasmid" || warn "plasmid comparison failed"

# =============================================================================
msg "06-EX (b): prophage detection"
# =============================================================================
if command -v phispy >/dev/null 2>&1; then
  for g in "reference:$REF_FA" "query:$QRY_FA"; do
    label="${g%%:*}"; fa="${g#*:}"
    d="$EX/prophage/$label"
    [ -s "$d/prophage.tsv" ] && continue
    mkdir -p "$d"
    run_logged "phispy_$label" bash -c \
      "cd '$d' && phispy --threads $THREADS --outputdir '$d/phispy_out' '$fa' \
         > '$d/phispy.log' 2>&1" || warn "phispy $label failed"
  done
else
  warn "phispy not installed -- prophage step skipped"
fi

# =============================================================================
msg "06-EX (c): CRISPR arrays and cas genes"
# =============================================================================
if command -v minced >/dev/null 2>&1; then
  for g in "reference:$REF_FA" "query:$QRY_FA"; do
    label="${g%%:*}"; fa="${g#*:}"
    run_logged "minced_$label" bash -c \
      "minced -minNR 3 -spacers '$fa' '$EX/crispr/${label}_crispr.txt'" || \
      warn "minced $label failed"
  done
fi
# cas genes are found by keyword in whichever annotation is available
python3 - "$REF_GFF" "$OUT/04_annotation/prokka_query/${QRY_NAME}.gff" \
           "$EX/crispr/cas_genes.tsv" <<'PY' || true
import sys, os
ref_gff, qry_gff, out = sys.argv[1], sys.argv[2], sys.argv[3]
KEYS = ("crispr", "cas1", "cas2", "cas3", "cas4", "cas5", "cas6", "cas7",
        "cas8", "cas9", "cas10", "csm", "cmr", "csy", "cse", "cascade")
rows = []
for label, path in (("reference", ref_gff), ("query", qry_gff)):
    if not os.path.exists(path):
        continue
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        hay = (attrs.get("product", "") + " " + attrs.get("gene", "")).lower()
        if any(k in hay for k in KEYS):
            rows.append((label, f[0], f[3], f[4], attrs.get("locus_tag", ""),
                         attrs.get("gene", ""), attrs.get("product", "")))
with open(out, "w") as fh:
    fh.write("genome\tcontig\tstart\tend\tlocus_tag\tgene\tproduct\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")
print(f"[crispr] {len(rows)} cas-like genes")
PY

# =============================================================================
msg "06-EX (d): mobile genetic elements / transposase landscape"
# =============================================================================
python3 - "$REF_GFF" "$OUT/04_annotation/prokka_query/${QRY_NAME}.gff" \
           "$EX/mobile/transposase_landscape.tsv" <<'PY' || true
import sys, os
ref_gff, qry_gff, out = sys.argv[1], sys.argv[2], sys.argv[3]
KEYS = ("transposase", "integrase", "recombinase", "insertion sequence",
        "conjugal", "relaxase", "mobilization", "resolvase", "is element")
with open(out, "w") as fh:
    fh.write("genome\tcontig\tstart\tend\tstrand\tlocus_tag\tgene\tproduct\telement_class\n")
    n = 0
    for label, path in (("reference", ref_gff), ("query", qry_gff)):
        if not os.path.exists(path):
            continue
        for line in open(path):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            hay = (attrs.get("product", "") + " " + attrs.get("gene", "")).lower()
            cls = next((k for k in KEYS if k in hay), None)
            if cls:
                fh.write("\t".join([label, f[0], f[3], f[4], f[6],
                                    attrs.get("locus_tag", ""),
                                    attrs.get("gene", ""),
                                    attrs.get("product", ""), cls]) + "\n")
                n += 1
    print(f"[mobile] {n} mobile-element genes")
PY

# =============================================================================
msg "06-EX (e): phylogenomic placement among B. fragilis genomes"
# =============================================================================
PHYLO_GENOMES="$WORK/phylo/genomes"
if [ ! -d "$PHYLO_GENOMES" ] || [ "$(ls -1 "$PHYLO_GENOMES" 2>/dev/null | wc -l)" -lt 5 ]; then
  msg "  downloading B. fragilis RefSeq assemblies ..."
  timeout 1800 python3 "$SCRIPTS/lib/fetch_bfragilis.py" \
    --outdir "$PHYLO_GENOMES" --max 25 --taxid 817 \
    > "$LOGDIR/fetch_genomes.log" 2>&1 || warn "genome download failed or timed out"
fi

if [ "$(ls -1 "$PHYLO_GENOMES" 2>/dev/null | wc -l)" -ge 5 ]; then
  # include the reference and the query in the panel
  mkdir -p "$WORK/phylo/panel"
  cp "$REF_FA" "$WORK/phylo/panel/reference_${REF_NAME}.fna"
  cp "$QRY_FA" "$WORK/phylo/panel/query_${QRY_NAME}.fna"
  cp "$PHYLO_GENOMES"/*.fna "$WORK/phylo/panel/" 2>/dev/null || true

  # --- ANI / Mash distance ------------------------------------------------
  if command -v fastANI >/dev/null 2>&1; then
    qlist="$WORK/phylo/query_list.txt"
    ls "$WORK/phylo/panel"/*.fna | grep -v "reference_${REF_NAME}" > "$qlist"
    run_logged fastani bash -c \
      "fastANI -q '$qlist' -r '$REF_FA' -o '$EX/phylogeny/fastANI_vs_reference.tsv' -t $THREADS" || \
      warn "fastANI failed"
  fi
  if command -v mash >/dev/null 2>&1; then
    run_logged mash_sketch bash -c \
      "mash sketch -p $THREADS -o '$WORK/phylo/mash_sketch' '$WORK/phylo/panel'/*.fna"
    run_logged mash_dist bash -c \
      "mash dist '$WORK/phylo/mash_sketch.msh' '$WORK/phylo/mash_sketch.msh' \
         > '$EX/phylogeny/mash_distances.tsv'" || true
  fi

  # --- core-genome SNP phylogeny via Prokka + Roary + MAFFT + IQ-TREE ----
  if command -v prokka >/dev/null 2>&1 && command -v roary >/dev/null 2>&1; then
    PROKDIR="$WORK/phylo/prokka"
    mkdir -p "$PROKDIR"
    if [ "$(ls -1 "$PROKDIR"/*.gff 2>/dev/null | wc -l)" -lt 5 ]; then
      msg "  annotating the panel with Prokka (this takes a few minutes) ..."
      timeout 2400 run_logged prokka_panel bash -c \
        "ls '$WORK/phylo/panel'/*.fna | xargs -P $((THREADS/4)) -I{} bash -c \
           'b=\$(basename {} .fna); prokka --quiet --outdir \"$PROKDIR/\$b\" \
              --prefix \"\$b\" --cpus 4 --kingdom Bacteria --gcode 11 --force {}'" || \
        warn "panel Prokka failed"
    fi
    if [ "$(ls -1 "$PROKDIR"/*/*.gff 2>/dev/null | wc -l)" -ge 5 ]; then
      ROARY_OUT="$WORK/phylo/roary"
      [ -s "$ROARY_OUT/core_gene_alignment.aln" ] || run_logged roary bash -c \
        "roary -p $THREADS -e --mafft -f '$ROARY_OUT' '$PROKDIR'/*/*.gff" || \
        warn "roary failed"
      cp "$ROARY_OUT/gene_presence_absence.csv" "$EX/pangenome/" 2>/dev/null || true
      cp "$ROARY_OUT/summary_statistics.txt"   "$EX/pangenome/" 2>/dev/null || true
      if [ -s "$ROARY_OUT/core_gene_alignment.aln" ]; then
        if command -v iqtree >/dev/null 2>&1 || command -v iqtree2 >/dev/null 2>&1; then
          IQ=$(command -v iqtree2 || command -v iqtree)
          run_logged iqtree bash -c \
            "'$IQ' -s '$ROARY_OUT/core_gene_alignment.aln' -m MFP -B 1000 \
                 -T $MAXTHREADS -pre '$EX/phylogeny/core_snp_tree' -quiet" || \
            warn "IQ-TREE failed"
        fi
        # also produce a SNP-only alignment for reference
        if command -v snp-sites >/dev/null 2>&1; then
          run_logged snpsites bash -c \
            "snp-sites -c -o '$EX/phylogeny/core_snp_alignment.fa' \
               '$ROARY_OUT/core_gene_alignment.aln'" || true
        fi
      fi
    fi
  else
    warn "prokka/roary not available -- core-genome tree skipped"
  fi
else
  warn "not enough genomes downloaded -- phylogenomics skipped"
fi

# =============================================================================
msg "06-EX (f): pan-genome statistics"
# =============================================================================
if [ -s "$EX/pangenome/summary_statistics.txt" ]; then
  python3 "$SCRIPTS/lib/pangenome_stats.py" \
    --roary-summary "$EX/pangenome/summary_statistics.txt" \
    --presence-absence "$EX/pangenome/gene_presence_absence.csv" \
    --outdir "$EX/pangenome" || warn "pan-genome stats failed"
fi

msg "06-Extra analyses complete -> $EX"
