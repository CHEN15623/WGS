#!/usr/bin/env bash
# =============================================================================
#  04_annotate.sh -- FUNCTIONAL PREDICTION AND VARIANT GENE ANNOTATION
# =============================================================================
#  (1) De-novo annotation of the query assembly (Prokka + Bakta).  Bakta adds
#      COG / KEGG / GO / EC / Pfam / UniRef cross-references, which is what
#      makes the "metabolism" and "functional prediction" questions answerable.
#  (2) Variant effect prediction with three independent engines:
#        - snpEff   (custom B. fragilis database built from the reference GFF)
#        - bcftools csq (haplotype-aware consequence calling)
#        - a purpose-written annotator that additionally *projects* the query
#          allele through the asm5 alignment and translates it, giving the true
#          protein-level consequence (frameshifts, premature stops, ...)
#  (3) Gene-content comparison between reference and query annotation.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

ANN="$OUT/04_annotation"
mkdir -p "$ANN/prokka_query" "$ANN/bakta_query" "$ANN/snpeff" \
         "$ANN/bcftools_csq" "$ANN/variant_effect" "$ANN/gene_content"

# Prefer the verified confident set; fall back to the merged high-confidence VCF
VCF="$OUT/03_variants/verified/variants_confident.vcf"
[ -s "$VCF" ] || VCF="$OUT/03_variants/merged/variants_highconf.vcf"
BAM_ASM="$WORK/bam/${QRY_NAME}_to_${REF_NAME}.asm5.bam"

# -----------------------------------------------------------------------------
msg "04-ANN (1a): Prokka annotation of the query assembly"
# -----------------------------------------------------------------------------
if command -v prokka >/dev/null 2>&1; then
  if [ ! -s "$ANN/prokka_query/${QRY_NAME}.gff" ]; then
    run_logged prokka_query \
      prokka --outdir "$ANN/prokka_query" --prefix "$QRY_NAME" \
             --cpus "$THREADS" --kingdom Bacteria --gcode 11 \
             --locustag "${QRY_NAME}_" --compliant --force \
             --addgenes "$QRY_FA" || warn "prokka failed"
  fi
  # also annotate the reference with the same tool so the two are comparable
  if [ ! -s "$ANN/prokka_reference/${REF_NAME}.gff" ]; then
    mkdir -p "$ANN/prokka_reference"
    run_logged prokka_reference \
      prokka --outdir "$ANN/prokka_reference" --prefix "$REF_NAME" \
             --cpus "$THREADS" --kingdom Bacteria --gcode 11 \
             --locustag "${REF_NAME}_" --compliant --force \
             --addgenes "$REF_FA" || warn "prokka reference failed"
  fi
else
  warn "prokka not installed"
fi

# -----------------------------------------------------------------------------
msg "04-ANN (1b): Bakta annotation of the query assembly (COG/KEGG/GO/EC)"
# -----------------------------------------------------------------------------
BAKTA_DB="$WORK/bakta_db"
if command -v bakta >/dev/null 2>&1; then
  if [ ! -d "$BAKTA_DB" ] || [ -z "$(ls -A "$BAKTA_DB" 2>/dev/null)" ]; then
    msg "  downloading the Bakta database (light) ..."
    run_logged bakta_download bash -c \
      "bakta_db download --output '$BAKTA_DB' --type light" || warn "bakta_db download failed"
  fi
  if [ -d "$BAKTA_DB" ] && [ -n "$(ls -A "$BAKTA_DB" 2>/dev/null)" ]; then
    if [ ! -s "$ANN/bakta_query/${QRY_NAME}.gff3" ]; then
      run_logged bakta_query \
        bakta --db "$BAKTA_DB" --output "$ANN/bakta_query" \
              --prefix "$QRY_NAME" --threads "$THREADS" \
              --force --compliant "$QRY_FA" || warn "bakta failed"
    fi
    if [ ! -s "$ANN/bakta_reference/${REF_NAME}.gff3" ]; then
      mkdir -p "$ANN/bakta_reference"
      run_logged bakta_reference \
        bakta --db "$BAKTA_DB" --output "$ANN/bakta_reference" \
              --prefix "$REF_NAME" --threads "$THREADS" \
              --force --compliant "$REF_FA" || warn "bakta reference failed"
    fi
  fi
else
  warn "bakta not installed"
fi

# -----------------------------------------------------------------------------
msg "04-ANN (1c): resolve which query annotation to use downstream"
# -----------------------------------------------------------------------------
# Prokka and Bakta are alternatives; whichever succeeded provides the query GFF
# and protein set used by the comparison steps.  Bakta is preferred when both
# are present because it carries COG/KEGG/GO/EC cross-references.
QRY_GFF=""
QRY_FAA_ANN=""
if [ -s "$ANN/bakta_query/${QRY_NAME}.gff3" ]; then
  QRY_GFF="$ANN/bakta_query/${QRY_NAME}.gff3"
  [ -s "$ANN/bakta_query/${QRY_NAME}.faa" ] && QRY_FAA_ANN="$ANN/bakta_query/${QRY_NAME}.faa"
fi
if [ ! -s "$QRY_GFF" ] && [ -s "$ANN/prokka_query/${QRY_NAME}.gff" ]; then
  QRY_GFF="$ANN/prokka_query/${QRY_NAME}.gff"
  [ -s "$ANN/prokka_query/${QRY_NAME}.faa" ] && QRY_FAA_ANN="$ANN/prokka_query/${QRY_NAME}.faa"
fi
if [ -z "$QRY_GFF" ] && [ -s "$ANN/transferred/query.gff3" ]; then
  QRY_GFF="$ANN/transferred/query.gff3"
  QRY_FAA_ANN="$ANN/transferred/query.faa"
fi
if [ -z "$QRY_GFF" ]; then
  warn "no query annotation available -- using the reference annotation"
  QRY_GFF="$REF_GFF"
fi
[ -n "$QRY_FAA_ANN" ] || QRY_FAA_ANN="$REF_FAA"
printf '%s\n' "$QRY_GFF" > "$ANN/query_annotation_source.txt"
echo "  query annotation : $QRY_GFF"
echo "  query proteins   : $QRY_FAA_ANN"
export QRY_GFF QRY_FAA_ANN

# -----------------------------------------------------------------------------
msg "04-ANN (2a): snpEff custom database + variant effect prediction"
# -----------------------------------------------------------------------------
if command -v snpEff >/dev/null 2>&1 && [ -s "$VCF" ]; then
  SE_DATA="$WORK/snpeff_data"
  SE_CFG="$WORK/snpeff.config"
  mkdir -p "$SE_DATA/$SNPEFF_DB"
  cp "$REF_FA"  "$SE_DATA/$SNPEFF_DB/sequences.fa"
  cp "$REF_GFF" "$SE_DATA/$SNPEFF_DB/genes.gff"
  cat > "$SE_CFG" <<EOF
data.dir = $SE_DATA
$SNPEFF_DB.genome : Bacteroides fragilis $REF_NAME
EOF
  # snpEff needs the sequence names to match the GFF exactly
  run_logged snpeff_build \
    snpEff build -gff3 -v -noCheckCds -noCheckProtein \
           -c "$SE_CFG" "$SNPEFF_DB" || warn "snpEff build failed"

  if [ -s "$SE_DATA/$SNPEFF_DB/snpEffectPredictor.bin" ]; then
    run_logged snpeff_annotate bash -c \
      "snpEff -c '$SE_CFG' -v -noStats -noLog '$SNPEFF_DB' '$VCF' \
         > '$ANN/snpeff/variants_snpeff.vcf' 2> '$LOGDIR/snpeff.err'"
    run_logged snpeff_genes \
      snpEff genes2bed -c "$SE_CFG" "$SNPEFF_DB" \
        > "$ANN/snpeff/genes.bed" || true
    # tabular summary of the effect classes snpEff found
    if [ -s "$ANN/snpeff/variants_snpeff.vcf" ]; then
      grep -v '^#' "$ANN/snpeff/variants_snpeff.vcf" | \
        awk -F'\t' '{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~/^ANN=/){split(a[i],b,"="); print b[2]}}' | \
        tr ',' '\n' | cut -d'|' -f2 | sort | uniq -c | sort -rn \
        > "$ANN/snpeff/snpeff_effect_counts.txt" || true
    fi
  fi
else
  warn "snpEff unavailable or no high-confidence VCF"
fi

# -----------------------------------------------------------------------------
msg "04-ANN (2b): bcftools csq (independent consequence caller)"
# -----------------------------------------------------------------------------
if [ -s "$VCF" ] && command -v bcftools >/dev/null 2>&1; then
  # bcftools csq wants phased genotypes; our haploid calls are already 1/1
  run_logged bcftools_csq bash -c \
    "bcftools csq -f '$REF_FA' -g '$REF_GFF' --phase a -Ov '$VCF' \
       > '$ANN/bcftools_csq/variants_csq.vcf' 2> '$LOGDIR/bcftools_csq.err'" || \
    warn "bcftools csq failed (often due to GFF/GFF3 attribute quirks)"
fi

# -----------------------------------------------------------------------------
msg "04-ANN (2c): allele-aware custom annotation (reference + query protein)"
# -----------------------------------------------------------------------------
HC="$OUT/03_variants/merged/variants_highconf.tsv"
[ -s "$HC" ] || HC="$OUT/03_variants/merged/variants_merged.tsv"

python3 "$SCRIPTS/lib/annotate_variants.py" \
  --variants "$HC" \
  --gff "$REF_GFF" \
  --ref "$REF_FA" \
  --bam "$BAM_ASM" \
  --outdir "$ANN/variant_effect"

# -----------------------------------------------------------------------------
msg "04-ANN (2d): full query annotation by alignment projection"
# -----------------------------------------------------------------------------
# The reference carries a curated RefSeq annotation.  At 99.99 % identity it
# transfers gene by gene through the alignment, which yields a complete query
# annotation (query coordinates + query sequence) without a de-novo gene caller.
python3 "$SCRIPTS/lib/transfer_annotation.py" \
  --ref-gff "$REF_GFF" --ref-fna "$REF_FA" --qry-fna "$QRY_FA" \
  --bam    "$WORK/bam/${QRY_NAME}_to_${REF_NAME}.asm5.bam" \
  --coords "$OUT/02_alignment/synteny/nucmer_coords_1to1.tsv" \
  --outdir "$ANN/transferred" || warn "annotation transfer failed"

# -----------------------------------------------------------------------------
msg "04-ANN (3): reference vs query gene-content comparison"
# -----------------------------------------------------------------------------
python3 "$SCRIPTS/lib/gene_content_compare.py" \
  --ref-gff "$REF_GFF" \
  --qry-gff "$QRY_GFF" \
  --ref-faa "$REF_FAA" \
  --qry-faa "$QRY_FAA_ANN" \
  --coords "$OUT/02_alignment/synteny/nucmer_coords_1to1.tsv" \
  --outdir "$ANN/gene_content" --threads "$THREADS" || \
  warn "gene-content comparison skipped"

msg "04-Annotation complete -> $ANN"
