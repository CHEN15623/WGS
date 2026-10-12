#!/usr/bin/env bash
# =============================================================================
#  01_prepare.sh -- step 0 : decompress inputs, build the workspace tree
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

msg "00-PREP: preparing inputs and workspace"

mkdir -p "$OUT"/{01_qc,02_alignment,03_variants,04_annotation,05_functional,06_extra,07_visualization,08_report}
mkdir -p "$WORK"/{ref,qry,bam,vcf,reads,cog_cache}
mkdir -p "$REF_DIR/unzipped" "$SEQ_DIR/unzipped" "$LOGDIR"

# ---- reference ------------------------------------------------------------
if [ ! -s "$REF_FA" ]; then
  msg "  decompressing the reference FASTA"
  zcat "$REF_DIR"/*.fna.gz > "$REF_FA"
fi
if [ ! -s "$REF_GFF" ]; then
  msg "  decompressing the reference GFF"
  zcat "$REF_DIR"/*.gff.gz > "$REF_GFF"
fi
if [ ! -s "$REF_FAA" ]; then
  msg "  decompressing the reference proteins"
  zcat "$REF_DIR"/*.faa.gz > "$REF_FAA"
fi

# ---- query ----------------------------------------------------------------
if [ ! -s "$QRY_FA" ]; then
  msg "  decompressing the query FASTA"
  zcat "$SEQ_DIR"/*.fna.gz > "$QRY_FA"
fi

# ---- sanity ---------------------------------------------------------------
for f in "$REF_FA" "$REF_GFF" "$REF_FAA" "$QRY_FA"; do
  [ -s "$f" ] || die "missing input: $f"
done

command -v samtools >/dev/null 2>&1 && samtools faidx "$REF_FA"
command -v samtools >/dev/null 2>&1 && samtools faidx "$QRY_FA"

# ---- record the environment for reproducibility ---------------------------
{
  echo "run_date: $(date '+%Y-%m-%d %H:%M:%S')"
  echo "host: $(hostname)"
  echo "cores: $(nproc)"
  echo "conda_env: ${CONDA_DEFAULT_ENV:-unknown}"
  echo "conda_prefix: ${CONDA_PREFIX:-unknown}"
  echo "--- tool versions ---"
  for t in samtools bcftools minimap2 bwa nucmer show-snps freebayes dnadiff \
           prokka bakta snpEff abricate amrfinder quast.py seqkit Rscript python3; do
    if command -v "$t" >/dev/null 2>&1; then
      v=$( { "$t" --version 2>&1 || "$t" -version 2>&1 || true; } | head -1 )
      echo "$t: $v"
    fi
  done
  echo "--- input checksums (md5) ---"
  md5sum "$REF_FA" "$REF_GFF" "$REF_FAA" "$QRY_FA"
} > "$OUT/08_report/software_versions.txt" 2>/dev/null || true

msg "00-PREP complete"
echo "  reference: $(grep -c '^>' "$REF_FA") sequence(s), $(grep -v '^>' "$REF_FA" | tr -d '\n' | wc -c) bp"
echo "  query    : $(grep -c '^>' "$QRY_FA") sequence(s), $(grep -v '^>' "$QRY_FA" | tr -d '\n' | wc -c) bp"
