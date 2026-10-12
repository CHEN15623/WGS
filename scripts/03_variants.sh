#!/usr/bin/env bash
# =============================================================================
#  03_variants.sh -- VARIANT DETECTION AND STATISTICS
# =============================================================================
#  Variants are called by four independent routes and then reconciled:
#    1. MUMmer `show-snps`   on the unique-anchor NUCmer delta   (assembly)
#    2. `bcftools call`      on the minimap2 asm5 contig BAM      (assembly)
#    3. `freebayes`          on the same BAM                      (assembly)
#    4. `bcftools call`      on the SIM-read BAM (deep coverage)  (reads)
#  The merged high-confidence set requires agreement of >=2 independent routes,
#  and every variant is annotated with the mapping quality observed in the
#  contig alignment so that calls inside repeats can be recognised.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

VAR="$OUT/03_variants"
mkdir -p "$VAR/mummer" "$VAR/bcftools_asm" "$VAR/freebayes_asm" \
         "$VAR/bcftools_sim" "$VAR/stats" "$VAR/merged"

need samtools bcftools show-snps
[ -s "$REF_FA.fai" ] || samtools faidx "$REF_FA"

BAM_ASM="$WORK/bam/${QRY_NAME}_to_${REF_NAME}.asm5.bam"
BAM_SIM="$WORK/bam/${QRY_NAME}_SIM_to_${REF_NAME}.bam"
NUC="$WORK/nucmer"

# =============================================================================
msg "03-VAR (1): MUMmer show-snps (assembly-to-assembly)"
# =============================================================================
if [ -s "${NUC}.1to1.delta" ]; then
  # IMPORTANT: `delta-filter` resolves every reference position to exactly one
  # query position, so on a *filtered* delta the alignment-overlap columns
  # [R]/[Q] that show-snps uses for `-C` are all zero -- `-C` would be far too
  # strict (it would drop almost everything).  The 1-to-1 delta is therefore the
  # correct input, and three views are produced:
  #   .1to1.tsv           every difference along the 1-to-1 path, [R]/[Q] kept
  #   .unique.tsv         strict unique-alignment subset (-C)
  #   .repeatfiltered.tsv -C applied to the RAW delta (finds real repeats)
  # Columns of `show-snps -rlTH`:
  #   P1 REF QRY P2 BUFF DIST R Q LENR LENQ FRM ? TAG1 TAG2
  run_logged show_snps \
    show-snps -rlTH "${NUC}.1to1.delta" > "$VAR/mummer/mummer_snps.1to1.tsv"
  run_logged show_snps_strict \
    show-snps -rlTHC "${NUC}.1to1.delta" > "$VAR/mummer/mummer_snps.unique.tsv"
  run_logged show_snps_repeatfiltered \
    show-snps -rlTHC "${NUC}.delta" > "$VAR/mummer/mummer_snps.repeatfiltered.tsv"
  # indels / differences reported by show-diff
  # show-diff categorises structural differences as GAP/DUP/BRK/JMP/INV/SEQ;
  # it only accepts -f -h -H -q -r (no -T/-l).
  run_logged show_diff \
    show-diff -rH "${NUC}.1to1.delta" > "$VAR/mummer/mummer_structural_diff.txt" || true
fi

# =============================================================================
msg "03-VAR (2): bcftools + freebayes on the contig (asm5) alignment"
# =============================================================================
if [ -s "$BAM_ASM" ]; then
  [ -s "$BAM_ASM.bai" ] || samtools index -@ 8 "$BAM_ASM"

  # ---- bcftools ---------------------------------------------------------
  if [ ! -s "$VAR/bcftools_asm/bcftools_asm.vcf.gz" ]; then
    run_logged bcftools_call_asm bash -c \
      "bcftools mpileup -f '$REF_FA' --no-BAQ -d 100000 -a FORMAT/AD,FORMAT/DP,INFO/AD \
          -Ou '$BAM_ASM' \
       | bcftools call -m -v --ploidy 1 -Ou \
       | bcftools norm -f '$REF_FA' -m -both -Oz -o '$VAR/bcftools_asm/bcftools_asm.vcf.gz'"
  fi
  [ -s "$VAR/bcftools_asm/bcftools_asm.vcf.gz" ] && \
    bcftools index -f -t "$VAR/bcftools_asm/bcftools_asm.vcf.gz"

  # ---- freebayes --------------------------------------------------------
  if [ ! -s "$VAR/freebayes_asm/freebayes_asm.vcf.gz" ]; then
    # NOTE: a contig alignment has no meaningful base qualities (minimap2 writes
    # '*'), so --min-base-quality must be 0 here or every observation is dropped.
    run_logged freebayes_asm bash -c \
      "freebayes -f '$REF_FA' -p 1 --min-mapping-quality 20 --min-base-quality 0 \
          --min-alternate-count 1 --min-alternate-fraction 0.5 --no-population-priors \
          '$BAM_ASM' 2> '$LOGDIR/freebayes_asm.err' \
       | bgzip -c > '$VAR/freebayes_asm/freebayes_asm.vcf.gz'" || \
    run_logged freebayes_asm_retry bash -c \
      "freebayes -f '$REF_FA' -p 1 --min-alternate-fraction 0.5 \
          '$BAM_ASM' 2> '$LOGDIR/freebayes_asm.err' \
       | bgzip -c > '$VAR/freebayes_asm/freebayes_asm.vcf.gz'" || true
  fi
fi

# =============================================================================
msg "03-VAR (3): read-based calling on the deep SIM-read alignment"
# =============================================================================
# Three read-based routes are produced:
#   bcftools_sim_snp    bcftools with -I (indel alleles suppressed)
#   bcftools_sim_indel  bcftools restricted to indels
#   freebayes_sim       freebayes
# Splitting the bcftools pass is not cosmetic: with deep read data the
# multiallelic caller can latch onto the `<*>` (spanning-deletion) allele and
# then fail to emit an obvious homozygous SNP.  Suppressing indel alleles with
# `-I` in the SNP pass removes that failure mode entirely.
SIMSNP="$VAR/bcftools_sim/bcftools_sim.snp.vcf.gz"
SIMIND="$VAR/bcftools_sim/bcftools_sim.indel.vcf.gz"
SIMFB="$VAR/freebayes_sim/freebayes_sim.vcf.gz"
mkdir -p "$VAR/freebayes_sim"

if [ -s "$BAM_SIM" ]; then
  [ -s "$BAM_SIM.bai" ] || samtools index -@ 8 "$BAM_SIM"

  # ---- bcftools, SNP pass (-I) ------------------------------------------
  if [ ! -s "$SIMSNP" ]; then
    run_logged bcftools_sim_snp bash -c \
      "bcftools mpileup -f '$REF_FA' -q 20 -Q 13 -I -d 100000 \
          -a FORMAT/AD,FORMAT/DP,INFO/AD -Ou '$BAM_SIM' \
       | bcftools call -m -v --ploidy 1 -Ou \
       | bcftools norm -f '$REF_FA' -m -both -Oz -o '$SIMSNP'"
  fi
  [ -s "$SIMSNP" ] && bcftools index -f -t "$SIMSNP"

  # ---- bcftools, indel pass ---------------------------------------------
  if [ ! -s "$SIMIND" ]; then
    run_logged bcftools_sim_indel bash -c \
      "bcftools mpileup -f '$REF_FA' -q 20 -Q 13 -d 100000 \
          -a FORMAT/AD,FORMAT/DP,INFO/AD -Ou '$BAM_SIM' \
       | bcftools call -m -v --ploidy 1 -Ou \
       | bcftools norm -f '$REF_FA' -m -both -Ou \
       | bcftools view -v indels -Oz -o '$SIMIND'"
  fi
  [ -s "$SIMIND" ] && bcftools index -f -t "$SIMIND"

  # ---- freebayes --------------------------------------------------------
  if [ ! -s "$SIMFB" ] && command -v freebayes >/dev/null 2>&1; then
    run_logged freebayes_sim bash -c \
      "freebayes -f '$REF_FA' -p 1 --min-mapping-quality 20 --min-base-quality 13 \
          --min-alternate-count 3 --min-alternate-fraction 0.7 \
          '$BAM_SIM' 2> '$LOGDIR/freebayes_sim.err' \
       | bgzip -c > '$SIMFB'" || warn "freebayes on SIM reads failed"
  fi
  [ -s "$SIMFB" ] && bcftools index -f -t "$SIMFB" 2>/dev/null || true

  # ---- a filtered, high-confidence read-derived SNP set -----------------
  if [ -s "$SIMSNP" ]; then
    run_logged bcftools_filter_sim bash -c \
      "bcftools filter -e 'QUAL<30 || FORMAT/DP<10' -s LowQual '$SIMSNP' \
         | bcftools view -f PASS -Oz -o '$VAR/bcftools_sim/bcftools_sim.highconf.vcf.gz'"
    [ -s "$VAR/bcftools_sim/bcftools_sim.highconf.vcf.gz" ] && \
      bcftools index -f -t "$VAR/bcftools_sim/bcftools_sim.highconf.vcf.gz"
  fi
fi

# =============================================================================
msg "03-VAR (4): merge callers into one reconciled variant table"
# =============================================================================
python3 "$SCRIPTS/lib/variant_merge.py" \
  --mummer      "$VAR/mummer/mummer_snps.1to1.tsv" \
  --asm-vcf     "bcftools_asm=$VAR/bcftools_asm/bcftools_asm.vcf.gz" \
  --asm-vcf     "freebayes_asm=$VAR/freebayes_asm/freebayes_asm.vcf.gz" \
  --read-vcf    "bcftools_sim_snp=$SIMSNP" \
  --read-vcf    "bcftools_sim_indel=$SIMIND" \
  --read-vcf    "freebayes_sim=$SIMFB" \
  --coords-all  "$OUT/02_alignment/synteny/nucmer_coords_all.tsv" \
  --bam-asm     "$BAM_ASM" \
  --bam-sim     "$BAM_SIM" \
  --ref         "$REF_FA" \
  --outdir      "$VAR/merged" \
  --threads     "$THREADS"

# =============================================================================
msg "03-VAR (4b): verification -- reconstruct ref/query alleles, split confident vs ambiguous"
# =============================================================================
mkdir -p "$VAR/verified"
python3 "$SCRIPTS/lib/variant_verification.py" \
  --merged "$VAR/merged/variants_merged.tsv" \
  --bam    "$BAM_ASM" \
  --ref    "$REF_FA" \
  --outdir "$VAR/verified" || warn "variant verification failed"

# =============================================================================
msg "03-VAR (5): variant statistics"
# =============================================================================
for v in "$VAR/bcftools_asm/bcftools_asm.vcf.gz" \
         "$VAR/freebayes_asm/freebayes_asm.vcf.gz" \
         "$SIMSNP" "$SIMIND" "$SIMFB"; do
  [ -s "$v" ] || continue
  base=$(basename "$v" .vcf.gz)
  run_logged "stats_$base" bash -c "bcftools stats '$v' > '$VAR/stats/${base}.bcftools_stats.txt'" || true
done

python3 "$SCRIPTS/lib/variant_stats.py" \
  --merged "$VAR/merged/variants_merged.tsv" \
  --confident "$VAR/verified/variants_confident.tsv" \
  --ref    "$REF_FA" \
  --gff    "$REF_GFF" \
  --outdir "$VAR/stats"

msg "03-Variants complete -> $VAR"
