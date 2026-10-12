#!/usr/bin/env bash
# =============================================================================
#  02_alignment.sh -- GENOME ALIGNMENT
# =============================================================================
#  Two complementary alignment strategies:
#    (A) ASSEMBLY-TO-ASSEMBLY whole-genome alignment with NUCmer/MUMmer4
#        -> delta, 1-to-1 filtered delta, coordinate blocks, dnadiff report.
#        This is the primary evidence layer for synteny and for variant calling.
#    (B) CONTIG-TO-REFERENCE read-style alignment with minimap2 -x asm5
#        -> sorted/indexed BAM, flagstat/stats/idxstats, depth profiles.
#        Gives per-base coverage and lets standard VCF tools be used.
#    (C) SIM-read-to-reference alignment (bwa mem + minimap2 sr) -> BAM used to
#        confirm the assembly-derived variants at high depth.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

ALN="$OUT/02_alignment"
mkdir -p "$ALN/wga" "$ALN/dnadiff" "$ALN/contig_mapping" \
         "$ALN/sim_read_mapping" "$ALN/coverage" "$ALN/synteny"

need samtools minimap2 nucmer

# -----------------------------------------------------------------------------
msg "02-ALN (0): index the reference"
# -----------------------------------------------------------------------------
[ -s "$REF_FA.fai" ] || samtools faidx "$REF_FA"
samtools faidx "$REF_FA"
if [ ! -s "$REF_FA.bwt" ] && command -v bwa >/dev/null 2>&1; then
  bwa index "$REF_FA" > "$LOGDIR/bwa_index.log" 2>&1
fi

# =============================================================================
msg "02-ALN (A): whole-genome alignment with NUCmer"
# =============================================================================
cd "$WORK"
NUC="$WORK/nucmer"
if [ ! -s "${NUC}.delta" ]; then
  # --maxmatch: exhaustive anchoring, best for closely related genomes
  run_logged nucmer \
    nucmer --maxmatch -c 100 -l 40 -t "$THREADS" -p "$NUC" "$REF_FA" "$QRY_FA"
fi

if [ -s "${NUC}.delta" ]; then
  # 1-to-1 alignment (best mutual hits) -- used for synteny blocks
  [ -s "${NUC}.1to1.delta" ] || run_logged delta_filter_1to1 \
    delta-filter -1 "${NUC}.delta" > "${NUC}.1to1.delta"
  # alignment blocks that are uniquely aligned -- used for variant calling
  [ -s "${NUC}.uniq.delta" ] || run_logged delta_filter_uniq \
    delta-filter -q -r "${NUC}.delta" > "${NUC}.uniq.delta"

  run_logged show_coords \
    show-coords -rclTHd -I 70 "${NUC}.1to1.delta" > "$ALN/synteny/nucmer_coords_1to1.tsv"
  run_logged show_coords_all \
    show-coords -rclTHd -I 70 "${NUC}.delta" > "$ALN/synteny/nucmer_coords_all.tsv"
  cp "${NUC}.delta" "$ALN/wga/reference_vs_query.delta"
  cp "${NUC}.1to1.delta" "$ALN/wga/reference_vs_query.1to1.delta"

  # alignment-level statistics (full pairwise alignment of the chromosome)
  QRY_CHR=$(awk '/^>/{print $1; exit}' "$QRY_FA" | cut -c2-)
  run_logged show_align bash -c \
    "show-aligns -r -w 60 '${NUC}.1to1.delta' '$REF_CHR' '$QRY_CHR' \
       | gzip -c > '$ALN/wga/example_alignment_chromosome.txt.gz'" || true
fi

# =============================================================================
msg "02-ALN (A2): structural rearrangements from the 1-to-1 blocks"
# =============================================================================
python3 "$SCRIPTS/lib/structural_variants.py" \
  --coords     "$ALN/synteny/nucmer_coords_1to1.tsv" \
  --coords-all "$ALN/synteny/nucmer_coords_all.tsv" \
  --ref "$REF_FA" --qry "$QRY_FA" --gff "$REF_GFF" \
  --outdir "$ALN/synteny" || warn "structural variant analysis failed"

# =============================================================================
msg "02-ALN (B): dnadiff whole-genome difference report"
# =============================================================================
if command -v dnadiff >/dev/null 2>&1; then
  mkdir -p "$WORK/dnadiff"
  [ -s "$WORK/dnadiff/ref_vs_qry.report" ] || run_logged dnadiff \
    dnadiff -p "$WORK/dnadiff/ref_vs_qry" "$REF_FA" "$QRY_FA" || true
  cp -r "$WORK/dnadiff/." "$ALN/dnadiff/" 2>/dev/null || true
  [ -f "$WORK/dnadiff/ref_vs_qry.report" ] && \
    cp "$WORK/dnadiff/ref_vs_qry.report" "$ALN/dnadiff/dnadiff.report" 2>/dev/null || true
fi

# =============================================================================
msg "02-ALN (C): contig-to-reference mapping (minimap2 asm5)"
# =============================================================================
CM="$ALN/contig_mapping"
BAM_ASM="$WORK/bam/${QRY_NAME}_to_${REF_NAME}.asm5.bam"

if [ ! -s "$BAM_ASM" ]; then
  run_logged minimap2_asm5 bash -c \
    "minimap2 -ax asm5 --eqx -t $THREADS '$REF_FA' '$QRY_FA' \
       | samtools sort -@ 8 -m 2G -o '$BAM_ASM' -"
fi
[ -s "$BAM_ASM.bai" ] || samtools index -@ 8 "$BAM_ASM"

run_logged samtools_flagstat   samtools flagstat -@ 8 "$BAM_ASM" > "$CM/flagstat.txt"
# NOTE: `samtools stats` emits per-cycle tables indexed by read length.  Here a
# "read" is an entire contig (5.2 Mb), which would produce a multi-gigabyte file
# of meaningless per-cycle rows.  Only the summary (SN) and per-sequence (CH)
# sections are kept.
run_logged samtools_stats bash -c \
  "samtools stats -@ 8 '$BAM_ASM' | grep -E '^#|^SN|^CH' > '$CM/samtools_stats.txt'"
run_logged samtools_idxstats   samtools idxstats -@ 8 "$BAM_ASM" > "$CM/idxstats.tsv"
run_logged samtools_coverage   samtools coverage     "$BAM_ASM" > "$CM/coverage_summary.tsv"
# per-base depth over the whole reference (all positions, incl. zero coverage)
run_logged samtools_depth      bash -c \
  "samtools depth -a -q 0 -Q 0 '$BAM_ASM' | gzip -c \
     > '$ALN/coverage/contig_depth_asm5.tsv.gz'"

# =============================================================================
msg "02-ALN (D): simulated-read mapping (depth-supported confirmation layer)"
# =============================================================================
SR="$ALN/sim_read_mapping"
R1="$WORK/reads/${QRY_NAME}_SIM_R1.fastq"
R2="$WORK/reads/${QRY_NAME}_SIM_R2.fastq"
BAM_SIM="$WORK/bam/${QRY_NAME}_SIM_to_${REF_NAME}.bam"

if [ -s "$R1" ] && [ ! -s "$BAM_SIM" ] && command -v bwa >/dev/null 2>&1; then
  run_logged bwa_mem_sim bash -c \
    "bwa mem -t $THREADS -R '@RG\tID:SIM\tSM:${QRY_NAME}_SIM\tPL:ILLUMINA' \
       '$REF_FA' '$R1' '$R2' 2> '$LOGDIR/bwa_mem_sim.err' \
       | samtools sort -@ 8 -m 2G -o '$BAM_SIM' -"
fi
if [ -s "$BAM_SIM" ]; then
  [ -s "$BAM_SIM.bai" ] || samtools index -@ 8 "$BAM_SIM"
  run_logged samtools_flagstat_sim samtools flagstat -@ 8 "$BAM_SIM" > "$SR/flagstat.txt"
  run_logged samtools_stats_sim bash -c \
    "samtools stats -@ 8 '$BAM_SIM' | grep -E '^#|^SN|^CH' > '$SR/samtools_stats.txt'"
  run_logged samtools_coverage_sim samtools coverage     "$BAM_SIM" > "$SR/coverage_summary.tsv"
  # mark duplicates so depth-based variant calling can be duplicate-aware
  if ! [ -s "$WORK/bam/${QRY_NAME}_SIM.markdup.bam" ]; then
    run_logged markdup bash -c \
      "samtools sort -n -@ 8 -m 2G -o /tmp/sim.nsrt.bam '$BAM_SIM' && \
       samtools fixmate -m -@ 8 /tmp/sim.nsrt.bam /tmp/sim.fix.bam && \
       samtools sort -@ 8 -m 2G -o /tmp/sim.pos.bam /tmp/sim.fix.bam && \
       samtools markdup -@ 8 -s /tmp/sim.pos.bam '$WORK/bam/${QRY_NAME}_SIM.markdup.bam'" || true
    [ -s "$WORK/bam/${QRY_NAME}_SIM.markdup.bam" ] && \
      samtools index -@ 8 "$WORK/bam/${QRY_NAME}_SIM.markdup.bam"
  fi
  run_logged samtools_depth_sim bash -c \
    "samtools depth -a -q 0 -Q 0 '$BAM_SIM' | gzip -c \
       > '$ALN/coverage/sim_read_depth.tsv.gz'"
fi

# =============================================================================
msg "02-ALN (E): coverage / uncovered-region summary"
# =============================================================================
python3 "$SCRIPTS/lib/coverage_summary.py" \
  "$ALN/coverage/contig_depth_asm5.tsv.gz" \
  "$ALN/coverage/coverage_report.tsv" || true

msg "02-Alignment complete -> $ALN"
