#!/usr/bin/env bash
# #############################################################################
#                                                                             #
#   WGS.sh -- complete whole-genome comparative analysis pipeline             #
#                                                                             #
#   Query     : Bacteroides fragilis CCUG 4856   (GCF_005706655.1)            #
#   Reference : Bacteroides fragilis NCTC 9343   (GCF_000025985.1)            #
#                                                                             #
#   Everything runs inside the conda environment `WGS`.                       #
#                                                                             #
#   USAGE                                                                     #
#     bash WGS.sh                 # run every step, skipping finished ones     #
#     bash WGS.sh --force         # re-run everything from scratch             #
#     bash WGS.sh --steps 1,3,4   # run only the listed steps                  #
#     bash WGS.sh --list          # show the steps                             #
#     bash WGS.sh --check         # verify the toolchain, run nothing          #
#                                                                             #
#   STEPS                                                                     #
#     0  prepare      decompress inputs, build the directory tree              #
#     1  qc           assembly QC, QUAST, in-silico read QC                    #
#     2  alignment    NUCmer WGA, dnadiff, contig mapping, coverage            #
#     3  variants     four-way variant calling, reconciliation, statistics     #
#     4  annotation   Prokka/Bakta, snpEff/csq, allele-aware variant effects   #
#     5  functional   AMR / virulence / metabolism screening and comparison    #
#     6  extra        phylogenomics, pan-genome, prophage, CRISPR, IS, plasmid #
#     7  visualize    every figure of the report                               #
#     8  report       assemble the final Markdown/HTML report and tables       #
#                                                                             #
#   Results land in /home/cz/WGS/result, logs in /home/cz/WGS/logs.            #
# #############################################################################

set -uo pipefail

# --------------------------------------------------------------------------- #
#  locate ourselves and the conda environment
# --------------------------------------------------------------------------- #
WGS_ROOT="/home/cz/WGS"
cd "$WGS_ROOT" || { echo "cannot cd to $WGS_ROOT"; exit 1; }

CONDA_SH="/home/cz/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="${CONDA_ENV:-WGS}"

if [ -f "$CONDA_SH" ]; then
  # shellcheck disable=SC1090
  source "$CONDA_SH"
  conda activate "$CONDA_ENV" 2>/dev/null || {
    echo "ERROR: conda environment '$CONDA_ENV' not found."; exit 1; }
fi

export WGS_ROOT
mkdir -p "$WGS_ROOT/logs" "$WGS_ROOT/result" "$WGS_ROOT/work" "$WGS_ROOT/scripts"

PIPELINE_LOG="$WGS_ROOT/logs/WGS_pipeline.log"

# --------------------------------------------------------------------------- #
#  arguments
# --------------------------------------------------------------------------- #
FORCE=0
STEPS_ARG=""
SHOW_LIST=0
CHECK_ONLY=0

while [ $# -gt 0 ]; do
  case "$1" in
    --force|-f) FORCE=1 ;;
    --steps)    STEPS_ARG="$2"; shift ;;
    --list|-l)  SHOW_LIST=1 ;;
    --check|-c) CHECK_ONLY=1 ;;
    -h|--help)  sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "unknown option: $1"; exit 2 ;;
  esac
  shift
done

STEPS=(prepare qc alignment variants annotation functional extra visualize report)

if [ "$SHOW_LIST" = "1" ]; then
  i=0
  for s in "${STEPS[@]}"; do echo "  $i  $s"; i=$((i+1)); done
  exit 0
fi

if [ -n "$STEPS_ARG" ]; then
  WANTED=",$(echo "$STEPS_ARG" | tr -d ' '),"
else
  WANTED=","
  i=0
  for s in "${STEPS[@]}"; do WANTED="$WANTED$i,"; i=$((i+1)); done
fi

want_step() {   # want_step <index>
  case "$WANTED" in *",$1,"*) return 0 ;; *) return 1 ;; esac
}

# --------------------------------------------------------------------------- #
#  toolchain verification
# --------------------------------------------------------------------------- #
CORE_TOOLS=(samtools bcftools bwa minimap2 nucmer show-snps python3)
OPT_TOOLS=(freebayes dnadiff prokka bakta snpEff abricate amrfinder quast.py \
           seqkit bandage mafft iqtree diamond Rscript)

check_tools() {
  echo "------------------------------------------------------------------"
  echo " toolchain check  (conda env: ${CONDA_NAME:-$CONDA_ENV})"
  echo "------------------------------------------------------------------"
  local missing_core=0
  for t in "${CORE_TOOLS[@]}"; do
    if command -v "$t" >/dev/null 2>&1; then
      printf '  [ ok ]  %-14s %s\n' "$t" "$(command -v "$t")"
    else
      printf '  [MISS]  %-14s REQUIRED\n' "$t"; missing_core=$((missing_core+1))
    fi
  done
  for t in "${OPT_TOOLS[@]}"; do
    if command -v "$t" >/dev/null 2>&1; then
      printf '  [ ok ]  %-14s\n' "$t"
    else
      printf '  [ -- ]  %-14s optional, step will be skipped\n' "$t"
    fi
  done
  echo "------------------------------------------------------------------"
  [ "$missing_core" -eq 0 ] && echo " core toolchain complete" \
                           || echo " WARNING: $missing_core required tool(s) missing"
  return "$missing_core"
}

if [ "$CHECK_ONLY" = "1" ]; then
  check_tools; exit 0
fi

# --------------------------------------------------------------------------- #
#  step runner
# --------------------------------------------------------------------------- #
declare -A STEP_STATUS=()
declare -A STEP_TIME=()

run_step() {
  local idx="$1" name="$2" script="$3"
  want_step "$idx" || { echo "-- skipping step $idx ($name)"; return 0; }

  echo ""
  echo "###################################################################"
  echo "#  STEP $idx  :  $name"
  echo "#  script    :  $script"
  echo "#  started   :  $(date '+%Y-%m-%d %H:%M:%S')"
  echo "###################################################################"

  if [ ! -f "$script" ]; then
    echo "!! script not found: $script"
    STEP_STATUS[$name]="MISSING"; return 1
  fi

  local t0 t1 rc
  t0=$(date +%s)
  # each step logs itself; we also tee the console
  bash "$script" 2>&1 | tee -a "$PIPELINE_LOG.$name"
  rc=${PIPESTATUS[0]}
  t1=$(date +%s)
  STEP_TIME[$name]=$((t1 - t0))

  if [ "$rc" -eq 0 ]; then
    STEP_STATUS[$name]="OK"
    echo "--> STEP $idx ($name) finished OK in ${STEP_TIME[$name]}s"
  else
    STEP_STATUS[$name]="FAILED(rc=$rc)"
    echo "--> STEP $idx ($name) FAILED (rc=$rc) after ${STEP_TIME[$name]}s"
    echo "    continuing with the remaining steps; see logs/${name}.log"
  fi
  return 0
}

# --------------------------------------------------------------------------- #
#  step 0 : prepare
# --------------------------------------------------------------------------- #
step_prepare() {
  set -e
  mkdir -p "$WGS_ROOT"/{logs,scripts,work,result}
  mkdir -p "$WGS_ROOT"/result/{01_qc,02_alignment,03_variants,04_annotation,05_functional,06_extra,07_visualization,08_report}
  mkdir -p "$WGS_ROOT"/work/{ref,qry,bam,vcf,reads}
  mkdir -p "$WGS_ROOT"/reference/unzipped "$WGS_ROOT"/seq/unzipped

  if [ ! -s "$WGS_ROOT/reference/unzipped/reference.fna" ]; then
    echo "  decompressing the reference ..."
    zcat "$WGS_ROOT"/reference/*.fna.gz > "$WGS_ROOT/reference/unzipped/reference.fna"
  fi
  if [ ! -s "$WGS_ROOT/reference/unzipped/reference.gff" ]; then
    zcat "$WGS_ROOT"/reference/*.gff.gz > "$WGS_ROOT/reference/unzipped/reference.gff"
  fi
  if [ ! -s "$WGS_ROOT/reference/unzipped/reference.faa" ]; then
    zcat "$WGS_ROOT"/reference/*.faa.gz > "$WGS_ROOT/reference/unzipped/reference.faa"
  fi
  if [ ! -s "$WGS_ROOT/seq/unzipped/query.fna" ]; then
    echo "  decompressing the query ..."
    zcat "$WGS_ROOT"/seq/*.fna.gz > "$WGS_ROOT/seq/unzipped/query.fna"
  fi

  samtools faidx "$WGS_ROOT/reference/unzipped/reference.fna"
  echo "  inputs ready:"
  ls -l "$WGS_ROOT/reference/unzipped/reference.fna" \
        "$WGS_ROOT/reference/unzipped/reference.gff" \
        "$WGS_ROOT/seq/unzipped/query.fna"
}

# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #
{
echo "==================================================================="
echo " WGS comparative genomics pipeline"
echo " query     : $WGS_ROOT/seq/unzipped/query.fna"
echo " reference : $WGS_ROOT/reference/unzipped/reference.fna"
echo " started   : $(date '+%Y-%m-%d %H:%M:%S')"
echo " host      : $(hostname), $(nproc) cores"
echo "==================================================================="
} | tee -a "$PIPELINE_LOG"

check_tools || true

T_START=$(date +%s)

run_step 0 prepare    "$WGS_ROOT/scripts/01_prepare.sh"
run_step 1 qc         "$WGS_ROOT/scripts/01_qc.sh"
run_step 2 alignment  "$WGS_ROOT/scripts/02_alignment.sh"
run_step 3 variants   "$WGS_ROOT/scripts/03_variants.sh"
run_step 4 annotation "$WGS_ROOT/scripts/04_annotate.sh"
run_step 5 functional "$WGS_ROOT/scripts/05_functional.sh"
run_step 6 extra      "$WGS_ROOT/scripts/06_extra.sh"
run_step 7 visualize  "$WGS_ROOT/scripts/07_visualize.sh"
run_step 8 report     "$WGS_ROOT/scripts/08_report.sh"

T_END=$(date +%s)

{
echo ""
echo "==================================================================="
echo " PIPELINE SUMMARY"
echo "==================================================================="
printf ' %-12s %-16s %s\n' "STEP" "STATUS" "SECONDS"
for s in "${STEPS[@]}"; do
  [ "$s" = "prepare" ] && continue
  printf ' %-12s %-16s %s\n' "$s" "${STEP_STATUS[$s]:-not run}" "${STEP_TIME[$s]:-}"
done
echo "-------------------------------------------------------------------"
echo " total wall time: $((T_END - T_START)) s"
echo " results        : $WGS_ROOT/result"
echo " report         : $WGS_ROOT/result/08_report/WGS_report.html"
echo " finished       : $(date '+%Y-%m-%d %H:%M:%S')"
echo "==================================================================="
} | tee -a "$PIPELINE_LOG"
