#!/usr/bin/env bash
# =============================================================================
#  config.sh -- shared configuration for the WGS comparative-genomics pipeline
# =============================================================================
#  Project : Whole-genome comparison of Bacteroides fragilis CCUG 4856 (query)
#            against B. fragilis NCTC 9343 (reference)
#  Author  : automated pipeline
# =============================================================================

# ---------------------------------------------------------------- directories
export WGS_ROOT="/home/cz/WGS"
export REF_DIR="$WGS_ROOT/reference"
export SEQ_DIR="$WGS_ROOT/seq"
export OUT="$WGS_ROOT/result"
export WORK="$WGS_ROOT/work"
export LOGDIR="$WGS_ROOT/logs"
export SCRIPTS="$WGS_ROOT/scripts"

# ------------------------------------------------------- input files (unzipped)
export REF_FA="$REF_DIR/unzipped/reference.fna"
export REF_GFF="$REF_DIR/unzipped/reference.gff"
export REF_FAA="$REF_DIR/unzipped/reference.faa"
export QRY_FA="$SEQ_DIR/unzipped/query.fna"

# ------------------------------------------------------------- sample naming
export REF_NAME="NCTC9343"          # reference strain short name
export QRY_NAME="CCUG4856"          # query strain short name
export REF_ACC="GCF_000025985.1"
export QRY_ACC="GCF_005706655.1"
export REF_TAXID="272559"           # B. fragilis NCTC 9343 taxid
export SPECIES="Bacteroides_fragilis"
export SNPEFF_DB="Bfragilis_${REF_NAME}"

# --------------------------------------------------------------- compute
export THREADS="${THREADS:-32}"      # main thread count
export MAXTHREADS="${MAXTHREADS:-64}" # for the heaviest steps
export MEM_MB="${MEM_MB:-64000}"

# ----------------------------------------------------- simulated-read defaults
# The query is a *finished assembly*.  To provide a genuine read-level QC and a
# depth-based confirmation of the assembly-derived variant calls, the pipeline
# generates in-silico Illumina reads from the query assembly.  These reads are
# clearly labelled `SIM` everywhere so they are never mistaken for wet-lab data.
export SIM_COVERAGE=100
export SIM_READLEN=150
export SIM_INSERTSIZE=350

# ------------------------------------------------------------- helper output
msg()  { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33m[WARN] %s\033[0m\n' "$*" >&2; }
die()  { printf '\033[1;31m[ERROR] %s\033[0m\n' "$*" >&2; exit 1; }

# ------------------------------------------------------------------- logging
# run_logged <logname> <cmd...>
#
# Call sites use the natural idiom
#     run_logged <name> <cmd> [args] > <result-file>
# so the helper must NOT touch the command's stdout: it must reach the caller's
# redirection untouched.  Only stderr is captured -- it is both appended to
# logs/<name>.log and echoed to the terminal.  The status line goes to stderr so
# it can never contaminate a result file.
run_logged() {
  local name="$1"; shift
  mkdir -p "$LOGDIR"
  : > "$LOGDIR/$name.log"
  local rc=0
  "$@" 2> >(tee -a "$LOGDIR/$name.log" >&2) || rc=$?
  if [ "$rc" -eq 0 ]; then
    echo "   [ok]   $name" >&2
  else
    echo "   [FAIL] $name (rc=$rc) -- see $LOGDIR/$name.log" >&2
    return 1
  fi
}

# Check that a program exists in the active environment
need() {
  for p in "$@"; do
    command -v "$p" >/dev/null 2>&1 || die "required program not found: $p"
  done
}

# ------------------------------------------------------------ housekeeping
mkdir -p "$OUT" "$WORK" "$LOGDIR" \
         "$WORK/ref" "$WORK/qry" "$WORK/bam" "$WORK/vcf" "$WORK/reads"

# Reference sequence names/accessions useful downstream
export REF_CHR="NC_003228.3"        # 5,205,140 bp chromosome
export REF_PLASMID="NC_006873.1"    # 36,560 bp plasmid
