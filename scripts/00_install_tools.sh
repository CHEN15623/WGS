#!/usr/bin/env bash
# Install the WGS toolchain into conda env `WGS`, one group at a time so that a
# single bad package cannot abort the whole solve.  Each group's log is kept.
set -uo pipefail
source /home/cz/miniconda3/etc/profile.d/conda.sh

LOGDIR=/home/cz/WGS/logs/install
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/SUMMARY.txt"
: > "$SUMMARY"

inst() {
  local name="$1"; shift
  echo "### GROUP $name : $* ###"
  if mamba install -y -n WGS -c conda-forge -c bioconda "$@" > "$LOGDIR/$name.log" 2>&1; then
    echo "OK      $name" | tee -a "$SUMMARY"
  else
    echo "FAILED  $name" | tee -a "$SUMMARY"
    grep -iE 'nothing provides|PackagesNotFound|conflict|error|failed' "$LOGDIR/$name.log" \
      | head -5 | sed 's/^/         /'
  fi
}

echo "=== install groups started $(date) ==="

inst 01_core      bwa minimap2 bedtools seqtk pigz seqkit
inst 02_var       freebayes snp-sites mummer4
inst 03_quast     quast
inst 04_python    python=3.11 numpy pandas scipy matplotlib pysam biopython
inst 05_prokka    prokka
inst 06_bakta     bakta
inst 07_snpeff    snpeff
inst 08_utils     entrez-direct parallel wgsim
inst 09_amr       abricate ncbi-amrfinderplus blast
inst 10_r         r-base=4.3 r-ggplot2 r-dplyr r-tidyr r-circlize r-gggenes r-ggrepel r-scales r-patchwork r-reshape2 r-rcolorbrewer r-readr r-stringr r-purrr r-tibble
inst 11_phylo     mafft trimal iqtree mash fastANI
inst 12_pangen    roary cd-hit diamond
inst 13_extra     minced phispy
inst 14_mob       mob_suite

echo "=== install groups finished $(date) ==="
cat "$SUMMARY"
