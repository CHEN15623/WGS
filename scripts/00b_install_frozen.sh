#!/usr/bin/env bash
# Install the remaining WGS toolchain using --freeze-installed, which keeps the
# already-solved environment fixed and therefore resolves in seconds instead of
# minutes.  Each group is attempted on its own so one conflict cannot block the
# rest; anything that still fails is reported at the end.
set -uo pipefail
source /home/cz/miniconda3/etc/profile.d/conda.sh

LOGDIR=/home/cz/WGS/logs/install
mkdir -p "$LOGDIR"
SUMMARY="$LOGDIR/SUMMARY2.txt"
: > "$SUMMARY"

inst() {
  local name="$1"; shift
  echo "### $name : $*"
  if mamba install -y -n WGS -c conda-forge -c bioconda --freeze-installed "$@" \
        > "$LOGDIR/f_$name.log" 2>&1; then
    echo "OK      $name" | tee -a "$SUMMARY"
  else
    echo "FAILED  $name" | tee -a "$SUMMARY"
    grep -iE 'nothing provides|PackagesNotFound|conflict|error' "$LOGDIR/f_$name.log" \
      | head -3 | sed 's/^/         /'
  fi
}

echo "=== frozen installs started $(date) ==="

inst prokka        prokka
inst bakta         bakta
inst snpeff        snpeff
inst amr           abricate ncbi-amrfinderplus
inst utils         entrez-direct parallel wgsim
inst phylo         mafft trimal iqtree mash fastani
inst pangen        roary cd-hit diamond
inst phispy        phispy
inst minced        minced
inst r             r-base=4.3 r-ggplot2 r-dplyr r-tidyr r-circlize r-gggenes r-ggrepel r-scales r-patchwork r-reshape2 r-rcolorbrewer r-readr r-stringr r-purrr r-tibble
inst r_extra       r-ggnewscale r-viridis r-gridextra r-cowplot r-ggsci
inst py_extra      python-kaleido plotly networkx scikit-learn statsmodels

echo "=== frozen installs finished $(date) ==="
cat "$SUMMARY"
echo "--- tool check ---"
conda run -n WGS bash -c '
for t in prokka bakta snpEff abricate amrfinder mafft iqtree mash fastANI \
         roary diamond cd-hit phispy minced Rscript quast.py fastqc fastp; do
  printf "%-12s %s\n" "$t" "$(command -v $t || echo MISSING)"
done'
