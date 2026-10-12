#!/usr/bin/env bash
# =============================================================================
#  05_functional.sh -- FUNCTIONAL PREDICTION AND ANNOTATION
#                       (antimicrobial resistance / virulence / metabolism)
# =============================================================================
#  Both genomes are screened with the same batteries so that every difference is
#  attributable to the strain and not to the tool:
#    * ABRicate         : CARD, ResFinder, NCBI-AMRFinder, ARG-ANNOT, VFDB,
#                         MEGARES, PlasmidFinder, ISfinder, Victors
#    * AMRFinderPlus    : NCBI curated AMR + virulence + stress + biocide
#    * curated scan     : B. fragilis-specific determinants (cfiA, nim, ermF,
#                         tetQ, cfxA, bft, gyrA/parC/rpoB QRDR, susC/susD, CPS)
#    * Bakta/Prokka     : COG / KEGG / GO / EC -> metabolic and functional
#                         category comparison between reference and query
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

FUN="$OUT/05_functional"
mkdir -p "$FUN/abricate" "$FUN/amrfinder" "$FUN/curated_scan" \
         "$FUN/cog_kegg" "$FUN/comparison"

ABRICATE_DBS="card resfinder ncbi argannot vfdb megares plasmidfinder isfinder victors"

# -----------------------------------------------------------------------------
msg "05-FUN (1): ABRicate screens on both genomes"
# -----------------------------------------------------------------------------
if command -v abricate >/dev/null 2>&1; then
  # make sure the databases are present locally
  if ! abricate --list 2>/dev/null | grep -q '^card'; then
    msg "  setting up ABRicate databases ..."
    run_logged abricate_setupdb abricate --setupdb || warn "abricate --setupdb failed"
  fi

  for genome in "reference:$REF_FA" "query:$QRY_FA"; do
    label="${genome%%:*}"; fasta="${genome#*:}"
    for db in $ABRICATE_DBS; do
      out="$FUN/abricate/${label}.${db}.tab"
      [ -s "$out" ] && continue
      # a database may legitimately be unavailable; record that instead of dying
      if abricate --db "$db" --minid 80 --mincov 80 "$fasta" > "$out" 2> "$LOGDIR/abricate_${label}_${db}.err"; then
        echo "   [ok]   abricate $label/$db ($(($(wc -l < "$out")-1)) hits)"
      else
        echo "" > "$out"; warn "abricate $label/$db failed"
      fi
    done
  done
else
  warn "abricate not installed"
fi

# -----------------------------------------------------------------------------
msg "05-FUN (2): NCBI AMRFinderPlus on both genomes"
# -----------------------------------------------------------------------------
if command -v amrfinder >/dev/null 2>&1; then
  # Locate the database that amrfinder_update installed.  `--organism` is NOT
  # used: Bacteroides fragilis has no curated point-mutation set in
  # AMRFinderPlus, and passing an unsupported name aborts the whole run.
  AMRDB=""
  for cand in "$CONDA_PREFIX"/bin/data/*/ "$CONDA_PREFIX"/share/amrfinderplus/data/*/; do
    [ -s "${cand}AMR.LIB" ] && AMRDB="$cand" && break
  done
  [ -n "$AMRDB" ] || AMRDB="$CONDA_PREFIX/bin/data/latest"
  if [ ! -s "${AMRDB}AMR.LIB" ]; then
    msg "  updating the AMRFinderPlus database ..."
    run_logged amrfinder_update amrfinder_update --force_update || \
      warn "amrfinder_update failed"
    for cand in "$CONDA_PREFIX"/bin/data/*/; do
      [ -s "${cand}AMR.LIB" ] && AMRDB="$cand" && break
    done
  fi
  echo "  AMRFinderPlus database: $AMRDB"
  for genome in "reference:$REF_FA" "query:$QRY_FA"; do
    label="${genome%%:*}"; fasta="${genome#*:}"
    out="$FUN/amrfinder/${label}.amrfinder.tsv"
    [ -s "$out" ] && continue
    run_logged "amrfinder_$label" bash -c \
      "amrfinder --database '$AMRDB' --nucleotide '$fasta' --plus \
         --output '$out' --threads $THREADS" || \
    run_logged "amrfinder_${label}_basic" bash -c \
      "amrfinder --database '$AMRDB' --nucleotide '$fasta' \
         --output '$out' --threads $THREADS" || warn "amrfinder $label failed"
  done
else
  warn "amrfinder not installed"
fi

# -----------------------------------------------------------------------------
msg "05-FUN (3): curated B. fragilis determinant scan (both proteomes)"
# -----------------------------------------------------------------------------
# whichever annotation succeeded (Bakta preferred, Prokka as fallback)
QRY_ANN_SRC="$OUT/04_annotation/query_annotation_source.txt"
QRY_GFF_ANN=$(cat "$QRY_ANN_SRC" 2>/dev/null || echo "")
QRY_FAA="$OUT/04_annotation/bakta_query/${QRY_NAME}.faa"
[ -s "$QRY_FAA" ] || QRY_FAA="$OUT/04_annotation/prokka_query/${QRY_NAME}.faa"
[ -s "$QRY_FAA" ] || QRY_FAA="$OUT/04_annotation/transferred/query.faa"
[ -s "$QRY_FAA" ] || QRY_FAA="$REF_FAA"
[ -n "$QRY_GFF_ANN" ] || QRY_GFF_ANN="$REF_GFF"
echo "  using query annotation: $QRY_GFF_ANN"

python3 "$SCRIPTS/lib/curated_determinants.py" \
  --ref-faa "$REF_FAA" \
  --ref-gff "$REF_GFF" \
  --qry-faa "$QRY_FAA" \
  --qry-gff "$QRY_GFF_ANN" \
  --qry-fna "$QRY_FA" \
  --ref-fna "$REF_FA" \
  --outdir "$FUN/curated_scan" \
  --threads "$THREADS" || warn "curated determinant scan failed"

# -----------------------------------------------------------------------------
msg "05-FUN (4): metabolic / functional category comparison (Bakta COG-KEGG-GO)"
# -----------------------------------------------------------------------------
python3 "$SCRIPTS/lib/functional_compare.py" \
  --ref-gff "$REF_GFF" \
  --qry-bakta "$OUT/04_annotation/bakta_query/${QRY_NAME}.gff3" \
  --qry-bakta-tsv "$OUT/04_annotation/bakta_query/${QRY_NAME}.tsv" \
  --ref-bakta-tsv "$OUT/04_annotation/bakta_reference/${REF_NAME}.tsv" \
  --qry-prokka "$QRY_GFF_ANN" \
  --abricate-dir "$FUN/abricate" \
  --amrfinder-dir "$FUN/amrfinder" \
  --outdir "$FUN/comparison" || warn "functional comparison failed"

# -----------------------------------------------------------------------------
msg "05-FUN (5): COG functional-category profiling (NCBI COG2020, via DIAMOND)"
# -----------------------------------------------------------------------------
msg "  reference proteome : $REF_FAA"
msg "  query proteome     : $QRY_FAA"
python3 "$SCRIPTS/lib/cog_profile.py" \
  --cog-dir "$WORK/cog" \
  --ref-faa "$REF_FAA" \
  --qry-faa "$QRY_FAA" \
  --outdir "$FUN/cog_kegg" --threads "$MAXTHREADS" || \
  warn "COG profiling skipped (database not available)"

# -----------------------------------------------------------------------------
msg "05-FUN (6): GO / EC profile comparison"
# -----------------------------------------------------------------------------
python3 "$SCRIPTS/lib/go_ec_profile.py" \
  --ref-gff "$REF_GFF" \
  --qry-gff "$QRY_GFF_ANN" \
  --outdir "$FUN/go_ec" || warn "GO/EC profiling failed"

msg "05-Functional analysis complete -> $FUN"
