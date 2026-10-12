#!/usr/bin/env bash
# =============================================================================
#  08_report.sh -- assemble the final Markdown + HTML report
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

REP="$OUT/08_report"
mkdir -p "$REP"

msg "08-REPORT: assembling the final report"

python3 "$SCRIPTS/lib/make_report.py" \
  --out       "$OUT" \
  --ref-name  "$REF_NAME" \
  --qry-name  "$QRY_NAME" \
  --ref-acc   "$REF_ACC" \
  --qry-acc   "$QRY_ACC" \
  --ref-chr   "$REF_CHR" \
  --ref-plas  "$REF_PLASMID" \
  --env-name  "${CONDA_DEFAULT_ENV:-WGS}" \
  --sim-cov   "$SIM_COVERAGE" \
  --sim-len   "$SIM_READLEN" \
  --sim-ins   "$SIM_INSERTSIZE" || die "report generation failed"

# a compact digest of every headline number, for quick inspection
{
  echo "WGS run summary -- $(date '+%Y-%m-%d %H:%M:%S')"
  echo "================================================================"
  echo
  echo "-- variant reconciliation (03_variants/merged/merge_summary.json)"
  python3 - <<'PY' 2>/dev/null || true
import json, os
p = os.environ.get("OUT", "") + "/03_variants/merged/merge_summary.json"
if os.path.exists(p):
    d = json.load(open(p))
    for k, v in d.items():
        if k != "per_caller":
            print(f"   {k}: {v}")
    print("   per_caller:", d.get("per_caller"))
PY
  echo
  echo "-- verified variants (03_variants/verified/variants_confident.tsv)"
  if [ -s "$OUT/03_variants/verified/variants_confident.tsv" ]; then
    awk -F'\t' 'NR>1{printf "   %s:%s %s>%s (%s) methods=%s\n",$1,$2,$3,$4,$5,$7}' \
      "$OUT/03_variants/verified/variants_confident.tsv"
  fi
  echo
  echo "-- affected genes (04_annotation/variant_effect/gene_impact.tsv)"
  if [ -s "$OUT/04_annotation/variant_effect/gene_impact.tsv" ]; then
    awk -F'\t' 'NR>1{printf "   %s  %s  %s / %s  %s\n",$1,$3,$11,$12,$21}' \
      "$OUT/04_annotation/variant_effect/gene_impact.tsv"
  fi
  echo
  echo "-- structural differences (02_alignment/synteny/structural_variants.tsv)"
  if [ -s "$OUT/02_alignment/synteny/structural_variants.tsv" ]; then
    awk -F'\t' 'NR>1{printf "   %-36s %s:%s-%s  %s bp\n",$1,$2,$3,$4,$8}' \
      "$OUT/02_alignment/synteny/structural_variants.tsv"
  fi
  echo
  echo "-- figures"
  ls -1 "$OUT/07_visualization"/*.png 2>/dev/null | sed 's|.*/|   |'
} > "$REP/RUN_SUMMARY.txt" 2>/dev/null || true

msg "08-Report complete -> $REP/WGS_report.html"
ls -la "$REP"/WGS_report.* 2>/dev/null || true
