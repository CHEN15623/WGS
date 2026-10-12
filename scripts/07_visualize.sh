#!/usr/bin/env bash
# =============================================================================
#  07_visualize.sh -- every figure of the report (pure matplotlib)
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

FIG="$OUT/07_visualization"
mkdir -p "$FIG"

msg "07-VIZ: generating figures"

python3 "$SCRIPTS/lib/make_figures.py" \
  --out      "$OUT" \
  --figdir   "$FIG" \
  --ref-fa   "$REF_FA" \
  --ref-gff  "$REF_GFF" \
  --ref-name "$REF_NAME" \
  --qry-name "$QRY_NAME" \
  --ref-chr  "$REF_CHR" || warn "figure generation reported errors"

n=$(ls -1 "$FIG"/*.png 2>/dev/null | wc -l)
msg "07-Visualisation complete: $n figure(s) -> $FIG"

# an index so the figures can be browsed
{
  echo "# Figure index"
  echo
  ls -1 "$FIG"/*.png 2>/dev/null | while read -r f; do
    b=$(basename "$f" .png)
    echo "- [$b]($b.png)"
  done
} > "$FIG/INDEX.md"
