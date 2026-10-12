#!/usr/bin/env bash
# =============================================================================
#  01_qc.sh -- QUALITY CONTROL
# =============================================================================
#  Two layers of QC:
#    (a) ASSEMBLY QC   -- the primary, real data.  Assembly statistics
#        (N50/L50/GC/ambiguous bases), QUAST single-assembly and
#        reference-based comparative reports.
#    (b) READ QC       -- QC of in-silico Illumina reads generated from the
#        query assembly.  These reads exist to (i) exercise a genuine
#        read-level QC step and (ii) give an independent, depth-supported
#        confirmation of the assembly-derived variant calls.  Everything
#        derived from them is tagged `SIM`.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

QCD="$OUT/01_qc"
mkdir -p "$QCD/assembly_qc" "$QCD/quast" "$QCD/reads_qc" \
         "$QCD/annotation_qc" "$WORK/reads"

# -----------------------------------------------------------------------------
msg "01-QC (a): assembly statistics"
# -----------------------------------------------------------------------------
python3 "$SCRIPTS/lib/assembly_stats.py" "$REF_FA" "reference_${REF_NAME}" "$QCD/assembly_qc"
python3 "$SCRIPTS/lib/assembly_stats.py" "$QRY_FA" "query_${QRY_NAME}"     "$QCD/assembly_qc"

# combined side-by-side table
python3 - "$QCD/assembly_qc/reference_${REF_NAME}_assembly_stats.tsv" \
           "$QCD/assembly_qc/query_${QRY_NAME}_assembly_stats.tsv" \
           "$QCD/assembly_qc/assembly_stats_comparison.tsv" <<'PY'
import sys
ref_f, qry_f, out = sys.argv[1], sys.argv[2], sys.argv[3]
def load(p):
    d = {}
    with open(p) as fh:
        next(fh)
        for line in fh:
            k, v = line.rstrip("\n").split("\t")
            d[k] = v
    return d
r, q = load(ref_f), load(qry_f)
with open(out, "w") as fh:
    fh.write("metric\treference\tquery\tdifference\n")
    for k in r:
        rv, qv = r[k], q.get(k, "")
        try:
            diff = f"{float(qv)-float(rv):+.4g}"
        except ValueError:
            diff = ""
        fh.write(f"{k}\t{rv}\t{qv}\t{diff}\n")
print(f"[qc] wrote {out}")
PY

# -----------------------------------------------------------------------------
msg "01-QC (b): QUAST assembly reports"
# -----------------------------------------------------------------------------
if command -v quast.py >/dev/null 2>&1; then
  # QUAST 5.x already uses Minimap2 internally; there is no --minimap2 flag.
  QUAST_ALIGNER=""

  # single-assembly reports
  run_logged quast_reference \
    quast.py -o "$QCD/quast/reference" -t "$THREADS" $QUAST_ALIGNER \
      --no-icarus --silent "$REF_FA" || true
  run_logged quast_query \
    quast.py -o "$QCD/quast/query" -t "$THREADS" $QUAST_ALIGNER \
      --no-icarus --silent "$QRY_FA" || true

  # reference-based comparative report: query vs reference
  run_logged quast_comparative \
    quast.py -o "$QCD/quast/comparative" -t "$THREADS" $QUAST_ALIGNER \
      -r "$REF_FA" --no-icarus --silent "$QRY_FA" || true

  # reference genes supplied so QUAST reports completeness of the annotation
  run_logged quast_comparative_genes \
    quast.py -o "$QCD/quast/comparative_genes" -t "$THREADS" $QUAST_ALIGNER \
      -r "$REF_FA" -g "$REF_GFF" --no-icarus --silent "$QRY_FA" || true

  for d in reference query comparative comparative_genes; do
    if [ -f "$QCD/quast/$d/report.tsv" ]; then
      cp "$QCD/quast/$d/report.tsv" "$QCD/quast/${d}_report.tsv"
    fi
    if [ -f "$QCD/quast/$d/report.txt" ]; then
      cp "$QCD/quast/$d/report.txt" "$QCD/quast/${d}_report.txt"
    fi
  done
else
  warn "quast.py not found -- skipping QUAST (assembly_stats.py output still produced)"
fi

# -----------------------------------------------------------------------------
msg "01-QC (c): in-silico reads from the query assembly + read-level QC"
# -----------------------------------------------------------------------------
READ1="$WORK/reads/${QRY_NAME}_SIM_R1.fastq"
READ2="$WORK/reads/${QRY_NAME}_SIM_R2.fastq"
STATS="$QCD/reads_qc/sim_reads_stats.json"

if [ ! -s "$READ1" ] && command -v wgsim >/dev/null 2>&1; then
  GSIZE=$(awk '/^>/{next}{n+=length($0)}END{print n}' "$QRY_FA")
  NREADS=$(python3 -c "print(int($SIM_COVERAGE*$GSIZE/(2*$SIM_READLEN)))")
  msg "  simulating $NREADS read pairs (~${SIM_COVERAGE}x, ${SIM_READLEN}bp, insert ${SIM_INSERTSIZE}bp)"
  # -e 0.005 : realistic 0.5 % base error rate (this also drives wgsim's coarse
  #            quality encoding, ~Phred 23 at this rate)
  # -r/-R/-X 0 : no extra mutations or indels -- the reads must faithfully
  #            represent the query assembly, otherwise the read-based variant
  #            calls would not be an independent check of the assembly ones
  run_logged wgsim \
    wgsim -N "$NREADS" -1 "$SIM_READLEN" -2 "$SIM_READLEN" \
          -d "$SIM_INSERTSIZE" -e 0.005 -S 20240512 -r 0 -R 0 -X 0 \
          "$QRY_FA" "$READ1" "$READ2" || warn "wgsim failed"
fi

if [ -s "$READ1" ]; then
  command -v seqkit >/dev/null 2>&1 && \
    run_logged seqkit_stats_reads seqkit stats -a -T "$READ1" "$READ2" || true

  if command -v fastqc >/dev/null 2>&1; then
    run_logged fastqc_reads \
      fastqc --threads "$THREADS" --outdir "$QCD/reads_qc" --quiet \
             "$READ1" "$READ2" || true
  else
    warn "fastqc not installed -- raw-read QC report skipped"
  fi

  if command -v fastp >/dev/null 2>&1; then
    run_logged fastp_reads \
      fastp -i "$READ1" -I "$READ2" \
            -o "$WORK/reads/${QRY_NAME}_SIM_R1.clean.fq.gz" \
            -O "$WORK/reads/${QRY_NAME}_SIM_R2.clean.fq.gz" \
            -j "$QCD/reads_qc/fastp.json" -h "$QCD/reads_qc/fastp.html" \
            -w "$THREADS" -q 20 -u 30 -l 50 || true
  fi
fi

# -----------------------------------------------------------------------------
msg "01-QC (d): reference annotation QC"
# -----------------------------------------------------------------------------
python3 - "$REF_GFF" "$QCD/annotation_qc/reference_annotation_summary.tsv" <<'PY'
import sys, collections
gff, out = sys.argv[1], sys.argv[2]
counts = collections.Counter()
products = 0
with open(gff) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        counts[f[2]] += 1
        if f[2] == "CDS" and "product=" in f[8]:
            products += 1
with open(out, "w") as fh:
    fh.write("feature_type\tcount\n")
    for k, v in counts.most_common():
        fh.write(f"{k}\t{v}\n")
    fh.write(f"CDS_with_product_annotation\t{products}\n")
print(f"[qc] annotation features: {dict(counts.most_common(6))}")
PY

msg "01-QC complete -> $QCD"
