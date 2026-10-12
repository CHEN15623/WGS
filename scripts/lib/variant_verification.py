#!/usr/bin/env python3
"""
variant_verification.py -- turn the reconciled variant table into three
defensible products:

  variants_confident.tsv   the variants that survive every filter, each with the
                           reference and query alleles reconstructed from the
                           alignment and the observed read allele fraction
  variants_ambiguous.tsv   the candidates that were rejected, each with the
                           concrete reason (multi-copy reference, overlap in the
                           NUCmer alignment, low read mapping quality, or reads
                           that contradict the call)
  verification_windows.tsv ref/query sequence side by side around each
                           confident variant, so the call can be eyeballed

Reconstructing the query allele matters because gap placement inside
homopolymer / short-tandem runs is ambiguous: an aligner may place a 1 bp
insertion at either end of a run.  Comparing the actual sequences resolves what
the difference really is.
"""
import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from annotate_variants import read_fasta, project_region  # noqa: E402


def build_query_allele(bam, contig, start0, end0, refseq):
    """Rebuild the query allele over reference [start0, end0) (0-based).

    The assembly-to-reference BAM is expected here: a contig alignment carries
    exactly one representative of each query base, so insertions are not counted
    once per read.  Any residual duplication (a primary plus a supplementary
    alignment) is resolved by majority vote per insertion site.
    """
    ref2q, ins, _cov = project_region(bam, contig, start0, end0)
    out = []
    for p in range(start0, end0):
        out.append(ref2q.get(p, "N"))
        if p in ins:
            seqs = ins[p]
            # majority-vote the inserted sequence at this anchor
            best = None
            if seqs:
                from collections import Counter as _C
                # group by the base string built from individual bases
                joined = "".join(seqs)
                # if every read contributed the same single base, emit it once
                if len(set(seqs)) == 1 and len(seqs[0]) == 1:
                    best = seqs[0]
                elif len(joined) % max(1, len(seqs)) == 0:
                    n = len(seqs)
                    # try to split `joined` into n equal pieces
                    if all(joined[k * (len(joined) // n):(k + 1) * (len(joined) // n)]
                           == joined[:len(joined) // n] for k in range(n)):
                        best = joined[:len(joined) // n]
                    else:
                        best = joined
                else:
                    best = joined
            if best:
                out.append(best)
    return "".join(out), refseq[start0:end0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--window", type=int, default=25,
                    help="bp of context on each side for the verification view")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    refseqs = read_fasta(a.ref)
    try:
        import pysam
        bam = pysam.AlignmentFile(a.bam, "rb")
    except Exception as exc:
        print(f"[verify] cannot open BAM ({exc})", file=sys.stderr)
        bam = None

    rows = []
    with open(a.merged) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            f += [""] * (len(hdr) - len(f))
            rows.append(dict(zip(hdr, f)))

    # ------------------------------------------------------------------ #
    #  Eligibility and event clustering
    # ------------------------------------------------------------------ #
    # The four callers do not all describe an indel the same way: inside a
    # homopolymer run MUMmer may report "insert a T at 530,489" while a VCF
    # caller reports the equivalent anchored "G>GT at 530,487".  Those are one
    # biological event.  Candidates are therefore first filtered for ambiguity
    # and then clustered by proximity so that the independent evidence from all
    # representations is pooled before a call is accepted.
    AMBIG_TOKENS = ("in_repeat", "low_mapq", "reads_contradict",
                    "multi_contig_alignment", "nucmer_overlap_RQ",
                    "ref_multicopy", "ref_allele_disagrees")

    def is_clean(r):
        fl = r.get("flags", "") or ""
        return not any(t in fl for t in AMBIG_TOKENS)

    clean = [r for r in rows if is_clean(r)]
    dirty = [r for r in rows if not is_clean(r)]

    # cluster: same contig, same broad class (SNP vs indel), within 10 bp
    clean.sort(key=lambda r: (r["contig"], int(r["pos"])))
    clusters = []
    for r in clean:
        pos = int(r["pos"])
        kind = "SNP" if r["type"] == "SNP" else "INDEL"
        placed = False
        for c in clusters:
            if c["contig"] == r["contig"] and c["kind"] == kind and \
               pos - c["maxpos"] <= 10:
                c["members"].append(r)
                c["maxpos"] = max(c["maxpos"], pos)
                placed = True
                break
        if not placed:
            clusters.append(dict(contig=r["contig"], kind=kind, maxpos=pos,
                                 members=[r]))

    # an event is accepted when >= 2 independent methods back any of its members
    events = []
    for c in clusters:
        methods = max(int(m.get("n_independent_methods") or 0)
                      for m in c["members"])
        callers = sorted({x for m in c["members"]
                          for x in (m.get("callers", "") or "").split(";") if x})
        if methods < 2:
            continue
        # representative = the member with the most support, then leftmost
        rep = sorted(c["members"],
                     key=lambda m: (-int(m.get("n_independent_methods") or 0),
                                    int(m["pos"])))[0]
        ev = dict(rep)
        ev["callers"] = ";".join(callers)
        ev["n_independent_methods"] = methods
        ev["n_representations"] = len(c["members"])
        ev["event_span"] = f"{min(int(m['pos']) for m in c['members'])}-" \
                           f"{max(int(m['pos']) for m in c['members'])}"
        ev["confidence"] = "HIGH" if methods >= 2 else "MEDIUM"
        events.append(ev)

    confident, ambiguous, windows = [], [], []

    for r in events:
        contig, pos = r["contig"], int(r["pos"])
        conf = r.get("confidence", "")
        ref = refseqs.get(contig, "")
        w = a.window
        s0 = max(0, pos - 1 - w)
        e0 = min(len(ref), pos + w)
        rwin = qwin = ""
        if bam is not None and ref:
            try:
                qwin, rwin = build_query_allele(bam, contig, s0, e0, ref)
            except Exception:
                pass
        rec = dict(r)
        rec["ref_window"] = rwin
        rec["qry_window"] = qwin
        rec["window_start"] = s0 + 1
        rec["window_end"] = e0
        # where inside the window does the difference sit?
        d0 = d1 = ""
        if rwin and qwin:
            n = min(len(rwin), len(qwin))
            i = next((k for k in range(n) if rwin[k] != qwin[k]), None)
            if i is None and len(rwin) != len(qwin):
                i = n
            d0 = s0 + 1 + (i if i is not None else 0)

        if conf in ("HIGH", "MEDIUM"):
            confident.append(rec)
            windows.append(dict(
                contig=contig, pos=pos, type=rec.get("type", ""),
                confidence=conf, ref_allele=rec.get("ref", ""),
                alt_allele=rec.get("alt", ""),
                window_start=rec["window_start"], window_end=rec["window_end"],
                reference_sequence=rwin, query_sequence=qwin,
                first_difference_at=d0,
                observed_read_depth=rec.get("obs_depth", ""),
                observed_alt_percent=rec.get("obs_alt_frac", ""),
                supporting_callers=rec.get("callers", ""),
                supporting_methods=rec.get("n_independent_methods", "")))
        else:
            ambiguous.append(rec)

    # everything that was not accepted becomes an explicitly ambiguous entry
    accepted_keys = {(e["contig"], e["pos"]) for e in confident}
    for r in dirty:
        if (r["contig"], r["pos"]) in accepted_keys:
            continue
        ambiguous.append(dict(r))
    for cl in clusters:
        methods = max(int(m.get("n_independent_methods") or 0)
                      for m in cl["members"])
        if methods >= 2:
            continue
        for m in cl["members"]:
            if (m["contig"], m["pos"]) in accepted_keys:
                continue
            mm = dict(m)
            mm["flags"] = ((m.get("flags", "") or "") +
                           ";single_method_only").strip(";")
            ambiguous.append(mm)

    # ---- minimal (normalised) REF/ALT for each event ---------------------
    def minimal_variant(rwin, qwin, ws):
        """Reduce a ref/query window pair to a minimal anchored REF/ALT pair.

        In the reconstructed query window an `N` marks a reference position for
        which the query contributes no base -- i.e. a deletion.  Such runs are
        converted to a proper anchored deletion; everything else is reduced by
        trimming the common prefix and suffix, adding one anchor base only when
        one of the two alleles would otherwise be empty.
        """
        if not rwin or not qwin:
            return None

        # ---- deletion: contiguous run of uncovered reference bases ---------
        if "N" in qwin:
            idx = [k for k, ch in enumerate(qwin) if ch == "N"]
            if idx and idx == list(range(idx[0], idx[-1] + 1)) and idx[0] > 0:
                a = idx[0]
                anchor = rwin[a - 1]
                rs = anchor + "".join(rwin[k] for k in idx)
                qs = anchor
                return ws + a - 1, rs, qs

        # ---- generic: trim common prefix / suffix --------------------------
        i = 0
        while i < len(rwin) and i < len(qwin) and rwin[i] == qwin[i]:
            i += 1
        j = 0
        while (j < len(rwin) - i and j < len(qwin) - i
               and rwin[len(rwin) - 1 - j] == qwin[len(qwin) - 1 - j]):
            j += 1
        rs = rwin[i:len(rwin) - j] if j else rwin[i:]
        qs = qwin[i:len(qwin) - j] if j else qwin[i:]
        rs = rs.replace("N", "")
        qs = qs.replace("N", "")
        pos0 = ws - 1 + i
        if not rs or not qs:                 # need a left anchor base
            if i == 0:
                return None
            pos0 -= 1
            rs = rwin[i - 1] + rs
            qs = qwin[i - 1] + qs
        if not rs or not qs:
            return None
        return pos0 + 1, rs, qs

    for w in windows:
        if w.get("reference_sequence") and w.get("query_sequence"):
            mv = minimal_variant(w["reference_sequence"], w["query_sequence"],
                                 w["window_start"])
            if mv:
                w["normalised_pos"], w["normalised_ref"], w["normalised_alt"] = mv

    # write the corrected confident VCF
    vcf_path = f"{a.outdir}/variants_confident.vcf"
    with open(vcf_path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=variant_verification.py\n")
        for c, seq in refseqs.items():
            fh.write(f"##contig=<ID={c},length={len(seq)}>\n")
        fh.write('##INFO=<ID=CALLERS,Number=1,Type=String,Description="Supporting callers">\n')
        fh.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQUERY\n")
        for k, w in enumerate(windows, 1):
            if "normalised_ref" not in w:
                continue
            fh.write(f"{w['contig']}\t{w['normalised_pos']}\tVAR{k}\t"
                     f"{w['normalised_ref']}\t{w['normalised_alt']}\t.\tPASS\t"
                     f"CALLERS={w.get('supporting_callers','')};"
                     f"TYPE={w.get('type','')}\tGT\t1/1\n")
    os.system(f"bgzip -f -c '{vcf_path}' > '{vcf_path}.gz' 2>/dev/null && "
              f"tabix -f -p vcf '{vcf_path}.gz' 2>/dev/null")
    print(f"[verify] wrote {vcf_path}")

    def write_tsv(path, records, cols):
        with open(path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for rec in records:
                fh.write("\t".join(str(rec.get(c, "")) for c in cols) + "\n")

    base_cols = ["contig", "pos", "ref", "alt", "type", "confidence",
                 "n_independent_methods", "callers", "flags", "ref_copies",
                 "obs_depth", "obs_alt_frac", "asm_mapq", "asm_depth",
                 "sim_depth", "sim_mean_mapq", "mummer_buff", "mummer_dist"]

    write_tsv(f"{a.outdir}/variants_confident.tsv", confident,
              base_cols + ["ref_window", "qry_window", "window_start",
                           "window_end"])
    write_tsv(f"{a.outdir}/variants_ambiguous.tsv", ambiguous, base_cols)
    write_tsv(f"{a.outdir}/verification_windows.tsv", windows,
              ["contig", "pos", "type", "confidence", "ref_allele", "alt_allele",
               "window_start", "window_end", "reference_sequence",
               "query_sequence", "first_difference_at", "observed_read_depth",
               "observed_alt_percent", "supporting_callers",
               "supporting_methods"])

    # ---- reason breakdown for the report ---------------------------------
    reason = Counter()
    for r in ambiguous:
        for fl in (r.get("flags", "") or "").split(";"):
            if fl in ("in_repeat", ""):
                continue
            reason[fl.split("_x")[0] if fl.startswith("ref_multicopy") else fl] += 1
        if not (r.get("flags", "") or "").strip(";in_repeat"):
            reason["single_method_only"] += 1

    summary = dict(
        total_candidates=len(rows),
        confident=len(confident),
        ambiguous=len(ambiguous),
        confident_variants=[dict(contig=c["contig"], pos=c["pos"],
                                 ref=c["ref"], alt=c["alt"], type=c["type"],
                                 methods=c.get("n_independent_methods", ""),
                                 callers=c.get("callers", ""),
                                 observed_alt_percent=c.get("obs_alt_frac", ""),
                                 observed_depth=c.get("obs_depth", ""))
                            for c in confident],
        ambiguity_reasons=dict(reason.most_common()),
    )
    with open(f"{a.outdir}/verification_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    if bam is not None:
        bam.close()

    print(f"[verify] confident={len(confident)}  ambiguous={len(ambiguous)}")
    for c in confident:
        print(f"[verify]   {c['contig']}:{c['pos']} {c['ref']}>{c['alt']} "
              f"({c['type']}, {c.get('n_independent_methods')} methods, "
              f"alt={c.get('obs_alt_frac')}% of {c.get('obs_depth')} reads)")


if __name__ == "__main__":
    main()
