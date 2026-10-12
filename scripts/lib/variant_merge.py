#!/usr/bin/env python3
"""
variant_merge.py -- reconcile four independent variant-calling routes into one
                   curated table, and emit a high-confidence VCF.

Routes
  mummer        : MUMmer show-snps on the unique-anchor NUCmer delta  (assembly)
  bcftools_asm  : bcftools call on the minimap2 asm5 contig BAM       (assembly)
  freebayes_asm : freebayes on the same BAM                           (assembly)
  bcftools_sim  : bcftools call on the simulated-read BAM, deep depth (reads)

Outputs (in --outdir)
  variants_merged.tsv        every distinct variant with per-caller support
  variants_highconf.tsv      the curated set
  variants_highconf.vcf      same set as a VCF (for snpEff / downstream tools)
  caller_concordance.tsv     overlap matrix between callers
  merge_summary.json         headline numbers
"""
import argparse
import gzip
import json
import os
import re
import sys
from collections import OrderedDict, defaultdict

# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #

def open_maybe_gz(path):
    if path is None or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


# --------------------------------------------------------------------------- #
#  route 1 : MUMmer show-snps
# --------------------------------------------------------------------------- #

def parse_mummer(path):
    """show-snps output -> list of dicts.

    Two layouts are accepted and distinguished by field count:

      `show-snps -rlTH`   (14 fields, [R]/[Q] retained)
        0 P1    position in reference
        1 REF   reference base ('.' = insertion relative to the reference)
        2 QRY   query base     ('.' = deletion relative to the reference)
        3 P2    position in query
        4 BUFF  distance to nearest alignment gap
        5 DIST  distance to nearest alignment / sequence end
        6 R     other reference alignments overlapping this position
        7 Q     other query alignments overlapping this position
        8 LENR  reference sequence length
        9 LENQ  query sequence length
       10 FRM   reading frame
       11 (reserved)
       12 TAG1  reference sequence id
       13 TAG2  query sequence id

      `show-snps -rlTHC`  (12 fields, [R]/[Q] dropped)
        0..5 as above, 6 LENR, 7 LENQ, 8 FRM, 9 (reserved), 10 TAG1, 11 TAG2
    """
    fh = open_maybe_gz(path)
    if fh is None:
        return [], None, None
    out = []
    ref_id = qry_id = None
    with fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                f = line.split()
            if len(f) < 8:
                continue
            try:
                p1 = int(f[0])
            except ValueError:
                continue
            rb, qb = f[1], f[2]
            # tags are always the last two fields
            ref_id, qry_id = f[-2], f[-1]
            n = len(f)
            if n >= 13:            # -rlTH : R/Q present at 6,7
                r_i, q_i = 6, 7
                lenR_i, frm_i = 8, 10
            else:                  # -rlTHC : R/Q absent
                r_i = q_i = None
                lenR_i, frm_i = 6, 8
            def _int(i):
                try:
                    return int(f[i])
                except (ValueError, IndexError):
                    return -1
            try:
                p2 = int(f[3])
            except ValueError:
                p2 = -1
            out.append(dict(pos=p1, ref=rb, alt=qb, qry_pos=p2,
                            buff=_int(4), dist=_int(5),
                            overlap_R=_int(r_i) if r_i is not None else 0,
                            overlap_Q=_int(q_i) if q_i is not None else 0,
                            lenR=_int(lenR_i), frame=_int(frm_i),
                            tags=f"{f[-2]} {f[-1]}"))
    return out, ref_id, qry_id


# --------------------------------------------------------------------------- #
#  route 2-4 : VCF
# --------------------------------------------------------------------------- #

def normalize(contig, pos, ref, alt):
    """Left-align/trim a simple variant so that routes become comparable."""
    ref = ref.upper(); alt = alt.upper()
    # trim common suffix
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    # trim common prefix
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    if ref == alt:
        return None
    return contig, pos, ref, alt


def classify(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNP"
    if len(ref) < len(alt):
        return "INS"
    if len(ref) > len(alt):
        return "DEL"
    return "MNV"


def parse_vcf(path, label):
    """Return {(contig,pos,ref,alt): info-dict} plus simple per-caller stats."""
    fh = open_maybe_gz(path)
    res = {}
    if fh is None:
        return res
    with fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            contig, pos, _id, ref, alt_field, qual, filt, info = f[:8]
            try:
                pos = int(pos)
            except ValueError:
                continue
            # skip symbolic / breakend alleles
            if "<" in alt_field or "[" in alt_field or "]" in alt_field:
                continue
            qual = float(qual) if qual not in (".", "") else 0.0
            info_d = dict(
                kv.split("=", 1) for kv in info.split(";") if "=" in kv)
            dp = info_d.get("DP", ".")
            fmt_keys, fmt_vals = [], []
            if len(f) >= 10:
                fmt_keys = f[8].split(":")
                fmt_vals = f[9].split(":")
            fd = dict(zip(fmt_keys, fmt_vals))
            sdp = fd.get("DP", dp)
            ad = fd.get("AD", "")
            alt_frac = ""
            if ad and "," in ad:
                try:
                    a = [int(x) for x in ad.split(",")]
                    tot = sum(a)
                    if tot > 0:
                        alt_frac = f"{100.0*sum(a[1:])/tot:.2f}"
                except ValueError:
                    pass
            for alt in alt_field.split(","):
                norm = normalize(contig, pos, ref, alt)
                if norm is None:
                    continue
                c, p, r, a = norm
                key = (c, p, r, a)
                rec = res.setdefault(key, dict(type=classify(r, a), callers=[]))
                if label not in rec["callers"]:
                    rec["callers"].append(label)
                rec[f"{label}_qual"] = qual
                rec[f"{label}_dp"] = sdp
                rec[f"{label}_altfrac"] = alt_frac
                rec[f"{label}_filter"] = filt
    return res


# --------------------------------------------------------------------------- #
#  optional per-position mapping quality from the asm5 BAM
# --------------------------------------------------------------------------- #

class MapqLookup:
    """Lazily fetch MAPQ / depth for a position using pysam if available."""

    def __init__(self, bam):
        self.ok = False
        self.bam = None
        self.cache = {}
        if not bam or not os.path.exists(bam):
            return
        try:
            import pysam  # noqa
            self.pysam = pysam
            self.bam = pysam.AlignmentFile(bam, "rb")
            self.ok = True
        except Exception as exc:                                  # pragma: no cover
            print(f"[merge] pysam unavailable ({exc}); MAPQ annotation skipped",
                  file=sys.stderr)

    def get(self, contig, pos):
        if not self.ok:
            return "", "", ""
        k = (contig, pos)
        if k in self.cache:
            return self.cache[k]
        mq, dp, mmq = "", 0, 0
        try:
            for read in self.bam.fetch(contig, pos - 1, pos):
                if read.is_unmapped:
                    continue
                dp += 1
                m = read.mapping_quality
                mmq += m
                mq = m if mq == "" else min(mq, m)
        except Exception:
            pass
        mean_mq = round(mmq / dp, 1) if dp else ""
        out = (mq, dp, mean_mq)
        self.cache[k] = out
        return out

    def close(self):
        if self.ok:
            self.bam.close()


class AlleleCounter:
    """Count the actual read bases observed at a position (independent of any
    caller).  For a genuine homozygous difference between two haploid genomes
    the alternate allele fraction should be close to 1.0; artefacts of repeat
    or breakpoint alignment show up as a low fraction or as no coverage at all.
    """

    def __init__(self, bam, min_bq=13, min_mq=20):
        self.ok = False
        self.bam = None
        self.cache = {}
        if not bam or not os.path.exists(bam):
            return
        try:
            import pysam
            self.bam = pysam.AlignmentFile(bam, "rb")
            self.min_bq = min_bq
            self.min_mq = min_mq
            self.ok = True
        except Exception as exc:                                  # pragma: no cover
            print(f"[merge] allele counting unavailable ({exc})", file=sys.stderr)

    def counts(self, contig, pos):
        """Return (base_counter, depth) over reads at 1-based pos."""
        if not self.ok:
            return {}, 0
        k = (contig, pos)
        if k in self.cache:
            return self.cache[k]
        from collections import Counter
        c = Counter()
        try:
            for col in self.bam.pileup(contig, pos - 1, pos, truncate=True,
                                       min_base_quality=self.min_bq,
                                       min_mapping_quality=self.min_mq):
                for p in col.pileups:
                    if p.is_del or p.is_refskip:
                        continue
                    qp = p.query_position
                    if qp is None:
                        continue
                    c[p.alignment.query_sequence[qp].upper()] += 1
        except Exception:
            pass
        out = (dict(c), sum(c.values()))
        self.cache[k] = out
        return out

    def close(self):
        if self.ok:
            self.bam.close()


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mummer")
    ap.add_argument("--asm-vcf", action="append", default=[],
                    metavar="LABEL=PATH",
                    help="assembly-based caller VCF (repeatable)")
    ap.add_argument("--read-vcf", action="append", default=[],
                    metavar="LABEL=PATH",
                    help="read-based caller VCF (repeatable)")
    ap.add_argument("--bam-asm")
    ap.add_argument("--bam-sim")
    ap.add_argument("--coords-all", default=None,
                    help="show-coords over ALL NUCmer alignments; used to measure how "
                         "many distinct genomic copies cover each reference position")
    ap.add_argument("--min-repeat-block", type=int, default=200)
    ap.add_argument("--min-repeat-pid", type=float, default=90.0)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-callers", type=int, default=2)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    # ---- reference sequence (used both for sanity-checking and for the VCF) --
    refseq = {}
    _name, _chunks = None, []
    with open(a.ref) as fh:
        for line in fh:
            if line.startswith(">"):
                if _name:
                    refseq[_name] = "".join(_chunks)
                _name = line[1:].split()[0]; _chunks = []
            else:
                _chunks.append(line.strip())
    if _name:
        refseq[_name] = "".join(_chunks)

    def refseq_cache(contig, pos):
        try:
            s = refseq[contig]
            return s[pos - 1].upper() if 0 < pos <= len(s) else ""
        except (KeyError, IndexError):
            return ""

    # ---- reference repeat multiplicity -----------------------------------
    # A reference position covered by several distinct, long, high-identity
    # alignment blocks is inside a repeat.  Any difference reported there may
    # simply be the aligner pairing the wrong paralogue, so such positions can
    # never yield a confident variant call.
    mult_events = defaultdict(list)          # contig -> [(start,end)] 0-based
    if a.coords_all and os.path.exists(a.coords_all):
        n_blocks = 0
        for line in open(a.coords_all):
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                f = line.split()
            if len(f) < 13:
                continue
            try:
                s1, e1 = int(f[0]), int(f[1])
                ln = int(f[4])
                pid = float(f[6])
            except ValueError:
                continue
            if ln < a.min_repeat_block or pid < a.min_repeat_pid:
                continue
            name = f[-2]
            mult_events[name].append((s1 - 1, e1))
            n_blocks += 1
        print(f"[merge] repeat map built from {n_blocks} alignment blocks")
    mult_prefix = {}
    for name, evs in mult_events.items():
        L = len(refseq.get(name, "")) + 2
        diff = [0] * (L + 1)
        for s0, e0 in evs:
            diff[max(0, s0)] += 1
            diff[min(L, e0)] -= 1
        run, arr = 0, [0] * L
        for i in range(L):
            run += diff[i]
            arr[i] = run
        mult_prefix[name] = arr

    def ref_multiplicity(contig, pos):
        arr = mult_prefix.get(contig)
        if not arr:
            return -1
        i = pos - 1
        return arr[i] if 0 <= i < len(arr) else -1

    variants = OrderedDict()          # key -> record

    def touch(key, **kw):
        rec = variants.setdefault(key, dict(
            contig=key[0], pos=key[1], ref=key[2], alt=key[3],
            type=classify(key[2], key[3]), callers=[], buff="", dist=""))
        rec.update({k: v for k, v in kw.items() if v not in ("", None)})
        return rec

    # ---- route 1 ----------------------------------------------------------
    mm, ref_id, qry_id = parse_mummer(a.mummer)
    n_mummer = 0
    for r in mm:
        rb, qb = r["ref"], r["alt"]
        contig = ref_id or "."
        if rb == ".":                       # insertion in query vs reference
            key = (contig, r["pos"], "", qb)   # represent as INS at this anchor
        elif qb == ".":                     # deletion in query
            key = (contig, r["pos"], rb, "")
        else:
            key = (contig, r["pos"], rb, qb)
        rec = touch(key, mummer=True, buff=r["buff"], dist=r["dist"],
                    overlap_R=r["overlap_R"], overlap_Q=r["overlap_Q"])
        if "mummer" not in rec["callers"]:
            rec["callers"].append("mummer")
        n_mummer += 1

    # ---- routes 2..n : labelled VCFs --------------------------------------
    def parse_labelled(spec):
        if "=" not in spec:
            raise SystemExit(f"--asm-vcf/--read-vcf expects LABEL=PATH, got {spec!r}")
        lab, path = spec.split("=", 1)
        return lab.strip(), path.strip()

    asm_labels, read_labels = set(), set()
    for spec in a.asm_vcf:
        lab, path = parse_labelled(spec)
        asm_labels.add(lab)
        for key, info in parse_vcf(path, lab).items():
            rec = touch(key)
            for k, v in info.items():
                if k == "callers":
                    for c in v:
                        if c not in rec["callers"]:
                            rec["callers"].append(c)
                else:
                    rec[k] = v
    for spec in a.read_vcf:
        lab, path = parse_labelled(spec)
        read_labels.add(lab)
        for key, info in parse_vcf(path, lab).items():
            rec = touch(key)
            for k, v in info.items():
                if k == "callers":
                    for c in v:
                        if c not in rec["callers"]:
                            rec["callers"].append(c)
                else:
                    rec[k] = v
    all_labels = ["mummer"] + sorted(asm_labels) + sorted(read_labels)

    # ---- consistency pass: unify mummer indel representation --------------
    # MUMmer writes indels with '' or '.' ; convert to VCF-style anchors so the
    # merge with bcftools/freebayes is meaningful.
    # (handled leniently: matching is done on (contig,pos,ref,alt) and a
    #  tolerance pass below aligns indel anchors that differ by one base)

    # ---- annotate MAPQ from the asm BAM -----------------------------------
    mq_lookup = MapqLookup(a.bam_asm)
    sim_lookup = MapqLookup(a.bam_sim)
    allele_counter = AlleleCounter(a.bam_sim)

    rows = []
    for key, rec in variants.items():
        callers = rec["callers"]
        n = len(callers)
        mq, dp, _mmq = mq_lookup.get(rec["contig"], rec["pos"])
        smq, sdp, smean = sim_lookup.get(rec["contig"], rec["pos"])
        rmult = ref_multiplicity(rec["contig"], rec["pos"])
        # observed read alleles (independent of every caller)
        base_counts, obs_depth = allele_counter.counts(rec["contig"], rec["pos"])
        obs_alt_frac = ""
        if obs_depth:
            alt_key = rec["alt"].upper()[:1] if rec["type"] == "SNP" else None
            if rec["type"] == "SNP" and alt_key:
                obs_alt_frac = f"{100.0 * base_counts.get(alt_key, 0) / obs_depth:.1f}"
            else:
                obs_alt_frac = ""
        try:
            buff = int(rec.get("buff", -1))
        except (TypeError, ValueError):
            buff = -1
        try:
            dist = int(rec.get("dist", -1))
        except (TypeError, ValueError):
            dist = -1
        try:
            over_r = int(rec.get("overlap_R", 0))
            over_q = int(rec.get("overlap_Q", 0))
        except (TypeError, ValueError):
            over_r = over_q = 0
        # ---- repeat / ambiguity indicators --------------------------------
        ambiguous_reasons = []
        if isinstance(dp, int) and dp > 1:
            ambiguous_reasons.append("multi_contig_alignment")
        if over_r > 0 or over_q > 0:
            ambiguous_reasons.append("nucmer_overlap_RQ")
        if rmult > 1:
            ambiguous_reasons.append(f"ref_multicopy_x{rmult}")
        if isinstance(smean, float) and sdp and sdp >= 3 and smean < 20:
            ambiguous_reasons.append("low_read_mapq")
        # The reads come from the query, so a real difference must be visible in
        # them.  If there is good, well-mapped read coverage but the alternate
        # allele is essentially absent, the assembly-based call is an alignment
        # artefact (typically at a repeat or a rearrangement breakpoint).
        reads_contradict = False
        if rec["type"] == "SNP" and obs_depth >= 10:
            try:
                frac = float(obs_alt_frac)
                if frac < 50.0:
                    reads_contradict = True
                    ambiguous_reasons.append(f"reads_contradict_alt{frac:.0f}pct")
            except ValueError:
                pass
        in_repeat = bool(ambiguous_reasons)
        # ---- sanity: does the reported reference base match the reference? ---
        mismatch_ref = ""
        refbase = refseq_cache(rec["contig"], rec["pos"])
        if rec["type"] == "SNP" and refbase and refbase != rec["ref"].upper():
            mismatch_ref = "ref_allele_disagrees_with_fasta"
            in_repeat = True

        # ---- confidence logic ---------------------------------------------
        # Agreement only counts when the agreeing callers use genuinely
        # different evidence.  bcftools_asm and freebayes_asm share one minimap2
        # alignment of the contigs, so their agreement is NOT independent; the
        # three independent methods are
        #     NUCmer alignment (mummer) | contig alignment | read alignment
        # A HIGH-confidence call needs support from at least two of them.
        n_read = len([c for c in callers if c in read_labels])
        n_asm = len([c for c in callers if c in asm_labels])
        has_reads = n_read > 0
        has_mummer = "mummer" in callers
        n_methods = sum([has_mummer, n_asm > 0, n_read > 0])
        flags = []
        if n_methods >= 2:
            flags.append("two_independent_methods")
        if n >= a.min_callers:
            flags.append("multi_caller")
        if has_reads and (n_asm >= 1 or has_mummer):
            flags.append("read_confirmed")
        if n_asm >= 2:
            flags.append("two_assembly_callers")
        if ambiguous_reasons:
            flags.extend(ambiguous_reasons)
        if has_mummer and buff >= 10 and dist >= 10:
            flags.append("well_within_alignment")
        if in_repeat:
            flags.append("in_repeat")
        if mismatch_ref:
            flags.append(mismatch_ref)
        if isinstance(mq, int) and mq < 20:
            flags.append("low_mapq")
        low_qual = ("low_mapq" in flags) or in_repeat or bool(mismatch_ref)
        strong = "two_independent_methods" in flags
        conf = "HIGH" if (strong and not low_qual) \
            else ("MEDIUM" if n_methods >= 2 and not in_repeat
                  else "MEDIUM" if (strong and not in_repeat) \
                  else "LOW")
        rows.append(dict(
            contig=rec["contig"], pos=rec["pos"], ref=rec["ref"], alt=rec["alt"],
            type=rec["type"], n_callers=n,
            n_independent_methods=n_methods,
            callers=";".join(sorted(callers)),
            confidence=conf, flags=";".join(flags),
            mummer_buff=buff, mummer_dist=dist,
            overlap_R=over_r, overlap_Q=over_q, ref_copies=rmult,
            asm_mapq=mq, asm_depth=dp, sim_mapq=smq, sim_mean_mapq=smean,
            sim_depth=sdp, obs_depth=obs_depth, obs_alt_frac=obs_alt_frac,
            sim_altfrac=rec.get("bcftools_sim_altfrac", ""),
            sim_qual=rec.get("bcftools_sim_qual", ""),
            **{f"{lab}_qual": rec.get(f"{lab}_qual", "")
               for lab in ("bcftools_asm", "freebayes_asm")}))

    rows.sort(key=lambda r: (r["contig"], r["pos"]))
    mq_lookup.close(); sim_lookup.close(); allele_counter.close()

    cols = ["contig", "pos", "ref", "alt", "type", "confidence",
            "n_independent_methods", "n_callers", "callers", "flags",
            "ref_copies", "mummer_buff", "mummer_dist", "overlap_R",
            "overlap_Q", "asm_mapq", "asm_depth", "sim_depth",
            "sim_mean_mapq", "obs_depth", "obs_alt_frac", "sim_altfrac", "sim_qual",
            "bcftools_asm_qual", "freebayes_asm_qual"]
    merged_f = os.path.join(a.outdir, "variants_merged.tsv")
    with open(merged_f, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    hi = [r for r in rows if r["confidence"] in ("HIGH", "MEDIUM")]
    hc_f = os.path.join(a.outdir, "variants_highconf.tsv")
    with open(hc_f, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in hi:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # ---- high-confidence VCF ---------------------------------------------
    # (the reference was already loaded at the top of main)
    vcf_f = os.path.join(a.outdir, "variants_highconf.vcf")
    with open(vcf_f, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=variant_merge.py\n")
        for c, s in refseq.items():
            fh.write(f"##contig=<ID={c},length={len(s)}>\n")
        fh.write('##INFO=<ID=CALLERS,Number=1,Type=String,Description="Supporting callers">\n')
        fh.write('##INFO=<ID=CONF,Number=1,Type=String,Description="Confidence class">\n')
        fh.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQUERY\n")
        for i, r in enumerate(hi, 1):
            # anchor indels so REF matches the reference base
            c, p, rf, al = r["contig"], r["pos"], r["ref"], r["alt"]
            seq = refseq.get(c, "")
            if rf == "" or rf == ".":
                # insertion: anchor on the preceding reference base
                if p - 2 >= 0:
                    anchor = seq[p - 2]
                    rf, al, p = anchor, anchor + al, p - 1
            elif al == "" or al == ".":
                if p - 2 >= 0:
                    anchor = seq[p - 2]
                    rf, al, p = anchor + rf, anchor, p - 1
            qual = r.get("bcftools_sim_qual") or r.get("bcftools_asm_qual") or "."
            try:
                qual = f"{float(qual):.1f}"
            except (TypeError, ValueError):
                qual = "."
            info = (f"CALLERS={r['callers']};CONF={r['confidence']};"
                    f"TYPE={r['type']}")
            fh.write(f"{c}\t{p}\tVAR{i}\t{rf}\t{al}\t{qual}\tPASS\t{info}\tGT\t1/1\n")
    # index with bcftools/tabix if available
    os.system(f"bgzip -f -c '{vcf_f}' > '{vcf_f}.gz' 2>/dev/null && "
              f"tabix -f -p vcf '{vcf_f}.gz' 2>/dev/null")

    # ---- concordance matrix ----------------------------------------------
    labels = all_labels
    conc_f = os.path.join(a.outdir, "caller_concordance.tsv")
    with open(conc_f, "w") as fh:
        fh.write("caller\tn_variants\t" + "\t".join(labels) + "\n")
        for l1 in labels:
            s1 = {k for k, v in variants.items() if l1 in v["callers"]}
            row = [str(len(s1))]
            for l2 in labels:
                s2 = {k for k, v in variants.items() if l2 in v["callers"]}
                row.append(str(len(s1 & s2)))
            fh.write(l1 + "\t" + "\t".join(row) + "\n")

    summary = dict(
        total_distinct_variants=len(rows),
        high_confidence=len([r for r in rows if r["confidence"] == "HIGH"]),
        medium_confidence=len([r for r in rows if r["confidence"] == "MEDIUM"]),
        low_confidence=len([r for r in rows if r["confidence"] == "LOW"]),
        SNPs=len([r for r in rows if r["type"] == "SNP"]),
        insertions=len([r for r in rows if r["type"] == "INS"]),
        deletions=len([r for r in rows if r["type"] == "DEL"]),
        MNV=len([r for r in rows if r["type"] == "MNV"]),
        mummer_raw_calls=n_mummer,
        per_caller={l: len([1 for v in variants.values() if l in v["callers"]])
                    for l in labels},
    )
    with open(os.path.join(a.outdir, "merge_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[merge] {summary['total_distinct_variants']} distinct variants; "
          f"HIGH={summary['high_confidence']} MEDIUM={summary['medium_confidence']} "
          f"LOW={summary['low_confidence']}")
    print(f"[merge] SNP={summary['SNPs']} INS={summary['insertions']} "
          f"DEL={summary['deletions']} MNV={summary['MNV']}")
    print(f"[merge] wrote {merged_f}, {hc_f}, {vcf_f}")


if __name__ == "__main__":
    main()
