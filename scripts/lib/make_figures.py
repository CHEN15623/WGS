#!/usr/bin/env python3
"""
make_figures.py -- every figure in the report.

All plots are pure matplotlib so no R/Bioconductor stack is required.  Each
figure is written as PNG (300 dpi) and, where it helps, also as a vector PDF.

Figure map
  fig01_assembly_qc          assembly statistics, cumulative contig curve
  fig02_read_qc              simulated-read length/quality summary
  fig03_synteny_dotplot      reference vs query collinearity (NUCmer blocks)
  fig04_structural_variants  linear rearrangement map: origin shift + inversion
  fig05_genome_circos        circular overview: GC, GC skew, CDS, variants, depth
  fig06_coverage             depth profile with uncovered regions
  fig07_variant_landscape    all candidates vs the verified set
  fig08_variant_confidence   why candidates were accepted or rejected
  fig09_substitution_spectrum  SNP substitution classes
  fig10_variant_verification sequence-level view of each verified variant
  fig11_gene_impact          ORF consequence per affected gene
  fig12_functional_profile   COG categories, reference vs query
  fig13_amr_determinants     resistance / virulence / metabolism determinant map
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict, OrderedDict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Wedge, FancyArrow
from matplotlib.lines import Line2D

plt.rcParams.update({
    "figure.dpi": 110,
    "savefig.dpi": 300,
    "font.size": 9,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.autolayout": False,
})

# ---------------------------------------------------------------- palette
C_REF = "#2166ac"
C_QRY = "#b2182b"
C_OK = "#1b7837"
C_WARN = "#d95f02"
C_BAD = "#999999"
CONF_COLOR = {"HIGH": "#1b7837", "MEDIUM": "#f4a582", "LOW": "#bbbbbb"}
IMPACT_COLOR = {"HIGH": "#b2182b", "MODERATE": "#ef8a62", "LOW": "#67a9cf",
                "MODIFIER": "#cccccc"}


# ====================================================================== utils
def read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    rows = []
    with open(path, errors="replace") as fh:
        hdr = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            f += [""] * (len(hdr) - len(f))
            rows.append(dict(zip(hdr, f)))
    return rows


def read_fasta(path):
    seqs, name, chunks = OrderedDict(), None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name:
        seqs[name] = "".join(chunks)
    return seqs


def save(fig, outdir, name, pdf=True):
    png = os.path.join(outdir, name + ".png")
    fig.savefig(png, bbox_inches="tight", facecolor="white")
    if pdf:
        try:
            fig.savefig(os.path.join(outdir, name + ".pdf"),
                        bbox_inches="tight", facecolor="white")
        except Exception:
            pass
    plt.close(fig)
    print(f"  [fig] {png}")


def load_json(path):
    if path and os.path.exists(path):
        try:
            return json.load(open(path))
        except Exception:
            return {}
    return {}


def load_cds(gff, types=("CDS",)):
    out = []
    if not gff or not os.path.exists(gff):
        return out
    for line in open(gff, errors="replace"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in types:
            continue
        d = {}
        for kv in f[8].split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                d.setdefault(k, v)
        out.append(dict(contig=f[0], start=int(f[3]), end=int(f[4]), strand=f[6],
                        type=f[2], locus=d.get("locus_tag", ""),
                        gene=d.get("gene", d.get("Name", "")),
                        product=d.get("product", "")))
    return out


def gc_windows(seq, win):
    n = max(1, len(seq) // win)
    gc, skew = np.zeros(n), np.zeros(n)
    for i in range(n):
        s = seq[i * win:(i + 1) * win].upper()
        g, c = s.count("G"), s.count("C")
        a, t = s.count("A"), s.count("T")
        gc[i] = 100.0 * (g + c) / max(1, g + c + a + t)
        skew[i] = (g - c) / (g + c) if (g + c) else 0.0
    return gc, skew


def read_depth(path, contig):
    """Read a `samtools depth` table; transparently handles .gz and falls back
    to the compressed name when only that exists."""
    if path and not os.path.exists(path) and os.path.exists(path + ".gz"):
        path = path + ".gz"
    pos, dep = [], []
    if not path or not os.path.exists(path):
        return np.array(pos), np.array(dep)
    import gzip
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0] != contig:
                continue
            try:
                pos.append(int(f[1])); dep.append(int(f[2]))
            except ValueError:
                continue
    return np.array(pos), np.array(dep)


# ================================================================== fig 01 QC
def fig_assembly_qc(ctx, outdir):
    r = load_json(os.path.join(ctx["out"], "01_qc/assembly_qc",
                               f"reference_{ctx['ref_name']}_assembly_stats.json"))
    q = load_json(os.path.join(ctx["out"], "01_qc/assembly_qc",
                               f"query_{ctx['qry_name']}_assembly_stats.json"))
    if not r or not q:
        return
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))

    # (a) headline metrics
    ax = axes[0]
    keys = ["total_length_bp", "N50_bp", "largest_contig_bp", "num_contigs",
            "GC_percent", "ambiguous_bases"]
    labels = ["total length\n(Mb)", "N50\n(Mb)", "largest contig\n(Mb)",
              "contigs", "GC %", "ambiguous\nbases"]
    rv = [r["total_length_bp"] / 1e6, r["N50_bp"] / 1e6,
          r["largest_contig_bp"] / 1e6, r["num_contigs"], r["GC_percent"],
          r["ambiguous_bases"]]
    qv = [q["total_length_bp"] / 1e6, q["N50_bp"] / 1e6,
          q["largest_contig_bp"] / 1e6, q["num_contigs"], q["GC_percent"],
          q["ambiguous_bases"]]
    x = np.arange(len(keys))
    ax.bar(x - 0.2, rv, 0.4, label=f"reference {ctx['ref_name']}", color=C_REF)
    ax.bar(x + 0.2, qv, 0.4, label=f"query {ctx['qry_name']}", color=C_QRY)
    for xi, (a, b) in enumerate(zip(rv, qv)):
        for xo, v in ((-0.2, a), (0.2, b)):
            ax.text(xi + xo, v, f"{v:g}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_yscale("log")
    ax.set_ylabel("value (log scale)")
    ax.set_title("(a) assembly quality metrics")
    ax.legend(fontsize=7.5, frameon=False)

    # (b) cumulative contig length
    ax = axes[1]
    for label, stats_name, col in ((f"reference {ctx['ref_name']}",
                                    f"reference_{ctx['ref_name']}", C_REF),
                                   (f"query {ctx['qry_name']}",
                                    f"query_{ctx['qry_name']}", C_QRY)):
        rows = read_tsv(os.path.join(ctx["out"], "01_qc/assembly_qc",
                                     f"{stats_name}_cumulative_length.tsv"))
        if not rows:
            continue
        xs = [1] + [int(x["contig_rank"]) for x in rows]
        ys = [0] + [float(x["cumulative_percent"]) for x in rows]
        ax.step(xs, ys, where="post", label=label, color=col, lw=2)
    ax.axhline(50, ls=":", color="k", lw=1)
    ax.text(0.98, 52, "N50", ha="right", fontsize=7, transform=ax.get_yaxis_transform())
    ax.set_xlabel("contig rank"); ax.set_ylabel("cumulative length (%)")
    ax.set_title("(b) cumulative contig curve")
    ax.set_ylim(0, 105); ax.legend(fontsize=7.5, frameon=False)

    # (c) annotation content
    ax = axes[2]
    rows = read_tsv(os.path.join(ctx["out"], "01_qc/annotation_qc",
                                 "reference_annotation_summary.tsv"))
    if rows:
        keep = [x for x in rows if x["feature_type"] in
                ("CDS", "tRNA", "rRNA", "tmRNA", "pseudogene", "gene", "exon")][:7]
        ax.barh([x["feature_type"] for x in keep][::-1],
                [float(x["count"]) for x in keep][::-1], color=C_REF)
        for i, x in enumerate(keep[::-1]):
            ax.text(float(x["count"]), i, f" {int(float(x['count']))}",
                    va="center", fontsize=7.5)
    ax.set_xlabel("features"); ax.set_title("(c) reference annotation content")
    fig.suptitle("Figure 1  Assembly and annotation quality control", y=1.03,
                 fontsize=12, fontweight="bold")
    save(fig, outdir, "fig01_assembly_qc")


# ============================================================== fig 02 read QC
def fig_read_qc(ctx, outdir):
    fp = load_json(os.path.join(ctx["out"], "01_qc/reads_qc/fastp.json"))
    stats = os.path.join(ctx["out"], "01_qc/reads_qc")
    files = [f for f in os.listdir(stats) if f.endswith(".tsv")] if os.path.isdir(stats) else []
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))

    ax = axes[0]
    if fp:
        s = fp.get("summary", {})
        bef = s.get("before_filtering", {})
        aft = s.get("after_filtering", {})
        keys = ["total_reads", "total_bases", "q30_rate", "gc_content"]
        labels = ["reads", "bases", "Q30 rate", "GC content"]
        bv, av = [], []
        for k in keys:
            bv.append(bef.get(k, 0)); av.append(aft.get(k, 0))
        x = np.arange(len(keys))
        ax.bar(x - 0.2, [bv[0], bv[1] / 1e6, bv[2] * 100, bv[3] * 100], 0.4,
               label="before", color="#cccccc")
        ax.bar(x + 0.2, [av[0], av[1] / 1e6, av[2] * 100, av[3] * 100], 0.4,
               label="after filtering", color=C_OK)
        ax.set_xticks(x)
        ax.set_xticklabels(["reads", "bases (Mb)", "Q30 (%)", "GC (%)"], fontsize=8)
        ax.legend(fontsize=8, frameon=False)
    ax.set_title("(a) simulated-read statistics (fastp)")
    ax.set_ylabel("value")

    ax = axes[1]
    seqkit = None
    for f in files:
        if "seqkit" in f:
            seqkit = read_tsv(os.path.join(stats, f)); break
    if not seqkit:
        # fall back to the pipeline-standard location
        for cand in os.listdir(ctx["out"]):
            pass
    if seqkit:
        rows = seqkit
        labels = [f"R1" if i == 0 else f"R2" for i in range(len(rows))]
        nums = [float(r.get("num_seqs", 0) or 0) / 1e6 for r in rows]
        lens = [float(r.get("avg_len", 0) or 0) for r in rows]
        ax.bar(labels, nums, color=C_QRY)
        for i, (n, l) in enumerate(zip(nums, lens)):
            ax.text(i, n, f"{n:.2f}M reads\n{l:.0f} bp", ha="center",
                    va="bottom", fontsize=8)
        ax.set_ylabel("read pairs (millions)")
    else:
        ax.text(0.5, 0.5, "seqkit summary not found", ha="center", va="center",
                transform=ax.transAxes)
    ax.set_title("(b) read count and length")
    fig.suptitle("Figure 2  Quality control of the in-silico reads used for "
                 "variant confirmation", y=1.04, fontsize=11, fontweight="bold")
    save(fig, outdir, "fig02_read_qc")


# ========================================================= fig 03 dot plot
def fig_synteny_dotplot(ctx, outdir):
    rows = read_tsv(os.path.join(ctx["out"], "02_alignment/synteny",
                                 "nucmer_coords_1to1.tsv"))
    if not rows:
        rows = read_tsv(os.path.join(ctx["out"], "02_alignment/synteny",
                                     "alignment_blocks.tsv"))
    blocks = []
    for r in rows:
        try:
            if "ref_start" in r:      # alignment_blocks.tsv layout
                rs, re_, qs, qe = (int(r["ref_start"]), int(r["ref_end"]),
                                   int(r["qry_start"]), int(r["qry_end"]))
                strand = r["strand"]
            else:                     # show-coords layout
                rs, re_, qs, qe = (int(r["NC_003228.3"]) if "NC_003228.3" in r
                                   else 0, 0, 0, 0)
                continue
        except (ValueError, KeyError):
            continue
        blocks.append((rs, re_, qs, qe, strand))
    if not rows:
        return
    # show-coords output is positional; re-parse it directly
    blocks = []
    p = os.path.join(ctx["out"], "02_alignment/synteny/nucmer_coords_1to1.tsv")
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            f = line.split()
        if len(f) < 13:
            continue
        try:
            s1, e1, s2, e2 = int(f[0]), int(f[1]), int(f[2]), int(f[3])
            qdir = int(f[-3])
        except ValueError:
            continue
        if f[-2] != ctx.get("ref_chr", "NC_003228.3"):
            continue
        blocks.append((s1, e1, s2, e2, "+" if qdir >= 0 else "-"))

    fig, ax = plt.subplots(figsize=(7.5, 7))
    for s1, e1, s2, e2, strand in blocks:
        col = C_OK if strand == "+" else C_QRY
        ax.plot([s1, e1], [s2, e2], color=col, lw=4 if (e1 - s1) > 1e5 else 2.2,
                alpha=0.85, solid_capstyle="butt")
    sv = read_tsv(os.path.join(ctx["out"], "02_alignment/synteny",
                               "structural_variants.tsv"))
    for r in sv:
        if r.get("type") == "inversion":
            try:
                ax.annotate("88 kb inversion",
                            xy=(int(r["ref_start"]), int(r["qry_end"])),
                            xytext=(0.35, 0.45), textcoords="axes fraction",
                            arrowprops=dict(arrowstyle="->", color="k", lw=1),
                            fontsize=9, fontweight="bold")
            except (ValueError, KeyError):
                pass
    ax.set_xlabel(f"reference {ctx['ref_name']} position (bp)")
    ax.set_ylabel(f"query {ctx['qry_name']} position (bp)")
    ax.set_title("Figure 3  Whole-genome alignment: collinearity of the two "
                 "chromosomes")
    ax.xaxis.set_major_formatter(lambda v, p: f"{v/1e6:.1f}M")
    ax.yaxis.set_major_formatter(lambda v, p: f"{v/1e6:.1f}M")
    ax.legend(handles=[Line2D([], [], color=C_OK, lw=3, label="forward block"),
                       Line2D([], [], color=C_QRY, lw=3, label="inverted block")],
              fontsize=8, loc="upper left", frameon=False)
    ax.text(0.02, 0.02,
            "the query is the reference re-linearised at a different origin;\n"
            "the only true rearrangement is an 88,340 bp inversion",
            transform=ax.transAxes, fontsize=7.5, style="italic", color="#444")
    ax.set_aspect("equal", adjustable="datalim")
    save(fig, outdir, "fig03_synteny_dotplot")


# ================================================= fig 04 structural variants
def fig_structural_variants(ctx, outdir):
    summ = load_json(os.path.join(ctx["out"], "02_alignment/synteny",
                                  "structural_summary.json"))
    origin = read_tsv(os.path.join(ctx["out"], "02_alignment/synteny",
                                   "origin_shift.tsv"))
    sv = read_tsv(os.path.join(ctx["out"], "02_alignment/synteny",
                               "structural_variants.tsv"))
    if not sv:
        return
    ref_len = 5205140
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(13, 5.2),
                                  gridspec_kw={"height_ratios": [1.25, 1]})

    # ---- rearranged reference segments --------------------------------
    palette = {"inversion": C_QRY, "whole_replicon_reverse_orientation": "#f4a582",
               "small_reverse_block": "#9ecae1", "inverted_repeat": "#fdae6b"}
    seen = set()
    for r in sv:
        if r["ref_contig"] != ctx.get("ref_chr", "NC_003228.3"):
            continue
        try:
            s, e = int(r["ref_start"]), int(r["ref_end"])
        except ValueError:
            continue
        t = r["type"]
        if s > e:
            s, e = e, s
        ax.add_patch(Rectangle((s, 0.42), max(1, e - s), 0.18,
                               color=palette.get(t, "#cccccc"),
                               label=None if t in seen else t))
        if (e - s) > 40000:
            ax.text((s + e) / 2, 0.51, f"{t.replace('_',' ')}\n{e-s:,} bp",
                    ha="center", va="center", fontsize=7.5, fontweight="bold")
        seen.add(t)
    ax.add_patch(Rectangle((0, 0.14), ref_len, 0.14, color="#dddddd"))
    ax.text(ref_len * 0.5, 0.21, "reference chromosome (collinear backbone)",
            ha="center", va="center", fontsize=8, color="#333")
    off = origin[0] if origin else {}
    txt = ("origin shift (deposited-sequence convention, NOT a biological "
           "difference):\nreference position 1 corresponds to query position "
           f"{int(off.get('dominant_offset_bp', 0)):,}")
    ax.text(0.008, 0.80, txt, transform=ax.transAxes, fontsize=8, va="top",
            bbox=dict(fc="#fff8e1", ec="#f0c36d", boxstyle="round,pad=0.4"))
    inv = [r for r in sv if r["type"] == "inversion"]
    if inv:
        i0 = int(inv[0]["ref_start"]); i1 = int(inv[0]["ref_end"])
        ax.annotate(f"inversion {i1-i0:,} bp", xy=((i0 + i1) / 2, 0.42),
                    xytext=((i0 + i1) / 2, 0.80), ha="center",
                    arrowprops=dict(arrowstyle="->", color=C_QRY, lw=1.4),
                    fontsize=9, color=C_QRY, fontweight="bold")
    ax.set_xlim(-ref_len * 0.02, ref_len * 1.02)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.xaxis.set_major_formatter(lambda v, p: f"{v/1e6:.0f} Mb")
    ax.set_xlabel("reference position")
    ax.set_title("Figure 4  Structural differences between the two chromosomes")
    ax.grid(False)

    # ---- alignment blocks as a query-vs-reference ribbon ---------------
    blocks = read_tsv(os.path.join(ctx["out"], "02_alignment/synteny",
                                   "alignment_blocks.tsv"))
    for b in blocks:
        try:
            rs, re_ = int(b["ref_start"]), int(b["ref_end"])
            qs, qe = int(b["qry_start"]), int(b["qry_end"])
        except (ValueError, KeyError):
            continue
        col = C_OK if b["strand"] == "+" else C_QRY
        ax2.plot([rs, re_], [0, 0], color=col, lw=6, solid_capstyle="butt",
                 alpha=0.9)
    ax2.set_xlim(0, ref_len); ax2.set_yticks([])
    ax2.set_xlabel("reference position")
    ax2.xaxis.set_major_formatter(lambda v, p: f"{v/1e6:.0f} Mb")
    ax2.set_title("1-to-1 alignment blocks along the reference "
                  "(green = same orientation, red = inverted)", fontsize=9)
    ax2.grid(False)
    fig.tight_layout()
    save(fig, outdir, "fig04_structural_variants")


# ============================================================ fig 05 circos
def _pol(angle_deg, r):
    a = np.radians(angle_deg)
    return r * np.cos(a), r * np.sin(a)


def fig_circos(ctx, outdir):
    ref = read_fasta(ctx["ref_fa"])
    contig = ctx.get("ref_chr", "NC_003228.3")
    if contig not in ref:
        return
    seq = ref[contig]
    L = len(seq)
    win = 10000
    gc, skew = gc_windows(seq, win)

    cds = [c for c in load_cds(ctx["ref_gff"]) if c["contig"] == contig]
    vars_conf = read_tsv(os.path.join(ctx["out"], "03_variants/verified",
                                      "variants_confident.tsv"))
    vars_amb = read_tsv(os.path.join(ctx["out"], "03_variants/verified",
                                     "variants_ambiguous.tsv"))
    depth_pos, depth = read_depth(os.path.join(ctx["out"], "02_alignment/coverage",
                                               "contig_depth_asm5.tsv"), contig)
    sim_pos, sim_depth = read_depth(os.path.join(ctx["out"], "02_alignment/coverage",
                                                 "sim_read_depth.tsv"), contig)

    fig = plt.figure(figsize=(13, 12))
    ax = fig.add_axes([0.03, 0.03, 0.94, 0.9])
    ax.set_aspect("equal"); ax.axis("off")

    def ang(pos):
        return 90.0 - 360.0 * (pos / L)

    R = 1.0
    # ---- CDS ring ------------------------------------------------------
    for c in cds:
        a0, a1 = ang(c["start"]), ang(c["end"])
        if a1 > a0:
            a0, a1 = a1, a0
        inward = c["strand"] == "-"
        r0, r1 = (0.90, 0.955) if not inward else (0.845, 0.90)
        th = np.linspace(a0, a1, max(2, int(abs(a1 - a0)) + 2))
        xs = [np.cos(np.radians(t)) * r for t in th for r in (r0, r1)]
        ys = [np.sin(np.radians(t)) * r for t in th for r in (r0, r1)]
        for t in th:
            x0, y0 = _pol(t, r0); x1, y1 = _pol(t, r1)
            ax.plot([x0, x1], [y0, y1], color="#2c7fb8" if not inward else "#7fcdbb",
                    lw=0.6, alpha=0.9)

    # ---- GC content ring ----------------------------------------------
    def track(values, r0, r1, color, vmin=None, vmax=None):
        vals = np.asarray(values, dtype=float)
        vmin = np.nanmin(vals) if vmin is None else vmin
        vmax = np.nanmax(vals) if vmax is None else vmax
        norm = (vals - vmin) / (vmax - vmin) if vmax > vmin else vals * 0
        x, y = [], []
        for i, v in enumerate(norm):
            pos = (i + 0.5) * win
            t = ang(min(pos, L))
            xx, yy = _pol(t, r0 + v * (r1 - r0))
            x.append(xx); y.append(yy)
        # close the circle
        x.append(x[0]); y.append(y[0])
        ax.plot(x, y, color=color, lw=1.1)

    track(gc, 0.72, 0.82, "#238b45", vmin=float(np.nanmin(gc)), vmax=float(np.nanmax(gc)))
    ax.text(0, 0.83, f"GC content  ({np.nanmin(gc):.0f}–{np.nanmax(gc):.0f} %)",
            ha="center", va="bottom", fontsize=7.5, color="#238b45")
    track(np.abs(skew), 0.60, 0.70, "#6a51a3", vmin=0.0, vmax=float(np.nanmax(np.abs(skew))) or 1)
    ax.text(0, 0.71, "GC skew |(G−C)/(G+C)|", ha="center", va="bottom",
            fontsize=7.5, color="#6a51a3")

    # ---- coverage ring -------------------------------------------------
    if len(depth):
        # decimate to ~2000 points
        idx = np.linspace(0, len(depth) - 1, min(2000, len(depth))).astype(int)
        x, y = [], []
        dmax = max(2.0, float(np.percentile(depth, 99)))
        for i in idx:
            t = ang(depth_pos[i])
            v = min(1.0, depth[i] / dmax)
            xx, yy = _pol(t, 0.46 + v * 0.10)
            x.append(xx); y.append(yy)
        x.append(x[0]); y.append(y[0])
        ax.plot(x, y, color="#e6550d", lw=0.8)
        ax.text(0, 0.57, "contig-alignment depth", ha="center", va="bottom",
                fontsize=7.5, color="#e6550d")
    if len(sim_depth):
        idx = np.linspace(0, len(sim_depth) - 1, min(2000, len(sim_depth))).astype(int)
        x, y = [], []
        dmax = max(1.0, float(np.percentile(sim_depth, 99)))
        for i in idx:
            t = ang(sim_pos[i])
            v = min(1.0, sim_depth[i] / dmax)
            xx, yy = _pol(t, 0.32 + v * 0.10)
            x.append(xx); y.append(yy)
        x.append(x[0]); y.append(y[0])
        ax.plot(x, y, color="#756bb1", lw=0.6, alpha=0.9)
        ax.text(0, 0.43, "read depth (SIM)", ha="center", va="bottom",
                fontsize=7.5, color="#756bb1")

    # ---- variant scatter ----------------------------------------------
    r_var = 0.24
    for rows, col, size, lab in ((vars_amb, "#bbbbbb", 6, "candidate (ambiguous)"),
                                 (vars_conf, None, 55, "verified")):
        for v in rows:
            try:
                pos = int(v["pos"])
            except ValueError:
                continue
            if v.get("contig") != contig:
                continue
            t = ang(pos)
            xx, yy = _pol(t, r_var)
            c = CONF_COLOR.get(v.get("confidence", ""), "#333333") if col is None else col
            ax.scatter([xx], [yy], s=size, color=c,
                       edgecolor="k" if col is None else "none",
                       linewidth=0.5, zorder=5, alpha=0.95)

    # ---- labelled verified variants -----------------------------------
    for v in vars_conf:
        try:
            pos = int(v["pos"])
        except ValueError:
            continue
        if v.get("contig") != contig:
            continue
        t = ang(pos)
        xx, yy = _pol(t, r_var - 0.045)
        x2, y2 = _pol(t, r_var - 0.13)
        ax.annotate(f"{v['ref']}>{v['alt']} @{pos:,}", xy=(xx, yy), xytext=(x2, y2),
                    fontsize=7.5, ha="center", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="#555"),
                    bbox=dict(fc="white", ec="#999", alpha=0.85,
                              boxstyle="round,pad=0.2"))

    # ---- inversion arc --------------------------------------------------
    sv = read_tsv(os.path.join(ctx["out"], "02_alignment/synteny",
                               "structural_variants.tsv"))
    for s in sv:
        if s.get("type") != "inversion" or s.get("ref_contig") != contig:
            continue
        try:
            a, b = int(s["ref_start"]), int(s["ref_end"])
        except ValueError:
            continue
        th = np.linspace(ang(a), ang(b), 200)
        x = np.cos(np.radians(th)) * 0.115
        y = np.sin(np.radians(th)) * 0.115
        ax.plot(x, y, color=C_QRY, lw=2.5)
        ax.text(np.mean(x), np.mean(y), f"inversion\n{b-a:,} bp", ha="center",
                va="center", fontsize=8, color=C_QRY, fontweight="bold")

    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.3)
    ax.set_title(f"Figure 5  Circular overview of the {ctx['ref_name']} "
                 f"chromosome ({contig}, {L:,} bp)\n"
                 "outer ring: CDS (blue = forward, teal = reverse); then GC "
                 "content, GC skew and coverage;\ninner dots: variant "
                 "candidates (grey) and the verified variants (green)",
                 fontsize=11)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color=C_OK, markeredgecolor="k",
               label="verified variant"),
        Line2D([], [], marker="o", ls="", color="#bbbbbb", label="ambiguous candidate"),
        Line2D([], [], color="#2c7fb8", lw=3, label="CDS forward"),
        Line2D([], [], color="#7fcdbb", lw=3, label="CDS reverse")],
        loc="lower left", fontsize=8, frameon=False)
    save(fig, outdir, "fig05_genome_circos")


# ============================================================ fig 06 coverage
def fig_coverage(ctx, outdir):
    contig = ctx.get("ref_chr", "NC_003228.3")
    p1, d1 = read_depth(os.path.join(ctx["out"], "02_alignment/coverage",
                                     "contig_depth_asm5.tsv"), contig)
    p2, d2 = read_depth(os.path.join(ctx["out"], "02_alignment/coverage",
                                     "sim_read_depth.tsv"), contig)
    if not len(p1) and not len(p2):
        return
    fig, axes = plt.subplots(2, 1, figsize=(13, 5.6), sharex=True,
                             gridspec_kw={"height_ratios": [1, 2]})
    ax = axes[0]
    if len(p2):
        step = max(1, len(p2) // 4000)
        ax.fill_between(p2[::step] / 1e6, d2[::step], color="#756bb1", alpha=0.7,
                        lw=0)
        ax.set_ylabel("read depth\n(SIM)")
    ax.set_ylim(0, max(10, np.percentile(d2, 99.5) * 1.15) if len(d2) else 10)
    ax.set_title("Figure 6  Coverage of the reference genome")

    ax = axes[1]
    if len(p1):
        step = max(1, len(p1) // 6000)
        ax.fill_between(p1[::step] / 1e6, d1[::step], color="#e6550d", alpha=0.75,
                        lw=0, label="contig alignment")
    ax.set_ylim(0, 3)
    ax.set_ylabel("depth\n(contig)")
    ax.set_xlabel(f"reference {ctx['ref_name']} position (Mb)")

    # mark uncovered regions and verified variants
    unc = read_tsv(os.path.join(ctx["out"], "02_alignment/coverage",
                                "uncovered_regions.tsv"))
    for u in unc:
        if u["contig"] != contig:
            continue
        try:
            s, e = int(u["start"]), int(u["end"])
        except ValueError:
            continue
        if e - s >= 5:
            ax.axvspan(s / 1e6, e / 1e6, color="#d62728", alpha=0.25)
    vc = read_tsv(os.path.join(ctx["out"], "03_variants/verified",
                               "variants_confident.tsv"))
    for v in vc:
        try:
            pos = int(v["pos"]) / 1e6
        except ValueError:
            continue
        ax.axvline(pos, color=C_OK, lw=1.2)
        ax.annotate(f"{v['ref']}>{v['alt']}", xy=(pos, 2.4), rotation=90,
                    fontsize=7, color=C_OK, ha="right", va="top")
    inv = read_tsv(os.path.join(ctx["out"], "02_alignment/synteny",
                                "structural_variants.tsv"))
    for s in inv:
        if s.get("type") == "inversion":
            try:
                ax.axvspan(int(s["ref_start"]) / 1e6, int(s["ref_end"]) / 1e6,
                           color=C_QRY, alpha=0.13)
            except ValueError:
                pass
    ax.legend(handles=[Line2D([], [], color="#e6550d", lw=4, label="depth"),
                       Line2D([], [], color=C_OK, lw=1.4, label="verified variant"),
                       Line2D([], [], color=C_QRY, lw=6, alpha=0.3,
                              label="88 kb inversion")],
              fontsize=8, frameon=False, ncol=3, loc="upper right")
    fig.tight_layout()
    save(fig, outdir, "fig06_coverage")


# ==================================================== fig 07 variant landscape
def fig_variant_landscape(ctx, outdir):
    merged = read_tsv(os.path.join(ctx["out"], "03_variants/merged",
                                   "variants_merged.tsv"))
    conf = read_tsv(os.path.join(ctx["out"], "03_variants/verified",
                                 "variants_confident.tsv"))
    if not merged:
        return
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 4.6),
                                  gridspec_kw={"width_ratios": [3, 1]})
    contig = ctx.get("ref_chr", "NC_003228.3")
    for r in merged:
        if r["contig"] != contig:
            continue
        try:
            pos = int(r["pos"]); mult = int(r.get("ref_copies", "1") or 1)
        except ValueError:
            continue
        c = CONF_COLOR.get(r["confidence"], "#999")
        ax.scatter(pos / 1e6, mult + np.random.uniform(-0.18, 0.18),
                   s=10 if r["confidence"] != "HIGH" else 70, color=c,
                   edgecolor="k" if r["confidence"] == "HIGH" else "none",
                   linewidth=0.6, alpha=0.85, zorder=3)
    for r in conf:
        if r["contig"] != contig:
            continue
        try:
            pos = int(r["pos"])
        except ValueError:
            continue
        ax.annotate(f"{r['ref']}>{r['alt']}", xy=(pos / 1e6, 1.0),
                    xytext=(pos / 1e6, -0.75), fontsize=7.5, ha="center",
                    color=C_OK, fontweight="bold",
                    arrowprops=dict(arrowstyle="->", color=C_OK, lw=0.8))
    inv = read_tsv(os.path.join(ctx["out"], "02_alignment/synteny",
                                "structural_variants.tsv"))
    for s in inv:
        if s.get("type") == "inversion":
            try:
                ax.axvspan(int(s["ref_start"]) / 1e6, int(s["ref_end"]) / 1e6,
                           color=C_QRY, alpha=0.12, zorder=0)
            except ValueError:
                pass
    ax.set_xlabel("reference position (Mb)")
    ax.set_ylabel("number of genomic copies covering the position\n"
                  "(alignment multiplicity)")
    ax.set_title("(a) every variant candidate, coloured by confidence\n"
                 "(shaded band = 88 kb inversion)")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=CONF_COLOR["HIGH"],
                              markeredgecolor="k", label="HIGH"),
                       Line2D([], [], marker="o", ls="", color=CONF_COLOR["LOW"],
                              label="LOW / rejected")],
              fontsize=8, frameon=False, loc="upper right")

    # (b) confidence + reason breakdown
    reason = Counter()
    for r in merged:
        fl = (r.get("flags", "") or "")
        toks = [t for t in fl.split(";") if t and t != "in_repeat"]
        if not toks:
            reason["single method only"] += 1
            continue
        for t in toks:
            key = ("multi-copy reference" if t.startswith("ref_multicopy")
                   else "reads contradict call" if t.startswith("reads_contradict")
                   else "NUCmer overlap (R/Q)" if t.startswith("nucmer_overlap")
                   else "multiple contig alignments" if t.startswith("multi_contig")
                   else "low read mapping quality" if t == "low_mapq"
                   else t)
            reason[key] += 1
    items = reason.most_common(7)
    ax2.barh([k for k, _ in items][::-1], [v for _, v in items][::-1],
             color="#969696")
    for i, (k, v) in enumerate(items[::-1]):
        ax2.text(v, i, f" {v}", va="center", fontsize=7.5)
    ax2.set_xlabel("candidates flagged")
    ax2.set_title("(b) why candidates were rejected")
    fig.suptitle("Figure 7  Variant landscape and candidate filtering",
                 y=1.03, fontsize=12, fontweight="bold")
    fig.tight_layout()
    save(fig, outdir, "fig07_variant_landscape")


# ================================================= fig 08 variant confidence
def fig_variant_confidence(ctx, outdir):
    summ = load_json(os.path.join(ctx["out"], "03_variants/merged",
                                  "merge_summary.json"))
    conc = read_tsv(os.path.join(ctx["out"], "03_variants/merged",
                                 "caller_concordance.tsv"))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    ax = axes[0]
    if summ:
        vals = [summ.get("high_confidence", 0), summ.get("medium_confidence", 0),
                summ.get("low_confidence", 0)]
        labels = ["HIGH", "MEDIUM", "LOW"]
        w, _, at = ax.pie(vals, labels=None, autopct=lambda p: f"{int(round(p*sum(vals)/100))}",
                          colors=[CONF_COLOR["HIGH"], CONF_COLOR["MEDIUM"],
                                  CONF_COLOR["LOW"]],
                          wedgeprops=dict(edgecolor="white", lw=1.5),
                          textprops=dict(fontsize=8, color="white",
                                         fontweight="bold"))
        ax.legend(w, [f"{l} ({v})" for l, v in zip(labels, vals)], fontsize=8,
                  frameon=False, loc="center left", bbox_to_anchor=(0.95, 0.5))
        ax.set_title(f"(a) confidence of {sum(vals)} candidates")
        ax.grid(False)

    ax = axes[1]
    if summ:
        pc = summ.get("per_caller", {})
        if pc:
            k = list(pc.keys()); v = [pc[x] for x in k]
            ax.barh(k[::-1], v[::-1], color=C_REF)
            for i, vv in enumerate(v[::-1]):
                ax.text(vv, i, f" {vv}", va="center", fontsize=8)
            ax.set_xlabel("variants called")
    ax.set_title("(b) calls per method")
    ax.set_ylabel("")

    ax = axes[2]
    if conc:
        names = [r.get("caller", "") for r in conc]
        M = np.array([[float(r.get(n, 0) or 0) for n in names] for r in conc])
        im = ax.imshow(M, cmap="YlGnBu")
        ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=45,
                                                             ha="right", fontsize=7)
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=7)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                ax.text(j, i, int(M[i, j]), ha="center", va="center", fontsize=6.5,
                        color="white" if M[i, j] > M.max() * 0.6 else "black")
        ax.set_title("(c) caller concordance")
        ax.grid(False)
        fig.colorbar(im, ax=ax, fraction=0.046, label="shared variants")
    fig.suptitle("Figure 8  Variant-caller concordance and confidence",
                 y=1.02, fontsize=12, fontweight="bold")
    fig.tight_layout()
    save(fig, outdir, "fig08_variant_confidence")


# ================================================== fig 09 substitution spectrum
def fig_substitution_spectrum(ctx, outdir):
    rows = read_tsv(os.path.join(ctx["out"], "03_variants/stats",
                                 "snp_substitution_spectrum.tsv"))
    if not rows:
        return
    subs, cls = [], []
    mode = None
    for r in rows:
        if "substitution" in r:
            mode = "sub"; continue
        if "class" in r:
            mode = "cls"; continue
        if mode == "sub":
            subs.append((r.get("substitution", ""), float(r.get("count", 0) or 0)))
        elif mode == "cls":
            cls.append((r.get("class", ""), float(r.get("count", 0) or 0)))
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
    ax = axes[0]
    if subs:
        labels = [s for s, _ in subs]; vals = [v for _, v in subs]
        colors = ["#1b7837" if s[0] in "AG" and s[2] in "AG" or
                  s[0] in "CT" and s[2] in "CT" else "#b2182b" for s in labels]
        ax.bar(labels, vals, color=colors)
        ax.set_ylabel("count"); ax.set_title("(a) SNP substitution spectrum")
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7.5)
    ax = axes[1]
    if cls:
        ax.barh([c for c, _ in cls][::-1], [v for _, v in cls][::-1],
                color=C_REF)
        for i, (c, v) in enumerate(cls[::-1]):
            ax.text(v, i, f" {int(v)}", va="center", fontsize=7.5)
        ax.set_xlabel("count"); ax.set_title("(b) substitution classes")
    fig.suptitle("Figure 9  Nucleotide substitution spectrum of the variant "
                 "candidates", y=1.05, fontsize=12, fontweight="bold")
    fig.tight_layout()
    save(fig, outdir, "fig09_substitution_spectrum")


# =============================================== fig 10 variant verification
def fig_variant_verification(ctx, outdir):
    rows = read_tsv(os.path.join(ctx["out"], "03_variants/verified",
                                 "verification_windows.tsv"))
    if not rows:
        return
    rows = rows[:6]
    fig, axes = plt.subplots(len(rows), 1, figsize=(13, 2.0 * len(rows)))
    if len(rows) == 1:
        axes = [axes]
    for ax, r in zip(axes, rows):
        ref = r.get("reference_sequence", "")
        qry = r.get("query_sequence", "")
        ws = int(r.get("window_start", 1))
        try:
            pos = int(r.get("pos"))
        except ValueError:
            pos = ws
        n = max(len(ref), len(qry))
        # draw bases as text rows, highlighting differences
        for i, (b1, b2) in enumerate(zip(ref.ljust(n), qry.ljust(n))):
            same = (b1 == b2)
            x = i
            ax.add_patch(Rectangle((x, 1.05), 0.94, 0.85,
                                   color="#f0f0f0" if same else "#ffe0e0",
                                   ec="white", lw=0.4))
            ax.text(x + 0.47, 1.47, b1, ha="center", va="center", fontsize=6.5,
                    family="monospace", color="#333" if same else C_REF)
            ax.add_patch(Rectangle((x, 0.05), 0.94, 0.85,
                                   color="#f0f0f0" if same else "#e0f0e0",
                                   ec="white", lw=0.4))
            ax.text(x + 0.47, 0.47, b2, ha="center", va="center", fontsize=6.5,
                    family="monospace", color="#333" if same else C_OK,
                    fontweight="normal" if same else "bold")
        ax.set_xlim(-0.5, n + 0.5); ax.set_ylim(-0.15, 2.1)
        ax.set_yticks([]); ax.set_xticks([])
        ax.grid(False)
        ax.text(-1.2, 1.47, "REF", ha="right", va="center", fontsize=7.5,
                fontweight="bold", color=C_REF)
        ax.text(-1.2, 0.47, "QRY", ha="right", va="center", fontsize=7.5,
                fontweight="bold", color=C_OK)
        ax.set_title(f"{r.get('contig')}:{pos}  {r.get('type')}  "
                     f"({r.get('supporting_callers','')})   "
                     f"window {ws}–{ws + n - 1}   "
                     f"reads: {r.get('observed_read_depth','?')} "
                     f"alt {r.get('observed_alt_percent','') or 'n/a'}%",
                     fontsize=8, loc="left")
    fig.suptitle("Figure 10  Sequence-level verification of every accepted "
                 "variant (reference above, query below)",
                 y=1.0, fontsize=12, fontweight="bold")
    fig.tight_layout()
    save(fig, outdir, "fig10_variant_verification")


# ===================================================== fig 11 gene impact
def fig_gene_impact(ctx, outdir):
    rows = read_tsv(os.path.join(ctx["out"], "04_annotation/variant_effect",
                                 "gene_impact.tsv"))
    if not rows:
        return
    fig, axes = plt.subplots(len(rows), 1, figsize=(13, 1.85 * len(rows)))
    if len(rows) == 1:
        axes = [axes]
    for ax, r in zip(axes, rows):
        try:
            start, end = int(r["start"]), int(r["end"])
        except ValueError:
            continue
        L = end - start + 1
        status = r.get("reference_gene_status", "")
        n_v = int(r.get("n_variants", 1) or 1)
        col = C_QRY if status == "pseudogene" else C_REF
        ax.add_patch(Rectangle((start, 0.5), L, 0.22, color=col, alpha=0.35,
                               ec=col, lw=1))
        # variant marker
        vpos = r.get("variant_positions", "")
        for chunk in vpos.split(";"):
            if ":" not in chunk:
                continue
            try:
                p = int(chunk.split(":")[1].split(">")[0].rstrip("ACGTN"))
            except (ValueError, IndexError):
                continue
            ax.plot([p, p], [0.42, 0.80], color=C_OK, lw=2.5, solid_capstyle="butt")
        # ORF integrity bars
        rl = float(r.get("ref_protein_len", 0) or 0)
        ql = float(r.get("qry_protein_len", 0) or 0)
        rs = float(r.get("ref_internal_stops", 0) or 0)
        qs = float(r.get("qry_internal_stops", 0) or 0)
        ax.barh(1.02, max(rl, 1), height=0.14, left=start,
                color="#cccccc" if rs else C_REF, alpha=0.8)
        ax.barh(0.20, max(ql, 1) * (L / max(rl, 1)), height=0.14, left=start,
                color="#bdbdbd" if qs else C_OK, alpha=0.8)
        ax.text(start, 1.24, f"reference protein: {rl:.0f} aa, "
                f"{rs:.0f} internal stop(s)", fontsize=7, color="#444")
        ax.text(start, 0.045, f"query protein: {ql:.0f} aa, "
                f"{qs:.0f} internal stop(s) — {r.get('query_orf_status','')}",
                fontsize=7, color="#444")
        ax.set_ylim(-0.05, 1.45)
        ax.set_xlim(start - L * 0.05, end + L * 0.15)
        ax.set_yticks([]); ax.grid(False)
        ax.xaxis.set_major_formatter(lambda v, pd: f"{v/1000:.1f} kb")
        title = (f"{r.get('locus_tag')}  ·  {r.get('product')}  ·  "
                 f"{r.get('worst_effect')} ({r.get('worst_impact')})")
        if status == "pseudogene":
            title += "  ·  reference annotated as PSEUDOGENE"
        ax.set_title(title, fontsize=8.5, loc="left")
    fig.suptitle("Figure 11  Consequence of each verified variant for the "
                 "affected gene", y=1.0, fontsize=12, fontweight="bold")
    fig.tight_layout()
    save(fig, outdir, "fig11_gene_impact")


# ================================================ fig 12 functional profile
def fig_functional_profile(ctx, outdir):
    rows = read_tsv(os.path.join(ctx["out"], "05_functional/cog_kegg",
                                 "cog_category_comparison.tsv"))
    if not rows:
        return
    letters = [r["COG_letter"] for r in rows]
    names = [r["COG_category"] for r in rows]
    rv = [float(r["reference_count"]) for r in rows]
    qv = [float(r["query_count"]) for r in rows]
    meta = [r.get("metabolism", "no") == "yes" for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2),
                             gridspec_kw={"width_ratios": [2.1, 1]})
    ax = axes[0]
    order = np.argsort(rv)[::-1]
    y = np.arange(len(letters))
    ax.barh(y - 0.2, [rv[i] for i in order], 0.4, color=C_REF,
            label=f"reference {ctx['ref_name']}")
    ax.barh(y + 0.2, [qv[i] for i in order], 0.4, color=C_QRY,
            label=f"query {ctx['qry_name']}")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{letters[i]}  {names[i][:44]}" for i in order],
                       fontsize=7)
    for k, i in enumerate(order):
        if meta[i]:
            ax.text(-0.01, k, "●", transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=6, color="#1b7837")
    ax.invert_yaxis()
    ax.set_xlabel("proteins assigned to the category")
    ax.set_title("(a) COG functional categories (● = metabolism)")
    ax.legend(fontsize=8, frameon=False)

    ax = axes[1]
    diff = [qv[i] - rv[i] for i in order]
    colors = [C_OK if d > 0 else C_QRY if d < 0 else "#bbbbbb" for d in diff]
    ax.barh(y, diff, color=colors)
    ax.set_yticks(y); ax.set_yticklabels([letters[i] for i in order], fontsize=7)
    ax.invert_yaxis(); ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("query − reference (proteins)")
    ax.set_title("(b) difference per category")
    fig.suptitle("Figure 12  Functional (COG) profile of the two genomes",
                 y=1.0, fontsize=12, fontweight="bold")
    fig.tight_layout()
    save(fig, outdir, "fig12_functional_profile")


# ================================================= fig 13 AMR determinants
def fig_amr_determinants(ctx, outdir):
    cat = read_tsv(os.path.join(ctx["out"], "05_functional/curated_scan",
                                "determinant_presence_matrix.tsv"))
    amrf = read_tsv(os.path.join(ctx["out"], "05_functional/amrfinder",
                                 "query.amrfinder.tsv"))
    curated = read_tsv(os.path.join(ctx["out"], "05_functional/curated_scan",
                                    "determinant_comparison.tsv"))
    if not cat and not curated:
        return
    fig, axes = plt.subplots(1, 2, figsize=(14, 6),
                             gridspec_kw={"width_ratios": [1.25, 1]})

    # (a) curated determinants: hits per category
    ax = axes[0]
    bycat = defaultdict(Counter)
    for r in curated:
        bycat[r.get("category", "?")][r.get("determinant", "")] += 1
    cats = ["resistance", "virulence", "metabolism"]
    cmap = {"resistance": C_QRY, "virulence": C_WARN, "metabolism": "#1b7837"}
    ypos, ylabels = [], []
    y = 0
    for c in cats:
        items = sorted(bycat.get(c, {}).items(), key=lambda kv: -kv[1])
        for det, n in items:
            ax.barh(y, n, color=cmap[c], alpha=0.85)
            ax.text(n, y, f" {n}", va="center", fontsize=6.5)
            ylabels.append(det); ypos.append(y); y += 1
        y += 1.2
    ax.set_yticks(ypos); ax.set_yticklabels(ylabels, fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("reference genes matching the determinant")
    ax.set_title("(a) curated B. fragilis determinant scan")
    ax.legend(handles=[Line2D([], [], color=cmap[c], lw=6, label=c) for c in cats],
              fontsize=7.5, frameon=False, loc="lower right")

    # (b) AMRFinderPlus hits
    ax = axes[1]
    if amrf and len(amrf) > 1:
        hdr = list(amrf[0].keys())
        sym = "Element symbol" if "Element symbol" in hdr else hdr[5]
        cls = "Class" if "Class" in hdr else hdr[8]
        rows = amrf[1:]
        labels = [f"{r.get(sym,'')}\n({r.get(cls,'')})" for r in rows]
        vals = [float(r.get("% Identity to reference sequence", 100) or 100)
                for r in rows]
        ax.barh(labels, vals, color=C_QRY)
        for i, v in enumerate(vals):
            ax.text(v, i, f" {v:.1f}%", va="center", fontsize=7.5)
        ax.set_xlim(0, 118)
        ax.set_xlabel("% identity to the AMRFinderPlus reference sequence")
    else:
        ax.text(0.5, 0.5, "no AMRFinderPlus hits", ha="center", va="center",
                transform=ax.transAxes)
    ax.set_title("(b) NCBI AMRFinderPlus — acquired resistance genes")
    fig.suptitle("Figure 13  Resistance / virulence / metabolism determinants",
                 y=1.0, fontsize=12, fontweight="bold")
    fig.tight_layout()
    save(fig, outdir, "fig13_amr_determinants")


# ====================================================================== main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--ref-fa", required=True)
    ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--ref-name", default="NCTC9343")
    ap.add_argument("--qry-name", default="CCUG4856")
    ap.add_argument("--ref-chr", default="NC_003228.3")
    a = ap.parse_args()
    os.makedirs(a.figdir, exist_ok=True)
    ctx = dict(out=a.out, ref_fa=a.ref_fa, ref_gff=a.ref_gff,
               ref_name=a.ref_name, qry_name=a.qry_name, ref_chr=a.ref_chr)

    figures = [
        ("assembly QC", fig_assembly_qc),
        ("read QC", fig_read_qc),
        ("synteny dot plot", fig_synteny_dotplot),
        ("structural variants", fig_structural_variants),
        ("genome circos", fig_circos),
        ("coverage", fig_coverage),
        ("variant landscape", fig_variant_landscape),
        ("variant confidence", fig_variant_confidence),
        ("substitution spectrum", fig_substitution_spectrum),
        ("variant verification", fig_variant_verification),
        ("gene impact", fig_gene_impact),
        ("functional profile", fig_functional_profile),
        ("AMR determinants", fig_amr_determinants),
    ]
    ok = 0
    for name, fn in figures:
        try:
            fn(ctx, a.figdir)
            ok += 1
        except Exception as exc:
            print(f"  [fig] {name} FAILED: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
    print(f"[figures] {ok}/{len(figures)} figures written to {a.figdir}")


if __name__ == "__main__":
    main()
