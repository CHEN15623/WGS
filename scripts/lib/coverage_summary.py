#!/usr/bin/env python3
"""
coverage_summary.py -- turn a `samtools depth -a` table into a coverage report.

Input : TSV with columns  contig, position, depth
Output: TSV with per-contig mean/median depth, breadth at >=1x/>=5x/>=10x,
        number and total length of zero-coverage (uncovered) gaps, longest gap.
        Uncovered intervals are also written next to the report.

Usage: coverage_summary.py <depth.tsv> <out.tsv>
"""
import sys
import os
import gzip
import statistics
import collections


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    depth_f, out_f = sys.argv[1], sys.argv[2]

    per = collections.OrderedDict()   # contig -> list of depths
    _opener = gzip.open if depth_f.endswith(".gz") else open
    with _opener(depth_f, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                pos, d = int(f[1]), int(f[2])
            except ValueError:
                continue
            per.setdefault(f[0], []).append(d)

    gaps_f = os.path.join(os.path.dirname(out_f), "uncovered_regions.tsv")
    with open(out_f, "w") as out, open(gaps_f, "w") as gaps:
        out.write("contig\tlength_bp\tmean_depth\tmedian_depth\tmax_depth\t"
                  "breadth_ge1x_percent\tbreadth_ge5x_percent\t"
                  "breadth_ge10x_percent\tuncovered_bp\tuncovered_percent\t"
                  "n_uncovered_regions\tlongest_uncovered_bp\n")
        gaps.write("contig\tstart\tend\tlength_bp\n")
        for name, depths in per.items():
            n = len(depths)
            if n == 0:
                continue
            mean = sum(depths) / n
            med = statistics.median(depths)
            mx = max(depths)
            b1 = 100.0 * sum(1 for d in depths if d >= 1) / n
            b5 = 100.0 * sum(1 for d in depths if d >= 5) / n
            b10 = 100.0 * sum(1 for d in depths if d >= 10) / n
            un = [i for i, d in enumerate(depths) if d == 0]
            # collapse into contiguous intervals
            intervals = []
            if un:
                start = prev = un[0]
                for i in un[1:]:
                    if i == prev + 1:
                        prev = i
                    else:
                        intervals.append((start, prev))
                        start = prev = i
                intervals.append((start, prev))
            for s, e in intervals:
                gaps.write(f"{name}\t{s+1}\t{e+1}\t{e-s+1}\n")
            longest = max((e - s + 1 for s, e in intervals), default=0)
            out.write(f"{name}\t{n}\t{mean:.2f}\t{med:.1f}\t{mx}\t"
                      f"{b1:.2f}\t{b5:.2f}\t{b10:.2f}\t{len(un)}\t"
                      f"{100.0*len(un)/n:.3f}\t{len(intervals)}\t{longest}\n")
    print(f"[coverage] wrote {out_f} and {gaps_f}")


if __name__ == "__main__":
    main()
