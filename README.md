# WGS — *Bacteroides fragilis* CCUG 4856 vs NCTC 9343

Whole-genome comparison of ***Bacteroides fragilis*** strain **CCUG 4856**
(`GCF_005706655.1`) against the type/reference strain **NCTC 9343**
(`GCF_000025985.1`), covering quality control, genome alignment, structural
variation, variant detection, functional prediction, variant-to-gene annotation
(resistance / virulence / metabolism) and visualisation.

Everything is driven by a single entry point, **`WGS.sh`**, and runs inside the
conda environment `WGS`.

---

## Headline result

The two genomes are **extraordinarily close**: one collinear alignment covers
100 % of both replicons at **99.99 %** average identity. The analysis accepts
**four single-site differences** plus **one large structural rearrangement**.

| reference position (`NC_003228.3`) | change | type | independent methods |
|---|---|---|---|
| 530,487 | T → TT (extra T in a T-run) | insertion | 2 |
| 669,715 | TTA → T (T + A lost from a T₆A₆ run) | deletion | 2 |
| **683,862** | **G → A** | **SNP** | **3** |
| 1,709,303 | G → GT (extra T) | insertion | 2 |

Plus an **88,340 bp inversion** at `NC_003228.3:2,005,600–2,093,939`
(confirmed by a perfect reverse-complement alignment with no forward copy, so it
is not an inverted-repeat artefact).

### Genes affected

| locus tag | product | consequence |
|---|---|---|
| `BF9343_RS02670` | RagB/SusD family nutrient-uptake outer membrane protein | **A141T missense** (BLOSUM62 = 0; ORF intact in both strains) |
| `BF9343_RS02095` | BamA/OMP85 family outer membrane protein | reference is a **pseudogene** (55 internal stops); **query reading frame restored**, 878 aa intact |
| `BF9343_RS06655` | FtsX-like permease family protein | reference is a **pseudogene** (15 internal stops); **query reading frame restored**, 416 aa intact |

### Two things that are *not* biological differences

* **Origin shift (coordinate convention).** The query chromosome is the
  reference chromosome **re-linearised at a different origin** — the query is the
  reference rotated by 936,212 bp (98.29 % of the reference). A naive
  position-by-position comparison would report this as ~100 % of the genome
  differing. The alignment-based analysis here is immune to it.
* **Plasmid re-linearisation.** Both plasmids are exactly 36,560 bp and
  **sequence-identical**, but were deposited from different origins and on
  opposite strands. Allowing for the circular rotation and both orientations
  gives 100 % identity and zero differences. The shared plasmid carries 47 genes
  including a complete conjugative-transfer apparatus.

### What was *not* accepted

215 further candidate positions were generated and **rejected**, each with a
recorded reason, in `result/03_variants/verified/variants_ambiguous.tsv`. They
sit in multi-copy repeats (e.g. the ~5.8 kb rRNA-operon-like units at
≈3.20 Mb, ≈3.85 Mb and ≈4.66 Mb) or at rearrangement breakpoints, where short
reads cannot decide which paralogue is which. They are candidates, not
differences.

---

## Repository layout

```
WGS.sh                     the whole pipeline (8 steps, resumable)
scripts/
  config.sh                shared paths, sample names, thread counts
  01_prepare.sh            decompress inputs, build the tree
  01_qc.sh                 assembly QC (stats, QUAST), in-silico read QC
  02_alignment.sh          NUCmer WGA, dnadiff, contig mapping, coverage, SV
  03_variants.sh           4-way variant calling, reconciliation, verification
  04_annotate.sh           snpEff / bcftools csq, allele-aware annotation,
                           annotation transfer, gene-content comparison
  05_functional.sh         AMRFinderPlus, curated determinant scan, COG, GO/EC
  06_extra.sh              plasmid, prophage, CRISPR, mobile elements, phylogeny
  07_visualize.sh          all figures
  08_report.sh             final Markdown + HTML report
  lib/                     the analysis and plotting modules (Python)
reference/                 reference genome (NCTC 9343), FASTA/GFF/protein
seq/                       query genome (CCUG 4856)
result/                    all results, 01_qc … 08_report
logs/                      per-step logs and software versions
```

`work/` is **not** committed — it holds regenerable intermediates (BAMs, the
simulated reads, the downloaded genome panel and the ~1.4 GB COG2020 database,
≈5.3 GB in total). `WGS.sh` rebuilds all of it.

### Where to look first

| I want… | open |
|---|---|
| the complete report | **`result/08_report/WGS_report.html`** |
| the same as Markdown | `result/08_report/WGS_report.md` |
| headline numbers | `result/08_report/key_findings.tsv`, `RUN_SUMMARY.txt` |
| the verified variants | `result/03_variants/verified/variants_confident.tsv` |
| the rejected candidates | `result/03_variants/verified/variants_ambiguous.tsv` |
| side-by-side allele evidence | `result/03_variants/verified/verification_windows.tsv` |
| affected genes | `result/04_annotation/variant_effect/gene_impact.tsv` |
| every gene differ | `result/04_annotation/transferred/projection_table.tsv` |
| structural differences | `result/02_alignment/synteny/structural_variants.tsv` |
| resistance / virulence / metabolism | `result/05_functional/` |
| all figures | `result/07_visualization/` |

---

## Reproducing the analysis

```bash
source /home/cz/miniconda3/etc/profile.d/conda.sh
conda activate WGS
cd /home/cz/WGS

bash WGS.sh                 # every step, skipping work that is already done
bash WGS.sh --force         # re-run from scratch
bash WGS.sh --steps 3,4     # only variant calling and annotation
bash WGS.sh --check         # verify the toolchain, run nothing
bash WGS.sh --list          # list the steps
```

Results land in `result/`, logs in `logs/`. Each step is idempotent and
regenerates only what is missing.

### Toolchain

Required: `samtools`, `bcftools`, `bwa`, `minimap2`, `nucmer`/`show-snps`
(MUMmer4), `freebayes`, `dnadiff`, `python3` with `pysam`/`numpy`/`matplotlib`,
`quast.py`, `seqkit`, `wgsim`, `snpEff`, `amrfinder`, `diamond`, `mash`.

Optional (the pipeline degrades gracefully and says so in the log): `prokka`,
`bakta`, `abricate`, `phispy`, `mafft`, `iqtree`, `roary`, `Rscript`.

---

## Method notes

**Four independent variant routes** are reconciled rather than trusted
individually: MUMmer `show-snps` on the 1-to-1 NUCmer delta (assembly),
`bcftools call` and `freebayes` on the minimap2 `asm5` contig BAM
(assembly), and `bcftools`/`freebayes` on simulated reads.

A call is only accepted when **at least two of the three genuinely independent
methods** agree — NUCmer alignment, contig alignment, read alignment. Note that
`bcftools_asm` and `freebayes_asm` share one alignment, so their agreement is
*not* independent evidence. A call is additionally rejected when the reference
position is covered by more than one genomic copy, when the NUCmer alignment
reports overlapping alternatives, when the reads map with low confidence, or
when **the reads contradict the call** (the alternate allele is essentially
absent).

**Reads.** The query is a finished assembly, so no wet-lab reads exist. Illumina
reads were simulated from the query with `wgsim` (~100× , 150 bp, 350 bp insert,
0.5 % error) to provide a real read-level QC layer and, more importantly, an
independent confirmation of the assembly-derived calls. Everything derived from
them is labelled `SIM`.

**Gap placement is ambiguous in homopolymer runs**, so the aligner's first
description of an indel is not necessarily the biologically accurate one. The
pipeline reconstructs both alleles from the alignment and re-reduces them to a
minimal REF/ALT pair; that is why the indels are reported as an extra T, a
T+A deletion, and an extra T rather than at the coordinates a caller happened to
choose.

**Annotation.** Because the genomes are collinear at 99.99 % identity, the
complete query annotation is produced by projecting every reference gene through
the alignment — query coordinates *and* query sequence — which yields an intact
query proteome of 4,231 proteins, 4,226 of them identical to the reference.

---

## Caveats

1. **Variant confirmation uses simulated reads**, which are exact copies of the
   query assembly. They cannot reveal errors already present in that assembly;
   orthogonal sequencing of the real strain remains the gold standard.
2. **215 candidate differences in repeats and at breakpoints are unresolved.**
   They are listed explicitly and must not be reported as differences; long-read
   or clone-level data would be needed.
3. **The origin shift is not a biological difference** and is reported separately
   so it is not mistaken for one.
4. **De-novo annotation with Bakta/Prokka could not be completed** in the build
   environment: Bakta's database is hosted on Zenodo, which was unreachable, and
   Prokka has an unsatisfiable Perl dependency in this conda environment.
   Functional annotation therefore rests on the reference RefSeq annotation
   transferred through the alignment, plus AMRFinderPlus and NCBI COG2020. This
   is appropriate at this level of identity but should be revisited if a Bakta
   database becomes available.
5. **ABRicate could not be installed** (same Perl conflict); AMRFinderPlus plus
   the curated 39-determinant scan cover the same ground with curated databases.
   AMRFinderPlus found *cepA* (CepA family extended-spectrum class A
   β-lactamase, cephalosporin resistance) in both strains.

---

## Software

Exact versions and input checksums: `result/08_report/software_versions.txt`.
Every step writes its own log under `logs/`.
