
[1;34m==> 01-QC (a): assembly statistics[0m
[assembly_stats] reference_NCTC9343: 2 contigs, 5,241,700 bp, N50=5,205,140, GC=43.11%
[assembly_stats] query_CCUG4856: 2 contigs, 5,241,693 bp, N50=5,205,133, GC=43.11%
[qc] wrote /home/cz/WGS/result/01_qc/assembly_qc/assembly_stats_comparison.tsv

[1;34m==> 01-QC (b): QUAST assembly reports[0m
WARNING: Python locale settings can't be changed
/home/cz/miniconda3/envs/WGS/bin/quast.py -o /home/cz/WGS/result/01_qc/quast/reference -t 32 --no-icarus --silent /home/cz/WGS/reference/unzipped/reference.fna


System information:
  OS: Linux-7.0.0-30-generic-x86_64-with-glibc2.39 (linux_64)
  Python version: 3.11.17
  CPUs number: 192

Started: 2026-10-12 17:27:38

Logging to /home/cz/WGS/result/01_qc/quast/reference/quast.log
NOTICE: Output directory already exists and looks like a QUAST output dir. Existing results can be reused (e.g. previously generated alignments)!

CWD: /home/cz/WGS
Main parameters: 
  MODE: default, threads: 32, min contig length: 500, min alignment length: 65, min alignment IDY: 95.0, \
  ambiguity: one, min local misassembly length: 200, min extensive misassembly length: 1000

Contigs:
  Pre-processing...
  /home/cz/WGS/reference/unzipped/reference.fna ==> reference

2026-10-12 17:27:38
Running Basic statistics processor...
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done.

NOTICE: Genes are not predicted by default. Use --gene-finding or --glimmer option to enable it.

2026-10-12 17:27:39
Creating large visual summaries...
This may take a while: press Ctrl-C to skip this step..
  1 of 1: Creating PDF with all tables and plots...
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done

2026-10-12 17:27:39
RESULTS:
  Text versions of total report are saved to /home/cz/WGS/result/01_qc/quast/reference/report.txt, report.tsv, and report.tex
  Text versions of transposed total report are saved to /home/cz/WGS/result/01_qc/quast/reference/transposed_report.txt, transposed_report.tsv, and transposed_report.tex
  PDF version (tables and plots) is saved to /home/cz/WGS/result/01_qc/quast/reference/report.pdf
  Log is saved to /home/cz/WGS/result/01_qc/quast/reference/quast.log

Finished: 2026-10-12 17:27:39
Elapsed time: 0:00:01.489216
NOTICEs: 2; WARNINGs: 1; non-fatal ERRORs: 0

Thank you for using QUAST!
   [ok]   quast_reference
WARNING: Python locale settings can't be changed
/home/cz/miniconda3/envs/WGS/bin/quast.py -o /home/cz/WGS/result/01_qc/quast/query -t 32 --no-icarus --silent /home/cz/WGS/seq/unzipped/query.fna


System information:
  OS: Linux-7.0.0-30-generic-x86_64-with-glibc2.39 (linux_64)
  Python version: 3.11.17
  CPUs number: 192

Started: 2026-10-12 17:27:39

Logging to /home/cz/WGS/result/01_qc/quast/query/quast.log
NOTICE: Output directory already exists and looks like a QUAST output dir. Existing results can be reused (e.g. previously generated alignments)!

CWD: /home/cz/WGS
Main parameters: 
  MODE: default, threads: 32, min contig length: 500, min alignment length: 65, min alignment IDY: 95.0, \
  ambiguity: one, min local misassembly length: 200, min extensive misassembly length: 1000

Contigs:
  Pre-processing...
  /home/cz/WGS/seq/unzipped/query.fna ==> query

2026-10-12 17:27:40
Running Basic statistics processor...
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done.

NOTICE: Genes are not predicted by default. Use --gene-finding or --glimmer option to enable it.

2026-10-12 17:27:41
Creating large visual summaries...
This may take a while: press Ctrl-C to skip this step..
  1 of 1: Creating PDF with all tables and plots...
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done

2026-10-12 17:27:41
RESULTS:
  Text versions of total report are saved to /home/cz/WGS/result/01_qc/quast/query/report.txt, report.tsv, and report.tex
  Text versions of transposed total report are saved to /home/cz/WGS/result/01_qc/quast/query/transposed_report.txt, transposed_report.tsv, and transposed_report.tex
  PDF version (tables and plots) is saved to /home/cz/WGS/result/01_qc/quast/query/report.pdf
  Log is saved to /home/cz/WGS/result/01_qc/quast/query/quast.log

Finished: 2026-10-12 17:27:41
Elapsed time: 0:00:01.471759
NOTICEs: 2; WARNINGs: 1; non-fatal ERRORs: 0

Thank you for using QUAST!
   [ok]   quast_query
WARNING: Python locale settings can't be changed
/home/cz/miniconda3/envs/WGS/bin/quast.py -o /home/cz/WGS/result/01_qc/quast/comparative -t 32 -r /home/cz/WGS/reference/unzipped/reference.fna --no-icarus --silent /home/cz/WGS/seq/unzipped/query.fna


System information:
  OS: Linux-7.0.0-30-generic-x86_64-with-glibc2.39 (linux_64)
  Python version: 3.11.17
  CPUs number: 192

Started: 2026-10-12 17:27:41

Logging to /home/cz/WGS/result/01_qc/quast/comparative/quast.log
NOTICE: Output directory already exists and looks like a QUAST output dir. Existing results can be reused (e.g. previously generated alignments)!

CWD: /home/cz/WGS
Main parameters: 
  MODE: default, threads: 32, min contig length: 500, min alignment length: 65, min alignment IDY: 95.0, \
  ambiguity: one, min local misassembly length: 200, min extensive misassembly length: 1000

Reference:
  /home/cz/WGS/reference/unzipped/reference.fna ==> reference

Contigs:
  Pre-processing...
  /home/cz/WGS/seq/unzipped/query.fna ==> query

2026-10-12 17:27:42
Running Basic statistics processor...
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done.

2026-10-12 17:27:43
Running Contig analyzer...
WARNING: Version of installed minimap2 differs from its version in the QUAST package (2.28). Please make sure that you use an actual version of software.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done.

2026-10-12 17:27:44
Running NA-NGA calculation...
Done.

2026-10-12 17:27:44
Running Genome analyzer...
  NOTICE: No file with genomic features were provided. Use the --features option if you want to specify it.

  NOTICE: No file with operons were provided. Use the -O option if you want to specify it.
Done.

NOTICE: Genes are not predicted by default. Use --gene-finding or --glimmer option to enable it.

2026-10-12 17:27:45
Creating large visual summaries...
This may take a while: press Ctrl-C to skip this step..
  1 of 1: Creating PDF with all tables and plots...
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done

2026-10-12 17:27:46
RESULTS:
  Text versions of total report are saved to /home/cz/WGS/result/01_qc/quast/comparative/report.txt, report.tsv, and report.tex
  Text versions of transposed total report are saved to /home/cz/WGS/result/01_qc/quast/comparative/transposed_report.txt, transposed_report.tsv, and transposed_report.tex
  PDF version (tables and plots) is saved to /home/cz/WGS/result/01_qc/quast/comparative/report.pdf
  Log is saved to /home/cz/WGS/result/01_qc/quast/comparative/quast.log

Finished: 2026-10-12 17:27:46
Elapsed time: 0:00:04.767222
NOTICEs: 4; WARNINGs: 2; non-fatal ERRORs: 0

Thank you for using QUAST!
   [ok]   quast_comparative
WARNING: Python locale settings can't be changed
/home/cz/miniconda3/envs/WGS/bin/quast.py -o /home/cz/WGS/result/01_qc/quast/comparative_genes -t 32 -r /home/cz/WGS/reference/unzipped/reference.fna -g /home/cz/WGS/reference/unzipped/reference.gff --no-icarus --silent /home/cz/WGS/seq/unzipped/query.fna


System information:
  OS: Linux-7.0.0-30-generic-x86_64-with-glibc2.39 (linux_64)
  Python version: 3.11.17
  CPUs number: 192

Started: 2026-10-12 17:27:46

Logging to /home/cz/WGS/result/01_qc/quast/comparative_genes/quast.log
NOTICE: Output directory already exists and looks like a QUAST output dir. Existing results can be reused (e.g. previously generated alignments)!

CWD: /home/cz/WGS
Main parameters: 
  MODE: default, threads: 32, min contig length: 500, min alignment length: 65, min alignment IDY: 95.0, \
  ambiguity: one, min local misassembly length: 200, min extensive misassembly length: 1000

Reference:
  /home/cz/WGS/reference/unzipped/reference.fna ==> reference

Contigs:
  Pre-processing...
  /home/cz/WGS/seq/unzipped/query.fna ==> query

2026-10-12 17:27:47
Running Basic statistics processor...
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done.

2026-10-12 17:27:48
Running Contig analyzer...
WARNING: Version of installed minimap2 differs from its version in the QUAST package (2.28). Please make sure that you use an actual version of software.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done.

2026-10-12 17:27:49
Running NA-NGA calculation...
Done.

2026-10-12 17:27:49
Running Genome analyzer...
  NOTICE: No file with operons were provided. Use the -O option if you want to specify it.
Done.

NOTICE: Genes are not predicted by default. Use --gene-finding or --glimmer option to enable it.

2026-10-12 17:27:50
Creating large visual summaries...
This may take a while: press Ctrl-C to skip this step..
  1 of 1: Creating PDF with all tables and plots...
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
findfont: Failed to find font weight medium for DejaVu Sans, now using 400.
Done

2026-10-12 17:27:51
RESULTS:
  Text versions of total report are saved to /home/cz/WGS/result/01_qc/quast/comparative_genes/report.txt, report.tsv, and report.tex
  Text versions of transposed total report are saved to /home/cz/WGS/result/01_qc/quast/comparative_genes/transposed_report.txt, transposed_report.tsv, and transposed_report.tex
  PDF version (tables and plots) is saved to /home/cz/WGS/result/01_qc/quast/comparative_genes/report.pdf
  Log is saved to /home/cz/WGS/result/01_qc/quast/comparative_genes/quast.log

Finished: 2026-10-12 17:27:51
Elapsed time: 0:00:05.123942
NOTICEs: 3; WARNINGs: 2; non-fatal ERRORs: 0

Thank you for using QUAST!
   [ok]   quast_comparative_genes

[1;34m==> 01-QC (c): in-silico reads from the query assembly + read-level QC[0m
processed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
[1A[Jprocessed files:  0 / 2 [------] ETA: 0s
file	format	type	num_seqs	sum_len	min_len	avg_len	max_len	Q1	Q2	Q3	sum_gap	N50	N50_num	Q20(%)	Q30(%)	AvgQual	GC(%)	sum_n
/home/cz/WGS/work/reads/CCUG4856_SIM_R1.fastq	FASTQ	DNA	1747231	262084650	150	150.0	150	150	150	150	0	150	873616	100	0	23.00	43.15	0
[1A[Jprocessed files:  1 / 2 [==>---] ETA: 0s
/home/cz/WGS/work/reads/CCUG4856_SIM_R2.fastq	FASTQ	DNA	1747231	262084650	150	150.0	150	150	150	150	0	150	873616	100	0	23.00	43.14	0
[1A[Jprocessed files:  2 / 2 [======] ETA: 0s
[1A[Jprocessed files:  2 / 2 [] ETA: 0s. done
[1A[Jprocessed files:  2 / 2 [] ETA: 0s. done
   [ok]   seqkit_stats_reads
[WARNING] You can't use more than 4 threads per file - setting this many is a waste
   [ok]   fastqc_reads
Read1 before filtering:
total reads: 1747231
total bases: 262084650
Q20 bases: 262084650(100%)
Q30 bases: 0(0%)
Q40 bases: 0(0%)

Read2 before filtering:
total reads: 1747231
total bases: 262084650
Q20 bases: 262084650(100%)
Q30 bases: 0(0%)
Q40 bases: 0(0%)

Read1 after filtering:
total reads: 1747201
total bases: 262075930
Q20 bases: 262075930(100%)
Q30 bases: 0(0%)
Q40 bases: 0(0%)

Read2 after filtering:
total reads: 1747201
total bases: 262075930
Q20 bases: 262075930(100%)
Q30 bases: 0(0%)
Q40 bases: 0(0%)

Filtering result:
reads passed filter: 3494402
reads failed due to low quality: 0
reads failed due to too many N: 0
reads failed due to too short: 60
reads failed due to adapter dimer: 0
reads with adapter trimmed: 274
bases trimmed due to adapters: 15208

Duplication rate: 0.00228934%

Insert size peak (evaluated by paired-end reads): 266

JSON report: /home/cz/WGS/result/01_qc/reads_qc/fastp.json
HTML report: /home/cz/WGS/result/01_qc/reads_qc/fastp.html

fastp -i /home/cz/WGS/work/reads/CCUG4856_SIM_R1.fastq -I /home/cz/WGS/work/reads/CCUG4856_SIM_R2.fastq -o /home/cz/WGS/work/reads/CCUG4856_SIM_R1.clean.fq.gz -O /home/cz/WGS/work/reads/CCUG4856_SIM_R2.clean.fq.gz -j /home/cz/WGS/result/01_qc/reads_qc/fastp.json -h /home/cz/WGS/result/01_qc/reads_qc/fastp.html -w 32 -q 20 -u 30 -l 50 
fastp v1.4.0, time used: 3 seconds
   [ok]   fastp_reads

[1;34m==> 01-QC (d): reference annotation QC[0m
[qc] annotation features: {'gene': 4254, 'CDS': 4231, 'exon': 94, 'tRNA': 73, 'pseudogene': 70, 'rRNA': 19}

[1;34m==> 01-QC complete -> /home/cz/WGS/result/01_qc[0m
