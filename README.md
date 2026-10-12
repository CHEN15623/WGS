# WGS —— 脆弱拟杆菌 CCUG 4856 与 NCTC 9343 全基因组比较

将 ***Bacteroides fragilis*（脆弱拟杆菌）** 菌株 **CCUG 4856**
（`GCF_005706655.1`）与模式/参考菌株 **NCTC 9343**（`GCF_000025985.1`）进行全基因组比较，
内容涵盖质量控制、基因组对齐、结构变异、变异检测、功能预测、变异位点基因注释
（耐药 / 毒力 / 代谢）以及结果可视化。

全流程由单一入口脚本 **`WGS.sh`** 驱动，在 conda 环境 `WGS` 中运行。

---

## 核心结论

两个基因组**高度接近**：一条共线性比对覆盖双方 100% 的序列，平均一致性达 **99.99%**。
分析最终确认 **4 个单碱基位点差异**，外加 **1 处大规模结构重排**。

| 参考基因组位置（`NC_003228.3`） | 变化 | 类型 | 独立方法数 |
|---|---|---|---|
| 530,487 | T → TT（T 串中多 1 个 T） | 插入 | 2 |
| 669,715 | TTA → T（T₆A₆ 串中缺失 T 和 A） | 缺失 | 2 |
| **683,862** | **G → A** | **SNP** | **3** |
| 1,709,303 | G → GT（多 1 个 T） | 插入 | 2 |

此外还有一处 **88,340 bp 倒位**，位于 `NC_003228.3:2,005,600–2,093,939`。
该区域存在完美的反向互补比对，且在正向无对应拷贝，因此可排除"反向重复序列造成的假象"。

### 受影响的基因

| locus tag | 产物 | 后果 |
|---|---|---|
| `BF9343_RS02670` | RagB/SusD 家族营养摄取外膜蛋白 | **A141T 错义突变**（BLOSUM62 = 0，保守替换；两菌株 ORF 均完整） |
| `BF9343_RS02095` | BamA/OMP85 家族外膜蛋白 | 参考菌株为**假基因**（55 个内部终止密码子）；**查询菌株阅读框恢复**，878 aa 完整 |
| `BF9343_RS06655` | FtsX 类通透酶家族蛋白 | 参考菌株为**假基因**（15 个内部终止密码子）；**查询菌株阅读框恢复**，416 aa 完整 |

### 两项**并非**生物学差异的现象

* **起点偏移（坐标约定问题）**：查询菌株染色体是参考染色体**从不同起点重新线性化**的结果
  —— 相当于把参考旋转 936,212 bp（占参考的 98.29%）。若采用朴素的逐位比较，会误报
  "全基因组约 100% 不同"。本分析基于比对，天然不受该问题影响。
* **质粒重新线性化**：两个质粒长度同为 36,560 bp 且**序列完全一致**，只是沉积（deposit）
  时起点不同、链方向相反。在允许环状旋转和两个方向后比较，一致性为 100%、差异为 0。
  该共享质粒携带 47 个基因，包含一整套完整的接合转移系统。

### 被**排除**的候选位点

流程共产生 215 个候选位点并予以**排除**，每一个都在
`result/03_variants/verified/variants_ambiguous.tsv` 中记录了具体原因。
它们位于多拷贝重复区（例如 ≈3.20 Mb、≈3.85 Mb、≈4.66 Mb 处约 5.8 kb 的类 rRNA 操纵子单元）
或重排断点处，短读长测序无法判定差异属于哪一个旁系同源拷贝。
这些是**候选**，不是已确认的差异。

---

## 仓库结构

```
WGS.sh                     完整流程脚本（8 个步骤，可断点续跑）
scripts/
  config.sh                共享路径、样本名、线程数配置
  01_prepare.sh            解压输入、构建目录树
  01_qc.sh                 组装质量评估（统计量、QUAST）、模拟读数质控
  02_alignment.sh          NUCmer 全基因组比对、dnadiff、contig 比对、覆盖度、结构变异
  03_variants.sh           四路变异检测、结果整合、验证
  04_annotate.sh           snpEff / bcftools csq、等位基因感知注释、
                           注释投射、基因内容比较
  05_functional.sh         AMRFinderPlus、精选决定因子扫描、COG、GO/EC
  06_extra.sh              质粒、前噬菌体、CRISPR、可移动元件、系统发育
  07_visualize.sh          全部图表
  08_report.sh             生成最终 Markdown + HTML 报告
  lib/                     分析与绘图模块（Python）
reference/                 参考基因组（NCTC 9343）FASTA/GFF/蛋白序列
seq/                       查询基因组（CCUG 4856）
result/                    全部结果，01_qc … 08_report
logs/                      各步骤日志与软件版本
```

`work/` **未纳入版本控制** —— 其中存放可再生的中间产物（BAM 文件、模拟读数、下载的基因组
面板，以及约 1.4 GB 的 COG2020 数据库，合计约 5.3 GB）。`WGS.sh` 会按需重建全部内容。

### 建议的阅读顺序

| 我想看…… | 打开 |
|---|---|
| 完整报告 | **`result/08_report/WGS_report.html`** |
| 同上的 Markdown 版本 | `result/08_report/WGS_report.md` |
| 关键数字速览 | `result/08_report/key_findings.tsv`、`RUN_SUMMARY.txt` |
| 经验证的变异位点 | `result/03_variants/verified/variants_confident.tsv` |
| 被排除的候选位点 | `result/03_variants/verified/variants_ambiguous.tsv` |
| 等位基因并排证据 | `result/03_variants/verified/verification_windows.tsv` |
| 受影响基因 | `result/04_annotation/variant_effect/gene_impact.tsv` |
| 全部有差异的基因 | `result/04_annotation/transferred/projection_table.tsv` |
| 结构差异 | `result/02_alignment/synteny/structural_variants.tsv` |
| 耐药 / 毒力 / 代谢 | `result/05_functional/` |
| 全部图表 | `result/07_visualization/` |

---

## 复现分析

```bash
source /home/cz/miniconda3/etc/profile.d/conda.sh
conda activate WGS
cd /home/cz/WGS

bash WGS.sh                 # 运行全部步骤，自动跳过已完成的工作
bash WGS.sh --force         # 从头重新运行
bash WGS.sh --steps 3,4     # 只运行变异检测与注释
bash WGS.sh --check         # 仅检查工具链，不执行任何分析
bash WGS.sh --list          # 列出所有步骤
```

结果输出到 `result/`，日志输出到 `logs/`。每个步骤都是幂等的，只重新生成缺失的部分。

### 工具链

必需：`samtools`、`bcftools`、`bwa`、`minimap2`、`nucmer`/`show-snps`（MUMmer4）、
`freebayes`、`dnadiff`、`python3`（含 `pysam`/`numpy`/`matplotlib`）、`quast.py`、
`seqkit`、`wgsim`、`snpEff`、`amrfinder`、`diamond`、`mash`。

可选（缺失时流程会自动降级，并在日志中明确说明）：`prokka`、`bakta`、`abricate`、
`phispy`、`mafft`、`iqtree`、`roary`、`Rscript`。

---

## 方法说明

**四路独立的变异检测路线**互相印证，而非单独采信任何一路：MUMmer `show-snps`
基于 1-to-1 NUCmer delta（组装层面）、`bcftools call` 与 `freebayes` 基于 minimap2
`asm5` contig 比对（组装层面）、以及 `bcftools`/`freebayes` 基于模拟读数（测序层面）。

一个位点只有在**三种真正独立的方法**（NUCmer 比对、contig 比对、读数比对）中
**至少两种支持**时才会被采纳。需要说明的是，`bcftools_asm` 与 `freebayes_asm`
共用同一份比对结果，因此它们的一致**不构成独立证据**。此外，出现下列任一情况时该位点
也会被排除：参考位置被多于一个基因组拷贝覆盖；NUCmer 比对报告了重叠的替代比对；
覆盖该位点的读数比对置信度低；或者**读数直接与该位点矛盾**（替代等位基因基本不存在）。

**关于读数**：查询菌株是完成图（finished assembly），本身不存在湿实验读数。
因此使用 `wgsim` 从查询组装序列模拟 Illumina 读数（约 100×、150 bp、350 bp 插入片段、
0.5% 错误率），一是提供真实的读数级质控层，二是更重要的——为组装层面得到的变异
提供独立确认。所有由模拟读数派生的一切均标注为 `SIM`。

**同聚物区（homopolymer run）中的空位放置本身是歧义的**，因此比对器最先给出的
indel 描述未必是生物学上准确的描述。流程会从比对中重建两个等位基因，并重新化简为
最小的 REF/ALT 表示；这也是为什么 indel 最终被报告为"多一个 T""缺失 T+A""多一个 T"，
而不是某个 caller 恰好选定的坐标。

**关于注释**：由于两个基因组在 99.99% 一致性下保持共线性，查询基因组的完整注释是通过
**把参考的每一个基因沿比对投射**得到 —— 同时获得查询坐标**和**查询序列 ——
由此得到 4,231 条完整的查询蛋白，其中 4,226 条与参考完全相同。

---

## 注意事项与局限

1. **变异的确认使用了模拟读数**，它们是查询组装序列的精确副本，因此无法揭示该组装
   序列本身已存在的错误。对真实菌株进行正交测序（Sanger、Illumina 或长读长）仍是金标准。
2. **重复区与断点处的 215 个候选差异无法解析**。它们在文件中被明确列出，不应当作
   已确认的差异报告；需要长读长或克隆水平的数据才能判定。
3. **起点偏移不属于生物学差异**，报告中单独列出，以免被误读为真实差异。
4. **Bakta/Prokka 的从头注释未能在本环境中完成**：Bakta 的数据库托管于 Zenodo，
   本机无法访问；Prokka 在当前 conda 环境中存在无法满足的 Perl 依赖冲突。
   因此功能注释建立在"参考 RefSeq 注释经比对投射"之上，并辅以 AMRFinderPlus 与
   NCBI COG2020。在这一同一性水平下该方案是恰当的，但若日后可获得 Bakta 数据库，
   建议重做这一部分。
5. **ABRicate 无法安装**（同样是 Perl 冲突）；AMRFinderPlus 加上精选的 39 项
   决定因子扫描已用带版本管理的策展数据库覆盖了同样的范围。AMRFinderPlus 在两个
   菌株中均检出 *cepA*（CepA 家族超广谱 A 类 β-内酰胺酶，头孢菌素耐药）。

---

## 软件版本

精确版本号与输入文件校验和见 `result/08_report/software_versions.txt`。
每个步骤各自在 `logs/` 下写入日志。
