# 输入输出详细说明（GUI & CLI）

本文档详细说明 rnaseq-de-app 的**输入格式规范**、**界面控件含义**和**输出内容解读**，
适用于 Streamlit 网页界面与命令行两种方式。

---

## 1. 快速开始

### 网页界面（Streamlit）

```bash
pip install streamlit
streamlit run streamlit_app.py      # 浏览器打开 http://localhost:8501
```

侧边栏勾选「使用内置示例数据」→ 点「🚀 运行分析」，约 20 秒出结果。
云端部署版本（Streamlit Cloud）用法完全相同。

### 命令行（CLI）

```bash
python run_pipeline.py \
    --fastq-dir data/fastq \
    --samplesheet data/samplesheet.csv \
    --transcriptome data/ref/transcriptome.fa \
    --outdir results \
    --de-backend ttest \
    --cond-a control --cond-b case
```

---

## 2. 输入规范

### 2.1 样本分组表（CSV，必填）

四列，列名必须完全一致：

| 列名 | 含义 | 示例 |
|---|---|---|
| `sample` | 样本名（唯一，用于图表标注） | `case_1` |
| `fastq_1` | R1 文件名（与上传文件名一致） | `case_1_R1.fastq.gz` |
| `fastq_2` | R2 文件名 | `case_1_R2.fastq.gz` |
| `condition` | 分组名，**恰好 2 组** | `case` / `control` |

完整示例：

```csv
sample,fastq_1,fastq_2,condition
control_1,control_1_R1.fastq.gz,control_1_R2.fastq.gz,control
control_2,control_2_R1.fastq.gz,control_2_R2.fastq.gz,control
control_3,control_3_R1.fastq.gz,control_3_R2.fastq.gz,control
case_1,case_1_R1.fastq.gz,case_1_R2.fastq.gz,case
case_2,case_2_R1.fastq.gz,case_2_R2.fastq.gz,case
case_3,case_3_R1.fastq.gz,case_3_R2.fastq.gz,case
```

**规则与建议：**

- 每组 ≥ 3 个生物学重复（否则差异分析无统计意义）
- 上传的 fastq 文件名必须与分组表中写的**一字不差**（页面会自动校验并列出缺失文件）
- 界面提供「下载分组表模板」按钮
- 只支持**双端（paired-end）**测序数据；单端支持在 Roadmap 中

### 2.2 fastq 文件

- 格式：FASTQ，可为 `.fq` / `.fastq` / `.fq.gz` / `.fastq.gz`（推荐 gzip 压缩，上传更快）
- 网页端多选上传；建议先跑示例数据熟悉流程再传真实数据
- 定量后端为 k-mer 伪比对（21-mer），支持 read 含少量测序错误

### 2.3 参考转录本 fasta（必填）

- 每条序列为一个转录本，**header 必须含基因名**，格式：

  ```
  &gt;tx_00001 gene=gene_0001
  ACGT...
  ```

- 同一基因的多个转录本（isoform）写同一个 `gene=` 名，定量时按分数分配并聚合到基因水平
- 无 `gene=` 标注的序列按"单转录本基因"处理（基因名 = 转录本名）
- 示例数据：`data/ref/transcriptome.fa`

### 2.4 内置示例数据

勾选后使用仓库自带的 6 样本数据（control×3 vs case×3），特点：

- 120 基因 / 230 转录本，isoform 共享主体序列（模拟真实可变剪接）
- 30 个已知差异基因（2–8 倍），真值记录在 `data/ground_truth.json`
- 用于教学演示和验证部署是否正确，**预期结果**：DE 基因 ≥28 个、比对率 98.8%+

---

## 3. 参数说明（界面控件）

| 控件 | 作用 | 默认 |
|---|---|---|
| 使用内置示例数据 | 免上传直接试跑 | 勾选 |
| fastq 上传 / fasta 上传 / 分组表上传 | 见第 2 节 | — |
| 基准组 A（对照） | 参照组；log2FC 以此为基准 | 字母序第一组 |
| 处理组 B（case） | 实验组 | 另一组 |
| 差异分析后端 | `auto`：装了 pydeseq2 自动用 DESeq2（推荐，论文级方法），否则回退 t-test；可手动强制 | auto |
| 显著性阈值 padj | 差异基因判定线（FDR 校正后 p 值） | 0.10 |
| 最小 \|log2FC\| | 火山图红点标注范围，不影响结果表 | 1.0 |
| 🚀 运行分析 | 启动流水线：定量 → 差异分析 → 出图 | — |

&gt; **方向约定：所有 log2FoldChange = B 组 相对 A 组。** 正值 = B 上调，负值 = B 下调。

---

## 4. 输出说明

### 4.1 运行日志

运行过程中实时显示：每个样本的比对率、实际使用的差异分析后端。
比对率应 ≥ 90%；显著偏低说明 fasta 参考或文件配对有问题。

### 4.2 顶部指标卡

| 指标 | 含义 |
|---|---|
| DE 基因（padj&lt;阈值） | 显著差异基因总数 / 检测基因总数 |
| 上调 | log2FC &gt; 0（B 组升高） |
| 下调 | log2FC &lt; 0（B 组降低） |

### 4.3 差异基因表（📊 标签页）

按 padj 从小到大排序，仅显示显著基因：

| 列名 | 含义 | 解读 |
|---|---|---|
| `gene` | 基因名（来自 fasta 的 `gene=`） | — |
| `baseMean` | 标准化后所有样本的平均表达量 | 表达越高效应估计越稳 |
| `log2FoldChange` | B vs A 的倍数变化（log2） | +1 = 2 倍上调；+2 = 4 倍上调；-1 = 2 倍下调 |
| `pvalue` | 原始 p 值 | 仅供排序参考 |
| `padj` | BH/FDR 校正后 p 值 | **显著性判断看这个**，&lt; 0.1 通常认为显著 |

### 4.4 图表（🖼️ 标签页）

**火山图（volcano）**
- 横轴 log2FC，纵轴 −log10(p)，每点一个基因
- 灰点 = 不显著；红点 = padj &lt; 阈值 且 |log2FC| &gt; 设定值
- 看右上角（强上调）和左上角（强下调）；理想形态应左右大致对称

**MA 图**
- 横轴 log10(平均表达量)，纵轴 log2FC
- 用于检查系统性偏差：点应围绕 y=0 水平线对称分布；
  若整体向上/向下倾斜，提示标准化不彻底或存在批次效应

**PCA 图**
- 样本按表达谱的无监督聚类
- **质量检查关键图**：同组样本应聚成一簇，两组分开；
  有样本远离自己组别 → 该样本数据可能有问题（污染、批次、搞混）

**Top30 热图**
- 差异最显著的 30 个基因的表达模式（行标准化，红高蓝低）
- 应清晰呈现"按列分成两组"的模式

**比对率统计（展开）**
- 每个样本的 read 对总数、比对成功数、比对率

### 4.5 下载（💾 标签页）

| 文件 | 内容 |
|---|---|
| `de_results.csv` | **全部基因**的完整差异分析结果（含不显著基因），可导入 Excel/R |
| `raw_counts.csv` | 基因 × 样本的原始计数矩阵 |
| `normalized_counts.csv` | 中位数比值标准化后的计数矩阵（适合画表达量热图） |

### 4.6 命令行模式输出目录结构

```
results/
├── qc_mapping.csv              # 比对率统计
├── counts/
│   ├── raw_counts.csv          # 原始计数矩阵
│   ├── normalized_counts.csv   # 标准化计数矩阵
│   └── size_factors.csv        # 各样本标准化因子
├── de/
│   └── de_results.csv          # 差异分析完整结果
├── figures/
│   ├── volcano.png             # 火山图
│   ├── MA_plot.png             # MA 图
│   ├── PCA.png                 # PCA
│   └── heatmap_top30.png       # Top30 热图
└── run_info.json               # 运行参数记录
```

---

## 5. 结果质量自检清单

拿到结果先过一遍这四项：

1. ✅ 比对率各样本 ≥ 90%
2. ✅ PCA 中同组样本聚在一起
3. ✅ MA 图围绕 y=0 对称、无整体倾斜
4. ✅ DE 基因数不是 0（真实数据 3v3、padj&lt;0.1 下通常几十个到上千个）

任一项异常 → 先检查输入数据和分组，再考虑调参。

---

## 6. 常见问题（FAQ）

**Q: 网页报"分组表引用但未上传的文件"？**
分组表 `fastq_1/fastq_2` 列写的名字和实际上传的文件名不一致（含大小写、后缀）。
R1/R2 都要传。

**Q: 比对率很低？**
参考 fasta 与物种/样本不匹配，或 fasta header 格式不对。先跑内置示例数据对比。

**Q: DE 基因数为 0？**
真实差异可能不存在；或每组 &lt; 3 重复、测序深度太浅。可适当放宽 padj 到 0.2 观察趋势。

**Q: pydeseq2 和 ttest 结果差多少？**
本工具两种后端在验证数据上方向完全一致；正式分析建议 `pip install pydeseq2`
用 DESeq2（对离散 counts 建模更严格）。

**Q: 能把图换成可交互的（鼠标悬停看基因名）吗？**
可以，Roadmap 中有计划（Plotly 版图表），需要随时提 issue。
