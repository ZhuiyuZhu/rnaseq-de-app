# rnaseq-de-app

两组（pair-wise）RNA-seq 差异表达分析迷你流水线：**fastq 进，差异基因表 + 图表出**。
零重依赖（numpy / pandas / scipy / matplotlib 即可跑通），Windows 本地友好，自带可复现的合
成验证数据（模拟已知差异，自动校验流水线正确性）。

## 快速开始

```bash
pip install -r requirements.txt          # 可选: pip install pydeseq2 启用 DESeq2 后端
python scripts/make_test_data.py --out data --seed 42   # 生成自带已知差异的验证数据
python run_pipeline.py --fastq-dir data/fastq --samplesheet data/samplesheet.csv \
    --transcriptome data/ref/transcriptome.fa --outdir results
python scripts/validate_results.py       # 两层自检: 定量保真度 + DE 召回率
```

预期验证结果：**DE 基因召回率 ≥ 90%，误报率 ≈ 0%（FDR<0.1，seed=42 实测 93.3% / 0%）**。

## 输入

- `samplesheet.csv` 四列：`sample,fastq_1,fastq_2,condition`；两组设计，每组 ≥3 生物学重复
- 参考转录本 fasta，header 需含基因名：`>tx_00001 gene=gene_0001`
- fastq 支持 `.gz`；单端测序暂不支持（v0.2 计划）

## 架构（双后端可插拔）

```
fastq ──► [定量: python k-mer 伪比对(零依赖) | salmon] ──► 基因计数矩阵
      ──► [DE: PyDESeq2(推荐) | ttest 零依赖回退] ──► 火山图 / MA / PCA / 热图
```

- `--de-backend auto`：装了 pydeseq2 自动用 DESeq2，否则回退 t-test
- 所有图表与计数表落盘 `results/`，另附 `qc_mapping.csv` 比对率统计

## 自带验证数据的设计

- 120 基因 / 230 isoform 合成转录组；同基因 isoform **共享主体序列**（模拟真实可变剪接，
  考验定量的分数分配能力）
- 30 个已知 DE 基因（2–8x，选自高表达前 1/3），负二项抽样（disp=0.05），PE60，含 0.5% 测序错误
- `data/expected_counts.csv` = 每基因模拟真值计数（层1：与流水线 raw counts 相关性 = 1.0）
- `data/ground_truth.json` = DE 基因真实倍数（层2：召回率/误报率）

## 测试与 CI

```bash
pytest -q        # 端到端冒烟: 生成小数据→跑流水线→校验输出
```

push 到 GitHub 后 Actions 自动执行：生成数据 → 跑流水线 → 验证召回率。

## 本地 Windows 调试

```powershell
python -m venv .venv; .venv\Scripts\Activate.ps1
pip install -r requirements.txt
# (同"快速开始"三条命令)
```

## Roadmap

- [ ] Streamlit GUI（上传 fastq + 分组表 → 进度条 → 交互式火山图/热图，可部署 Streamlit Cloud）
- [ ] 富集分析（GSEA / clusterProfiler）
- [ ] HTML 报告（MultiQC 风格聚合）
- [ ] 单端测序支持

## License

MIT
