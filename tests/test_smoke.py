"""端到端冒烟测试: 生成小数据 -> 跑流水线 -> 校验输出。"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_pipeline_end2end(tmp_path):
    data = tmp_path / "data"
    res = tmp_path / "res"
    subprocess.run([sys.executable, "scripts/make_test_data.py",
                    "--out", str(data), "--reads", "4000", "--seed", "1"],
                   check=True, cwd=ROOT)
    subprocess.run([sys.executable, "run_pipeline.py",
                    "--fastq-dir", str(data / "fastq"),
                    "--samplesheet", str(data / "samplesheet.csv"),
                    "--transcriptome", str(data / "ref" / "transcriptome.fa"),
                    "--outdir", str(res), "--de-backend", "ttest",
                    "--cond-a", "control", "--cond-b", "case"],
                   check=True, cwd=ROOT)
    de = pd.read_csv(res / "de" / "de_results.csv", index_col=0)
    assert de["padj"].notna().any()
    assert (de["padj"] < 0.1).sum() >= 5, "检出的 DE 基因太少"
    # 真实 DE 基因的方向正确率应显著高于随机 (二项检验)
    truth = json.load(open(data / "ground_truth.json"))["fold_change_B_vs_A"]
    sub = de.loc[[g for g in truth if g in de.index]]
    sign_ok = (np.sign(sub["log2FoldChange"]) ==
               np.sign(np.log2(pd.Series(truth)[sub.index]))).mean()
    from scipy import stats
    assert stats.binomtest(int(sign_ok * len(sub)), len(sub), 0.5).pvalue < 0.01
    qc = pd.read_csv(res / "qc_mapping.csv")
    assert qc["mapping_rate"].min() > 0.5, "比对率过低"
    for f in ["volcano.png", "MA_plot.png", "PCA.png", "heatmap_top30.png"]:
        assert (res / "figures" / f).stat().st_size > 0
