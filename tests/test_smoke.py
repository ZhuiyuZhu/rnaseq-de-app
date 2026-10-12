"""端到端冒烟测试: 生成小数据 -> 跑流水线 -> 校验输出。"""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_pipeline_end2end(tmp_path):
    data = tmp_path / "data"
    res = tmp_path / "res"
    subprocess.run([sys.executable, "scripts/make_test_data.py",
                    "--out", str(data), "--reads", "2000", "--seed", "1"],
                   check=True, cwd=ROOT)
    subprocess.run([sys.executable, "run_pipeline.py",
                    "--fastq-dir", str(data / "fastq"),
                    "--samplesheet", str(data / "samplesheet.csv"),
                    "--transcriptome", str(data / "ref" / "transcriptome.fa"),
                    "--outdir", str(res), "--de-backend", "ttest"],
                   check=True, cwd=ROOT)
    de = pd.read_csv(res / "de" / "de_results.csv", index_col=0)
    assert de["padj"].notna().any()
    assert (de["padj"] < 0.05).sum() >= 5, "检出的 DE 基因太少"
    qc = pd.read_csv(res / "qc_mapping.csv")
    assert qc["mapping_rate"].min() > 0.5, "比对率过低"
    for f in ["volcano.png", "MA_plot.png", "PCA.png", "heatmap_top30.png"]:
        assert (res / "figures" / f).stat().st_size > 0
