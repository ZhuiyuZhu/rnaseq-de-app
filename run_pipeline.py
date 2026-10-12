#!/usr/bin/env python
"""RNApipe 命令行入口: fastq -> 基因计数 -> 差异分析 -> 图表"""
import argparse, logging, sys
from rnaseq_de.pipeline import run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fastq-dir", required=True, help="fastq 所在目录")
    p.add_argument("--samplesheet", required=True,
                   help="CSV: sample,fastq_1,fastq_2,condition")
    p.add_argument("--transcriptome", required=True, help="参考转录本 fasta")
    p.add_argument("--outdir", default="results")
    p.add_argument("--quant-backend", default="python",
                   choices=["python", "salmon"])
    p.add_argument("--de-backend", default="auto",
                   choices=["auto", "pydeseq2", "ttest"])
    p.add_argument("--cond-a", default=None, help="基准组名 (默认按字母序第一组)")
    p.add_argument("--cond-b", default=None, help="处理组名")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=1.0)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        stream=sys.stdout)
    run(a.fastq_dir, a.samplesheet, a.transcriptome, a.outdir,
        a.quant_backend, a.de_backend, a.cond_a, a.cond_b, a.alpha, a.lfc)


if __name__ == "__main__":
    main()
