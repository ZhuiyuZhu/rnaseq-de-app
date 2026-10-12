"""Pipeline orchestration: fastq -> quant -> gene counts -> DE -> figures."""
from __future__ import annotations
import json
import logging
from pathlib import Path
import pandas as pd

from .quant import (build_kmer_index, parse_tx2gene_from_fasta,
                    quant_sample_python, tx_to_gene_counts)
from .de import run_de, median_of_ratios
from . import plots

log = logging.getLogger("rnaseq_de")


def load_samplesheet(path):
    df = pd.read_csv(path)
    need = {"sample", "fastq_1", "fastq_2", "condition"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"samplesheet 缺少列: {missing}；需要 {sorted(need)}")
    if df["sample"].duplicated().any():
        raise ValueError("samplesheet 中 sample 名重复")
    return df


def run(fastq_dir, samplesheet, transcriptome, outdir,
        quant_backend="python", de_backend="auto",
        cond_a=None, cond_b=None, alpha=0.05, lfc=1.0):
    outdir = Path(outdir)
    (outdir / "counts").mkdir(parents=True, exist_ok=True)
    (outdir / "de").mkdir(parents=True, exist_ok=True)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)

    ss = load_samplesheet(samplesheet)
    conds = sorted(ss["condition"].unique())
    if len(conds) != 2:
        raise ValueError(f"本版本仅支持两组对比，当前有 {len(conds)} 组: {conds}")
    cond_a = cond_a or conds[0]
    cond_b = cond_b or conds[1]
    log.info("对比: %s (base) vs %s", cond_a, cond_b)

    tx_names, gene_of = parse_tx2gene_from_fasta(transcriptome)
    index, _ = build_kmer_index(transcriptome)

    count_tbl, qc_rows = {}, []
    for _, row in ss.iterrows():
        fq1 = Path(fastq_dir) / row["fastq_1"]
        fq2 = Path(fastq_dir) / row["fastq_2"]
        log.info("定量 %s ...", row["sample"])
        txc, stats = quant_sample_python(fq1, fq2, index, len(tx_names))
        genes, gc = tx_to_gene_counts(txc, tx_names, gene_of)
        count_tbl[row["sample"]] = gc
        qc_rows.append({"sample": row["sample"], **stats})
        log.info("  比对率 %.1f%%", stats["mapping_rate"] * 100)

    counts = pd.DataFrame(count_tbl, index=genes)
    conditions = ss.set_index("sample")["condition"]
    counts = counts[conditions.index]  # 对齐顺序

    pd.DataFrame(qc_rows).to_csv(outdir / "qc_mapping.csv", index=False)
    counts.to_csv(outdir / "counts" / "raw_counts.csv")
    norm, sf = median_of_ratios(counts.round())
    norm.to_csv(outdir / "counts" / "normalized_counts.csv")
    sf.rename("size_factor").to_csv(outdir / "counts" / "size_factors.csv")

    log.info("差异分析 (backend=%s) ...", de_backend)
    res, used_backend = run_de(counts, conditions, cond_a, cond_b, de_backend)
    log.info("实际使用 backend: %s", used_backend)
    res.to_csv(outdir / "de" / "de_results.csv")

    plots.volcano(res, outdir / "figures" / "volcano.png", alpha, lfc)
    plots.ma_plot(res, outdir / "figures" / "MA_plot.png", alpha)
    plots.pca_plot(norm, conditions, outdir / "figures" / "PCA.png")
    plots.heatmap_top(res, norm, outdir / "figures" / "heatmap_top30.png")
    (outdir / "run_info.json").write_text(json.dumps({
        "quant_backend": quant_backend, "de_backend": used_backend,
        "cond_a": cond_a, "cond_b": cond_b, "alpha": alpha, "lfc": lfc,
        "n_genes": int(counts.shape[0]), "n_samples": int(counts.shape[1]),
    }, indent=2))
    log.info("完成。结果在 %s", outdir)
    return res, norm
