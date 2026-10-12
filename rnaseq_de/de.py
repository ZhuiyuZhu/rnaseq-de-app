"""Differential-expression backends.

- `pydeseq2`: DESeq2 的官方 Python 复现（推荐，论文级方法）。
- `ttest`: 中位数比值标准化 + log2 转换 + Welch t 检验（零额外依赖，
  适合快速调试 / CI / 无 pydeseq2 的环境）。
"""
from __future__ import annotations
import numpy as np
import pandas as pd


def median_of_ratios(counts: pd.DataFrame) -> pd.DataFrame:
    """DESeq2 风格 size-factor 标准化。counts: genes x samples"""
    logc = np.log(counts.where(counts > 0))
    geo = logc.mean(axis=1)  # 只在所有样本都>0的基因上估计
    valid = geo.notna()
    sf = np.exp(logc.loc[valid].sub(geo[valid], axis=0).median(axis=0))
    return counts.div(sf, axis=1), sf


def _bh(pvals):
    p = np.asarray(pvals, dtype=float)
    n = np.sum(~np.isnan(p))
    out = np.full_like(p, np.nan)
    if n == 0:
        return out
    idx = np.where(~np.isnan(p))[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.clip(ranked, 0, 1)
    return out


def de_ttest(counts: pd.DataFrame, conditions: pd.Series,
             cond_a: str, cond_b: str) -> pd.DataFrame:
    from scipy import stats
    norm, _ = median_of_ratios(counts.round())
    logn = np.log2(norm + 1.0)
    A = logn[conditions[conditions == cond_a].index]
    B = logn[conditions[conditions == cond_b].index]
    res = pd.DataFrame(index=counts.index)
    res["baseMean"] = norm.mean(axis=1)
    res["log2FoldChange"] = B.mean(axis=1) - A.mean(axis=1)
    t, p = stats.ttest_ind(B.values, A.values, axis=1, equal_var=False)
    res["pvalue"] = p
    res["padj"] = _bh(p)
    return res.sort_values("padj")


def de_pydeseq2(counts: pd.DataFrame, conditions: pd.Series,
                cond_a: str, cond_b: str) -> pd.DataFrame:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as e:
        raise ImportError(
            "未安装 pydeseq2。请 `pip install pydeseq2`，"
            "或使用 --de-backend ttest 运行零依赖回退方案。") from e
    d = counts.round().astype(int).T  # samples x genes
    meta = pd.DataFrame({"condition": conditions}).loc[d.index]
    meta["condition"] = pd.Categorical(meta["condition"],
                                       categories=[cond_a, cond_b])
    dds = DeseqDataSet(counts=d, metadata=meta, design_factors="condition",
                       quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=("condition", cond_b, cond_a), quiet=True)
    stat.summary()
    res = stat.results_df.sort_values("padj")
    return res


def run_de(counts, conditions, cond_a, cond_b, backend="auto"):
    if backend == "auto":
        try:
            import pydeseq2  # noqa
            backend = "pydeseq2"
        except ImportError:
            backend = "ttest"
    if backend == "pydeseq2":
        return de_pydeseq2(counts, conditions, cond_a, cond_b), backend
    if backend == "ttest":
        return de_ttest(counts, conditions, cond_a, cond_b), backend
    raise ValueError(f"unknown backend: {backend}")
