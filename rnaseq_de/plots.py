"""Publication-style figures (matplotlib only)."""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def volcano(res: pd.DataFrame, out, alpha=0.05, lfc=1.0):
    p = -np.log10(res["pvalue"].clip(lower=1e-300))
    sig = (res["padj"] < alpha) & (res["log2FoldChange"].abs() > lfc)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(res.loc[~sig, "log2FoldChange"], p.loc[~sig], s=8, c="gray",
               alpha=0.5, label="ns")
    ax.scatter(res.loc[sig, "log2FoldChange"], p.loc[sig], s=10, c="crimson",
               alpha=0.8, label=f"DE (padj<{alpha}, |log2FC|>{lfc})")
    ax.axhline(-np.log10(0.05), ls="--", c="k", lw=0.8)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Volcano plot")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def ma_plot(res: pd.DataFrame, out, alpha=0.05):
    sig = res["padj"] < alpha
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(np.log10(res.loc[~sig, "baseMean"]),
               res.loc[~sig, "log2FoldChange"], s=8, c="gray", alpha=0.5)
    ax.scatter(np.log10(res.loc[sig, "baseMean"]),
               res.loc[sig, "log2FoldChange"], s=10, c="crimson", alpha=0.8)
    ax.axhline(0, ls="--", c="k", lw=0.8)
    ax.set_xlabel("log10(mean of normalized counts)")
    ax.set_ylabel("log2 fold change")
    ax.set_title("MA plot")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def pca_plot(norm: pd.DataFrame, conditions: pd.Series, out):
    X = np.log2(norm + 1.0).values
    X = X - X.mean(axis=1, keepdims=True)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    var = S**2 / (S**2).sum()
    pc1, pc2 = Vt[0], Vt[1]
    fig, ax = plt.subplots(figsize=(6, 5))
    for cond in conditions.unique():
        idx = [i for i, c in enumerate(conditions) if c == cond]
        ax.scatter(pc1[idx], pc2[idx], s=60, label=str(cond))
        for i in idx:
            ax.annotate(conditions.index[i], (pc1[i], pc2[i]), fontsize=7,
                        xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]*100:.1f}%)")
    ax.set_title("PCA of normalized counts")
    ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def heatmap_top(res, norm, out, n=30):
    top = res.dropna(subset=["padj"]).nsmallest(n, "padj").index
    X = np.log2(norm.loc[top] + 1.0)
    X = X.sub(X.mean(axis=1), axis=0)
    X = X.div(X.std(axis=1).replace(0, 1), axis=0)
    fig, ax = plt.subplots(figsize=(1 + 0.5 * X.shape[1], 8))
    im = ax.imshow(X.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_yticks(range(len(top)), top, fontsize=6)
    ax.set_xticks(range(X.shape[1]), X.columns, rotation=90, fontsize=7)
    fig.colorbar(im, ax=ax, label="row z-score")
    ax.set_title(f"Top {len(top)} DE genes")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
