#!/usr/bin/env python
"""用 ground_truth.json 验证差异分析结果: 召回率 + 背景误报率。"""
import json
import sys
import numpy as np
import pandas as pd

res_path = sys.argv[1] if len(sys.argv) > 1 else "results/de/de_results.csv"
truth_path = sys.argv[2] if len(sys.argv) > 2 else "data/ground_truth.json"
res = pd.read_csv(res_path, index_col=0)
truth = json.load(open(truth_path))["fold_change_B_vs_A"]
v = pd.DataFrame(
    [(g, np.log2(fc), res.loc[g, "log2FoldChange"], res.loc[g, "padj"])
     for g, fc in truth.items() if g in res.index],
    columns=["gene", "true_log2fc", "est_log2fc", "padj"])
ALPHA = 0.1  # RNA-seq 常用 FDR 截断
hit = (v["padj"] < ALPHA) & (np.sign(v["est_log2fc"]) == np.sign(v["true_log2fc"]))
print(f"DE 基因召回率 (padj<{ALPHA} 且方向正确): {hit.mean():.1%} ({hit.sum()}/{len(v)})")
nulls = res.index.difference(list(truth))
print(f"非 DE 基因误报率 (padj<{ALPHA}): {(res.loc[nulls, 'padj'] < ALPHA).mean():.1%} "
      f"({int((res.loc[nulls, 'padj'] < ALPHA).sum())}/{len(nulls)})")
ok = hit.mean() >= 0.8 and (res.loc[nulls, "padj"] < ALPHA).mean() <= 0.1
print("验证:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
