#!/usr/bin/env python
"""生成自带"已知差异"的模拟验证数据。

- 参考转录组: N_GENES 个基因，每基因 1-3 个 isoform；
  同基因 isoform 共享一段共同主体 + 不同的 5'/3' 端（模拟真实可变剪接）。
- 模拟 reads: 负二项分布抽样基因计数，片段均匀取样，PE60，含少量测序错误。
- ground truth: 随机挑 N_DE 个基因在 B 组上调/下调 (2~8x)，写入 ground_truth.json。

用法: python scripts/make_test_data.py --out data --reads 8000 --seed 42
"""
import argparse, gzip, json
import pandas as pd
from pathlib import Path
import numpy as np

BASES = np.array(list("ACGT"))
COMP = str.maketrans("ACGT", "TGCA")
N_GENES, N_DE = 120, 30


def rand_seq(rng, n):
    return "".join(rng.choice(BASES, size=n))


def make_transcriptome(rng, out_fa):
    lines, meta = [], []
    for g in range(N_GENES):
        gene = f"gene_{g:04d}"
        n_iso = rng.integers(1, 4)
        L = int(rng.integers(900, 2400))
        common = rand_seq(rng, max(L - 400, 400))
        for j in range(n_iso):
            tx = f"tx_{g:04d}_{j}"
            l5 = int(rng.integers(150, 250)) if j > 0 else 0
            l3 = int(rng.integers(150, 250)) if j > 0 else 0
            if len(common) + l5 + l3 < 800:
                l5 = l3 = 100
            seq = rand_seq(rng, l5) + common + rand_seq(rng, l3)
            lines.append(f">{tx} gene={gene}\n{seq}")
            meta.append({"tx": tx, "gene": gene, "len": len(seq)})
    Path(out_fa).write_text("\n".join(lines) + "\n")
    return meta


def revcomp(s):
    return s.translate(COMP)[::-1]


def true_expression(meta, rng, frag_len=250):
    """所有样本共享的真实表达水平(生物学重复只共享真值, 差异仅来自抽样)。"""
    genes = sorted({m["gene"] for m in meta})
    by_gene = {g: [m for m in meta if m["gene"] == g] for g in genes}
    expr = dict(zip(genes, rng.lognormal(mean=8.0, sigma=1.0,
                                         size=len(genes))))
    eff = {g: float(np.mean([max(m["len"] - frag_len + 1, 1)
                             for m in by_gene[g]])) for g in genes}
    return expr, eff


def simulate(meta, rng, group_fc, n_reads, expr, eff,
             read_len=60, frag_len=250, err=0.005):
    genes = sorted({m["gene"] for m in meta})
    by_gene = {g: [m for m in meta if m["gene"] == g] for g in genes}
    base = np.array([expr[g] * eff[g] / 1000.0 for g in genes])
    base *= n_reads / base.sum()                # 基准文库(无DE)只做一次归一化
    depth = rng.uniform(0.95, 1.05)             # 测序深度微抖动
    lam = base * group_fc * depth               # DE 改变组成比例, 不回压文库总量
    disp = 0.05
    counts = rng.negative_binomial(1 / disp, 1 / (1 + lam * disp))
    out = []
    for g, c in zip(genes, counts):
        txs = by_gene[g]
        w = np.array([m["len"] for m in txs], float)
        w /= w.sum()
        for _ in range(c):
            m = txs[rng.choice(len(txs), p=w)]
            seq = None
            for line_no, m2 in enumerate(meta):
                if m2 is m:
                    seq = SEQ_CACHE[line_no]
                    break
            if len(seq) < frag_len:
                continue
            p = int(rng.integers(0, len(seq) - frag_len + 1))
            r1, r2 = seq[p:p + read_len], revcomp(seq[p + frag_len - read_len:
                                                     p + frag_len])
            for r in (r1, r2):
                arr = np.array(list(r))
                mask = rng.random(len(arr)) < err
                if mask.any():
                    arr[mask] = rng.choice(BASES, size=mask.sum())
                    if r is r1:
                        r1 = "".join(arr)
                    else:
                        r2 = "".join(arr)
            out.append((r1, r2))
    return out, counts


def write_fastq(pairs, f1, f2, retries=3):
    """顺序写 R1 再写 R2, 关闭失败(网络文件系统抖动)时重试。"""
    for attempt in range(retries):
        try:
            with gzip.open(f1, "wt") as a:
                for i, (r1, _) in enumerate(pairs):
                    a.write(f"@read{i}/1\n{r1}\n+\n{'I' * len(r1)}\n")
            with gzip.open(f2, "wt") as b:
                for i, (_, r2) in enumerate(pairs):
                    b.write(f"@read{i}/2\n{r2}\n+\n{'I' * len(r2)}\n")
            return
        except OSError:
            if attempt == retries - 1:
                raise


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data")
    ap.add_argument("--reads", type=int, default=40000, help="每样本 read pairs")
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    out = Path(a.out)
    (out / "ref").mkdir(parents=True, exist_ok=True)
    (out / "fastq").mkdir(parents=True, exist_ok=True)

    meta = make_transcriptome(rng, out / "ref" / "transcriptome.fa")
    global SEQ_CACHE
    SEQ_CACHE = []
    cur = None
    for line in open(out / "ref" / "transcriptome.fa"):
        line = line.strip()
        if line.startswith(">"):
            if cur:
                SEQ_CACHE.append("".join(cur))
            cur = []
        else:
            cur.append(line)
    SEQ_CACHE.append("".join(cur))

    genes = sorted({m["gene"] for m in meta})
    expr, eff = true_expression(meta, rng)
    baseline = np.array([expr[g] * eff[g] for g in genes])
    cand = np.argsort(baseline)[-(len(genes) // 3):]  # 高表达的前 1/3
    de_idx = rng.choice(cand, size=N_DE, replace=False)
    fc = {g: 1.0 for g in genes}
    for i, gi in enumerate(de_idx):
        fc[genes[gi]] = float(2.0 ** rng.uniform(1, 3)) * (1 if i % 2 == 0 else -1)
    # 真实 fold change (B vs A): fc=+x -> x 倍上调; fc=-x -> 1/x 倍下调
    truth = {g: (f if f > 0 else 1.0 / abs(f)) for g, f in fc.items()
             if f != 1.0}
    (out / "ground_truth.json").write_text(json.dumps({
        "group_A": "control", "group_B": "case",
        "fold_change_B_vs_A": truth}, indent=2))

    case_vec = np.array([f if f > 0 else 1.0 / abs(f) for f in
                         (fc[g] for g in genes)])
    ones = np.ones(len(genes))
    expected = {}
    rows = []
    for cond in ["control", "case"]:
        gvec = case_vec if cond == "case" else ones  # 对照组不加倍数
        for rep in range(3):
            sample = f"{cond}_{rep+1}"
            print(f"模拟 {sample} ...")
            pairs, sim_counts = simulate(meta, rng, gvec, a.reads, expr, eff)
            expected[sample] = sim_counts
            f1 = f"{sample}_R1.fastq.gz"
            f2 = f"{sample}_R2.fastq.gz"
            write_fastq(pairs, out / "fastq" / f1, out / "fastq" / f2)
            rows.append({"sample": sample, "fastq_1": f1, "fastq_2": f2,
                         "condition": cond})
    import csv
    exp = pd.DataFrame(expected, index=genes)
    exp.to_csv(out / "expected_counts.csv")
    with open(out / "samplesheet.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2",
                                           "condition"])
        w.writeheader()
        w.writerows(rows)
    print("完成 ->", out)


if __name__ == "__main__":
    main()
