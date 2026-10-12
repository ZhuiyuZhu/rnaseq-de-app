"""Transcript quantification backends.

- `python`: dependency-free k-mer pseudo-aligner (kallisto/salmon 的简化版)。
  适用于小参考、调试和教学；对真实大数据请用 salmon。
- `salmon`: 调用 salmon 二进制（需要自行安装，推荐 Linux/macOS/conda）。
"""
from __future__ import annotations
import subprocess
from pathlib import Path
import numpy as np
from .io import read_fasta, read_fastq

KMER = 21
MAX_KMER_FREQ = 10  # 出现在过多转录本中的 k-mer 视为重复序列，忽略


def parse_tx2gene_from_fasta(fasta_path):
    """fasta header: '>tx_00001 gene=gene_0001' -> (tx_names, gene_of)"""
    tx_names, gene_of = [], []
    for header_line in _headers(fasta_path):
        parts = header_line.split()
        tx = parts[0]
        gene = None
        for p in parts[1:]:
            if p.startswith("gene="):
                gene = p.split("=", 1)[1]
        if gene is None:
            gene = tx  # 单转录本基因：基因名=转录本名
        tx_names.append(tx)
        gene_of.append(gene)
    return tx_names, gene_of


def _headers(fasta_path):
    from .io import open_maybe_gzip
    with open_maybe_gzip(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].strip()


COMP = str.maketrans("ACGT", "TGCA")


def build_kmer_index(fasta_path, k=KMER):
    """k-mer -> list[transcript_idx]，正/反两条链都建索引，并过滤高频 k-mer。

    R2 read 是转录本的反向互补序列，必须索引 revcomp 链才能比对。
    """
    raw = {}
    seqs = []
    for tx, seq in read_fasta(fasta_path):
        idx = len(seqs)
        seqs.append(seq)
        rc = seq.translate(COMP)[::-1]
        for s in (seq, rc):
            for i in range(len(s) - k + 1):
                raw.setdefault(s[i:i + k], []).append(idx)
    index = {km: tuple(v) for km, v in raw.items() if len(v) <= MAX_KMER_FREQ}
    return index, seqs


def _mate_votes(seq, index, k=KMER):
    """一条 read 的转录本投票：返回 {tx_idx: hit_count}"""
    hits = {}
    n = len(seq) - k + 1
    if n <= 0:
        return hits
    for i in range(n):
        for tx in index.get(seq[i:i + k], ()):
            hits[tx] = hits.get(tx, 0) + 1
    if not hits:
        return hits
    top = max(hits.values())
    if top < max(5, 0.3 * n):      # 至少有 5 个命中，且达到理论最大的一定比例
        return {}
    cutoff = 0.6 * top
    return {tx: c for tx, c in hits.items() if c >= cutoff}


def quant_sample_python(fq1, fq2, index, n_tx, k=KMER):
    """简化版 pair-end pseudo-alignment，分数计数（ties 均分）。

    基因层面分析时，同一基因各 isoform 间的均分不影响基因总计数。
    返回 (tx_counts, stats)。
    """
    counts = np.zeros(n_tx)
    n_pairs = n_mapped = 0
    for (n1, s1, q1), (n2, s2, q2) in zip(read_fastq(fq1), read_fastq(fq2)):
        n_pairs += 1
        v1 = _mate_votes(s1, index, k)
        v2 = _mate_votes(s2, index, k)
        if not v1 or not v2:
            continue
        common = set(v1) & set(v2)
        if not common:
            continue
        score = {t: v1.get(t, 0) + v2.get(t, 0) for t in common}
        best = max(score.values())
        winners = [t for t, sc in score.items() if sc == best]
        for t in winners:
            counts[t] += 1.0 / len(winners)
        n_mapped += 1
    stats = {"pairs": n_pairs, "mapped": n_mapped,
             "mapping_rate": n_mapped / max(n_pairs, 1)}
    return counts, stats


def quant_sample_salmon(fq1, fq2, transcriptome, outdir, threads=4):
    """调用 salmon 进行定量（需已安装 salmon 并加入 PATH）。"""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    idx = outdir.parent / "salmon_index"
    if not idx.exists():
        subprocess.run(["salmon", "index", "-t", str(transcriptome),
                        "-i", str(idx), "-k", str(KMER)], check=True)
    cmd = ["salmon", "quant", "-i", str(idx), "-l", "A",
           "-1", str(fq1), "-2", str(fq2), "-o", str(outdir),
           "--threads", str(threads), "--seqBias", "--gcBias"]
    subprocess.run(cmd, check=True)
    tpm = {}
    with open(outdir / "quant.sf") as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            tpm[f[0]] = float(f[3])
    return tpm


def tx_to_gene_counts(tx_counts, tx_names, gene_of):
    """把转录本计数(可含分数)聚合成基因计数。"""
    genes = list(dict.fromkeys(gene_of))
    gidx = {g: i for i, g in enumerate(genes)}
    gcounts = np.zeros(len(genes))
    for c, g in zip(tx_counts, gene_of):
        gcounts[gidx[g]] += c
    return genes, gcounts
