"""FASTA/FASTQ I/O helpers (gzip transparent)."""
from __future__ import annotations
import gzip
from pathlib import Path


def open_maybe_gzip(path, mode="rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)


def read_fasta(path):
    """Yield (name, seq) tuples. Header format: '>tx_00001 gene=gene_0001 ...'"""
    name, seq = None, []
    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.upper())
        if name is not None:
            yield name, "".join(seq)


def read_fastq(path):
    """Yield (name, seq, qual) tuples."""
    with open_maybe_gzip(path) as fh:
        while True:
            name = fh.readline()
            if not name:
                return
            seq = fh.readline().strip()
            fh.readline()  # '+'
            qual = fh.readline().strip()
            yield name.strip()[1:], seq, qual
