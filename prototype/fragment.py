#!/usr/bin/env python3
"""Fragment a genome into draft-assembly-like contigs, preserving all sequence.

Partitioning does not change sequence content, so the true ANI against any
reference is *exactly* unchanged. Any drift as N50 falls is therefore pure
fragmentation artifact, which makes this an exact control rather than an
approximation.

Two things real assemblers do that matter here:

- **Contig orientation is arbitrary.** Roughly half of a draft assembly's
  contigs are submitted reverse-complemented relative to the reference. Tags in
  those contigs are stored as the reverse complement of their homolog, so
  without strand-canonical hashing they cannot match at all — about half the
  shared tags silently vanish. `--forward-only` disables the flips so the two
  runs can be compared.
- **Contig order is arbitrary**, so anything that assumes reference-like
  ordering breaks.

Usage:
    python3 fragment.py <genome.fasta> <outdir> [--forward-only]
"""
import os
import sys

import numpy as np

# Target mean contig lengths. "0" means leave the genome intact.
CONTIG_SIZES = [0, 500_000, 200_000, 100_000, 50_000, 20_000, 10_000, 5_000]

_COMP = bytes.maketrans(b"ACGTN", b"TGCAN")


def read_fasta_concat(path):
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    s = "".join(seqs).upper()
    return "".join(c for c in s if c in "ACGT")


def revcomp(b: bytes) -> bytes:
    return b.translate(_COMP)[::-1]


def fragment(seq: bytes, mean_len: int, rng, flip: bool):
    """Cut into contigs of length ~Uniform(0.5, 1.5) x mean_len, no sequence lost."""
    if mean_len <= 0:
        return [seq]
    contigs, i = [], 0
    n = len(seq)
    while i < n:
        ln = int(rng.integers(mean_len // 2, mean_len * 3 // 2 + 1))
        contigs.append(seq[i : i + ln])
        i += ln
    if flip:
        contigs = [revcomp(c) if rng.random() < 0.5 else c for c in contigs]
        order = rng.permutation(len(contigs))
        contigs = [contigs[i] for i in order]
    return contigs


def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    half = sum(lengths) / 2
    acc = 0
    for l in lengths:
        acc += l
        if acc >= half:
            return l
    return 0


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    forward_only = "--forward-only" in sys.argv
    if len(args) != 2:
        print(__doc__)
        return 1
    src, outdir = args
    os.makedirs(outdir, exist_ok=True)

    seq = read_fasta_concat(src).encode()
    base = os.path.basename(src).replace(".fasta", "")
    print(f"source {base}: {len(seq):,} bp  flip={'no' if forward_only else 'yes'}")

    manifest = []
    for i, size in enumerate(CONTIG_SIZES):
        rng = np.random.default_rng(4242 + i)
        contigs = fragment(seq, size, rng, flip=not forward_only)
        lens = [len(c) for c in contigs]
        assert sum(lens) == len(seq), "fragmentation must not lose sequence"
        tag = "complete" if size == 0 else f"n{size // 1000}kb"
        path = os.path.join(outdir, f"{base}.{tag}.fasta")
        with open(path, "w") as fh:
            for j, c in enumerate(contigs):
                fh.write(f">contig_{j + 1}\n")
                d = c.decode()
                for k in range(0, len(d), 80):
                    fh.write(d[k : k + 80] + "\n")
        manifest.append((tag, len(contigs), n50(lens), path))
        print(f"  {tag:<10} contigs={len(contigs):<6} N50={n50(lens):>9,} bp")

    with open(os.path.join(outdir, "manifest.tsv"), "w") as fh:
        fh.write("tag\tn_contigs\tn50\tpath\n")
        for tag, nc, n5, path in manifest:
            fh.write(f"{tag}\t{nc}\t{n5}\t{path}\n")
    print(f"\nwrote {len(manifest)} assemblies to {outdir}")
    print("True ANI is identical across all of them — no sequence was changed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
