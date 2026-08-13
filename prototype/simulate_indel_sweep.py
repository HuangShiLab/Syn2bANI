#!/usr/bin/env python3
"""Indel-rate sweep at fixed substitution divergence (true ANI 95%).

Complements simulate.py's ANI ladder: here the substitution count is identical
in every genome and only the number of 200-2000 bp deletions varies, so any
change in the estimate isolates the effect of gap arithmetic / chaining across
indels. True ANI is exactly 1 - n_sub / L for every row; the deleted fraction
is the AF ground truth.

Usage:
    python3 simulate_indel_sweep.py <reference.fasta> <outdir>
"""
import os
import sys

import numpy as np

from simulate import mutate, read_fasta_single

INDEL_RATES = [0.0, 0.5, 1.0, 2.0, 4.0]  # deletions per 100 kb
ANI = 0.95


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 1
    ref_path, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    seq = "".join(read_fasta_single(ref_path)).upper()
    seq = "".join(c for c in seq if c in "ACGT")
    seq_u8 = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    n = seq_u8.size

    ref_out = os.path.join(outdir, "ref.fasta")
    with open(ref_out, "w") as fh:
        fh.write(">ref\n")
        b = seq_u8.tobytes().decode()
        for i in range(0, len(b), 80):
            fh.write(b[i : i + 80] + "\n")

    rows = []
    for i, rate in enumerate(INDEL_RATES):
        rng = np.random.default_rng(7000 + i)
        mut, n_sub = mutate(seq_u8, ANI, rng, inversion=True, indel_rate=rate)
        true_ani = 1.0 - n_sub / n
        kept = mut.size / n
        name = f"q_indel{rate:g}"
        path = os.path.join(outdir, name + ".fasta")
        with open(path, "w") as fh:
            fh.write(f">{name}\n")
            b = mut.tobytes().decode()
            for j in range(0, len(b), 80):
                fh.write(b[j : j + 80] + "\n")
        rows.append((name, rate, true_ani, kept))
        print(f"  {name}: true_ani={true_ani:.6f} deleted={1-kept:.4%} len={mut.size:,}")

    with open(os.path.join(outdir, "manifest.tsv"), "w") as fh:
        fh.write("name\tindel_rate\ttrue_ani\tkept_frac\n")
        for name, rate, true_ani, kept in rows:
            fh.write(f"{name}\t{rate}\t{true_ani:.6f}\t{kept:.4f}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
