#!/usr/bin/env python3
"""Generate mutated genomes at exactly known ANI from a reference FASTA.

Ground truth is exact by construction: we substitute a counted number of bases,
so true ANI = 1 - n_subs / genome_length. Optionally add an inversion and
indels to exercise chaining and gap arithmetic.

Usage:
    python3 simulate.py <reference.fasta> <outdir>
"""
import os
import sys

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)

# ANI levels to generate. Dense in the strain range, plus the mid-ANI band the
# HPC validation used (0.85-0.95) so the two regimes can be compared directly.
ANI_LEVELS = [0.85, 0.88, 0.90, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999]


def read_fasta_single(path):
    """Read a FASTA file, concatenating ALL contigs (unlike the Rust extractor)."""
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
    return seqs


def mutate(seq_u8, ani, rng, inversion=False, indel_rate=0.0):
    """Apply exactly round((1-ani)*L) substitutions, then optional SV/indels.

    Returns (mutated_uint8_array, n_substitutions).
    """
    out = seq_u8.copy()
    n = out.size
    n_sub = int(round((1.0 - ani) * n))
    if n_sub > 0:
        pos = rng.choice(n, size=n_sub, replace=False)
        # Draw a replacement base that differs from the original.
        cur = out[pos]
        offset = rng.integers(1, 4, size=n_sub, dtype=np.int64)
        idx = np.searchsorted(BASES, cur)
        idx = np.where(idx >= 4, 0, idx)
        out[pos] = BASES[(idx + offset) % 4]

    if inversion:
        # Invert a 400 kb segment at ~1/3 of the genome to test collinear
        # chaining across an orientation flip.
        lo = n // 3
        hi = min(lo + 400_000, n)
        seg = out[lo:hi][::-1]
        comp = {b"A"[0]: b"T"[0], b"C"[0]: b"G"[0], b"G"[0]: b"C"[0], b"T"[0]: b"A"[0]}
        lut = np.arange(256, dtype=np.uint8)
        for k, v in comp.items():
            lut[k] = v
        out[lo:hi] = lut[seg]

    if indel_rate > 0.0:
        # Delete a handful of 200-2000 bp blocks to make gap arithmetic disagree.
        n_indel = max(1, int(indel_rate * n / 100_000))
        starts = np.sort(rng.choice(n - 3000, size=n_indel, replace=False))
        keep = np.ones(out.size, dtype=bool)
        for s in starts:
            ln = int(rng.integers(200, 2000))
            keep[s : s + ln] = False
        out = out[keep]

    return out, n_sub


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 1
    ref_path, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)

    contigs = read_fasta_single(ref_path)
    seq = "".join(contigs).upper()
    seq = "".join(c for c in seq if c in "ACGT")
    seq_u8 = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    print(f"reference: {ref_path}  length={seq_u8.size:,}")

    ref_out = os.path.join(outdir, "ref.fasta")
    with open(ref_out, "w") as fh:
        fh.write(">ref\n")
        b = seq_u8.tobytes().decode()
        for i in range(0, len(b), 80):
            fh.write(b[i : i + 80] + "\n")

    manifest = []
    for i, ani in enumerate(ANI_LEVELS):
        rng = np.random.default_rng(1000 + i)
        mut, n_sub = mutate(seq_u8, ani, rng, inversion=True, indel_rate=0.0)
        true_ani = 1.0 - n_sub / seq_u8.size
        name = f"q_ani{ani:.4f}"
        path = os.path.join(outdir, name + ".fasta")
        with open(path, "w") as fh:
            fh.write(f">{name}\n")
            b = mut.tobytes().decode()
            for j in range(0, len(b), 80):
                fh.write(b[j : j + 80] + "\n")
        manifest.append((name, true_ani, path))
        print(f"  {name}: true_ani={true_ani:.6f}  n_sub={n_sub:,}  len={mut.size:,}")

    with open(os.path.join(outdir, "manifest.tsv"), "w") as fh:
        fh.write("name\ttrue_ani\tpath\n")
        for name, true_ani, path in manifest:
            fh.write(f"{name}\t{true_ani:.6f}\t{path}\n")
    print(f"\nwrote {len(manifest)} query genomes + ref.fasta to {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
