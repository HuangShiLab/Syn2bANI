#!/usr/bin/env python3
"""Isolate the accessory-genome confound at CONSTANT true ANI.

Every query genome here has exactly the same core divergence (default 95% ANI
over homologous regions). Only the fraction of the genome replaced by
non-homologous ("accessory") sequence varies. A correct ANI estimator must
return the same value for all of them; one that is confounded by shared-content
fraction will drift monotonically with F.

Usage:
    python3 simulate_accessory.py <reference.fasta> <outdir> [core_ani]
"""
import os
import sys

import numpy as np

from simulate import BASES, read_fasta_single

ACCESSORY_FRACTIONS = [0.0, 0.10, 0.20, 0.30, 0.40, 0.50]
N_BLOCKS = 5


def substitute(seq_u8, ani, rng):
    out = seq_u8.copy()
    n = out.size
    n_sub = int(round((1.0 - ani) * n))
    pos = rng.choice(n, size=n_sub, replace=False)
    cur = out[pos]
    idx = np.searchsorted(BASES, cur)
    idx = np.where(idx >= 4, 0, idx)
    out[pos] = BASES[(idx + rng.integers(1, 4, size=n_sub, dtype=np.int64)) % 4]
    return out, n_sub


def replace_accessory(seq_u8, frac, rng, n_blocks=None):
    """Replace `frac` of the genome with shuffled sequence in N_BLOCKS chunks.

    Shuffling preserves base composition (so enzyme sites still occur at the
    normal rate and the accessory region produces its own tags, as real
    accessory genes do) while destroying all homology to the reference.
    """
    out = seq_u8.copy()
    n = out.size
    if frac <= 0:
        return out, 0
    nb = n_blocks if n_blocks is not None else N_BLOCKS
    block = int(frac * n / nb)
    total = 0
    # Evenly spaced blocks, kept away from the ends.
    stride = n // (nb + 1)
    for b in range(nb):
        start = stride * (b + 1) - block // 2
        start = max(0, min(start, n - block))
        seg = out[start : start + block].copy()
        rng.shuffle(seg)
        out[start : start + block] = seg
        total += block
    return out, total


def write_fasta(path, name, arr):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        b = arr.tobytes().decode()
        for i in range(0, len(b), 80):
            fh.write(b[i : i + 80] + "\n")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    ref_path, outdir = sys.argv[1], sys.argv[2]
    core_ani = float(sys.argv[3]) if len(sys.argv) > 3 else 0.95
    os.makedirs(outdir, exist_ok=True)

    seq = "".join(read_fasta_single(ref_path)).upper()
    seq = "".join(c for c in seq if c in "ACGT")
    seq_u8 = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    write_fasta(os.path.join(outdir, "ref.fasta"), "ref", seq_u8)
    print(f"reference length={seq_u8.size:,}  core_ani={core_ani}")

    manifest = []
    for i, frac in enumerate(ACCESSORY_FRACTIONS):
        rng = np.random.default_rng(7000 + i)
        mut, n_sub = substitute(seq_u8, core_ani, rng)
        mut, n_acc = replace_accessory(mut, frac, rng)
        # True ANI over the homologous core is unchanged by the replacement.
        homologous = seq_u8.size - n_acc
        true_ani = core_ani
        name = f"acc{frac:.2f}"
        path = os.path.join(outdir, name + ".fasta")
        write_fasta(path, name, mut)
        manifest.append((name, true_ani, path, frac))
        print(
            f"  {name}: core_ani={true_ani:.4f}  accessory={n_acc:,} bp "
            f"({n_acc/seq_u8.size:.1%})  homologous={homologous:,} bp"
        )

    with open(os.path.join(outdir, "manifest.tsv"), "w") as fh:
        fh.write("name\ttrue_ani\tpath\n")
        for name, true_ani, path, _ in manifest:
            fh.write(f"{name}\t{true_ani:.6f}\t{path}\n")
    print(f"\nwrote {len(manifest)} genomes to {outdir}")
    print("Expected: a confounded estimator drifts with accessory fraction;")
    print("a chain-restricted one stays flat at the core ANI.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
