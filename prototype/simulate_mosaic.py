#!/usr/bin/env python3
"""Simulate mosaic divergence — the case the uniform-rate simulations never test.

`simulate.py` mutates every site at the same rate, which is why it reports MAE
around 0.1%. Real genome pairs are mosaics: conserved core under purifying
selection next to far more divergent segments. That is exactly what the gamma
rate-heterogeneity model in `mle.rs` claims to handle, and it has never been
checked against ground truth, because generating it needs this script.

Ground truth stays exact: substitutions are counted, so true ANI over the whole
genome is `1 - n_subs / length` no matter how they are distributed.

Two regimes:

* **gamma** — per-block rate multipliers drawn from Gamma(alpha, alpha), which is
  precisely the model `estimate_heterogeneous` fits. If it cannot recover the
  truth here, the estimator is wrong rather than the data.
* **bimodal** — a conserved fraction at one rate and the rest at a much higher
  rate, i.e. deliberately *not* the assumed distribution. This is the
  misspecification test.

Usage:
    python3 simulate_mosaic.py <reference.fasta> <outdir> [block_kb]
"""
import os
import sys

import numpy as np

from simulate import BASES, read_fasta_single
from simulate_accessory import write_fasta

# (label, regime, mean ANI, shape or conserved-fraction)
CASES = [
    ("gamma_a0.5_ani95", "gamma", 0.95, 0.5),
    ("gamma_a1.0_ani95", "gamma", 0.95, 1.0),
    ("gamma_a2.0_ani95", "gamma", 0.95, 2.0),
    ("gamma_a0.5_ani90", "gamma", 0.90, 0.5),
    ("gamma_a1.0_ani90", "gamma", 0.90, 1.0),
    ("gamma_a1.0_ani98", "gamma", 0.98, 1.0),
    ("bimodal_70core_ani95", "bimodal", 0.95, 0.70),
    ("bimodal_50core_ani95", "bimodal", 0.95, 0.50),
    ("bimodal_70core_ani90", "bimodal", 0.90, 0.70),
]


def mutate_blockwise(seq_u8, rates, block, rng):
    """Apply per-block substitution rates. Returns (mutated, total substitutions).

    `rates[i]` is the per-site substitution probability for block i.
    """
    out = seq_u8.copy()
    n = out.size
    total = 0
    for i, rate in enumerate(rates):
        lo = i * block
        hi = min(lo + block, n)
        if hi <= lo or rate <= 0:
            continue
        span = hi - lo
        k = int(round(rate * span))
        if k <= 0:
            continue
        k = min(k, span)
        pos = lo + rng.choice(span, size=k, replace=False)
        cur = out[pos]
        idx = np.searchsorted(BASES, cur)
        idx = np.where(idx >= 4, 0, idx)
        out[pos] = BASES[(idx + rng.integers(1, 4, size=k, dtype=np.int64)) % 4]
        total += k
    return out, total


def block_rates(regime, mean_ani, param, n_blocks, rng):
    """Per-block substitution probabilities with the requested mean."""
    mean_rate = 1.0 - mean_ani
    if regime == "gamma":
        alpha = param
        mult = rng.gamma(shape=alpha, scale=1.0 / alpha, size=n_blocks)
    else:
        conserved = param
        # Conserved blocks evolve 10x slower than the divergent remainder; the
        # mixture is then rescaled to the requested mean.
        is_cons = rng.random(n_blocks) < conserved
        mult = np.where(is_cons, 0.1, 1.0)
    mult = mult / mult.mean()               # exact mean multiplier of 1
    rates = np.clip(mean_rate * mult, 0.0, 0.75)
    return rates


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    ref_path, outdir = sys.argv[1], sys.argv[2]
    block = int(float(sys.argv[3]) * 1000) if len(sys.argv) > 3 else 5000

    os.makedirs(outdir, exist_ok=True)
    seq = "".join(read_fasta_single(ref_path)).upper()
    seq = "".join(c for c in seq if c in "ACGT")
    seq_u8 = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    write_fasta(os.path.join(outdir, "ref.fasta"), "ref", seq_u8)
    n_blocks = (seq_u8.size + block - 1) // block
    print(f"reference {seq_u8.size:,} bp, {n_blocks} blocks of {block:,} bp")

    manifest = []
    for i, (label, regime, mean_ani, param) in enumerate(CASES):
        rng = np.random.default_rng(90000 + i)
        rates = block_rates(regime, mean_ani, param, n_blocks, rng)
        mut, n_sub = mutate_blockwise(seq_u8, rates, block, rng)
        true_ani = 1.0 - n_sub / seq_u8.size
        # Name carries the truth so the evaluator can read it back.
        name = f"q_ani{true_ani:.4f}__{label}"
        write_fasta(os.path.join(outdir, name + ".fasta"), name, mut)
        manifest.append((name, true_ani, regime, param))
        print(f"  {label:<22} true_ani={true_ani:.4f}  "
              f"rate spread p10-p90 = {np.percentile(rates,10):.4f}-{np.percentile(rates,90):.4f}")

    with open(os.path.join(outdir, "manifest.tsv"), "w") as fh:
        fh.write("name\ttrue_ani\tregime\tparam\n")
        for name, true_ani, regime, param in manifest:
            fh.write(f"{name}\t{true_ani:.6f}\t{regime}\t{param}\n")
    print(f"\nwrote {len(manifest)} genomes to {outdir}")
    print("True ANI is exact in every case; only the spatial distribution differs.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
