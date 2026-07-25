#!/usr/bin/env python3
"""Diagnostic: does the residual bias scale with chain count or accessory amount?

Holds the accessory FRACTION fixed and varies only how many blocks it is split
into. More blocks -> more chain boundaries, but identical total non-homologous
content and identical true ANI. If bias tracks block count, the cause is
chain-boundary censoring (chain ends are successful anchors by construction, so
failing tags just outside the span get dropped from the denominator).
"""
import os
import sys

import numpy as np

from simulate import read_fasta_single
from simulate_accessory import substitute, replace_accessory, write_fasta

REF = sys.argv[1] if len(sys.argv) > 1 else "mg1655.fasta"
OUT = sys.argv[2] if len(sys.argv) > 2 else "simblocks"
CORE_ANI = 0.95
FRAC = 0.20
BLOCK_COUNTS = [1, 2, 5, 10, 20, 40]

os.makedirs(OUT, exist_ok=True)
seq = "".join(read_fasta_single(REF)).upper()
seq = "".join(c for c in seq if c in "ACGT")
seq_u8 = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
write_fasta(os.path.join(OUT, "ref.fasta"), "ref", seq_u8)

rows = []
for i, nb in enumerate(BLOCK_COUNTS):
    # Same seed for every block count so the substitution pattern is identical
    # and only the block geometry differs.
    rng = np.random.default_rng(31337)
    mut, n_sub = substitute(seq_u8, CORE_ANI, rng)
    mut, n_acc = replace_accessory(mut, FRAC, rng, n_blocks=nb)
    name = f"nb{nb:02d}"
    write_fasta(os.path.join(OUT, name + ".fasta"), name, mut)
    rows.append((name, nb, n_acc / seq_u8.size))
    print(f"  {name}: blocks={nb:2d}  accessory={n_acc/seq_u8.size:.1%}")

with open(os.path.join(OUT, "manifest.tsv"), "w") as fh:
    fh.write("name\ttrue_ani\tblocks\n")
    for name, nb, _ in rows:
        fh.write(f"{name}\t{CORE_ANI:.6f}\t{nb}\n")
print(f"\ntrue ANI is {CORE_ANI*100:.3f} for all {len(rows)} genomes")
