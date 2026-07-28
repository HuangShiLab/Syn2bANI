#!/usr/bin/env python3
"""Measure each enzyme's systematic bias against exactly-known truth.

Generates replicate genomes at one true ANI with independent mutation draws, so
anything that survives averaging is systematic rather than sampling noise. The
divergence is uniform per site, which means the enzymes cannot legitimately
disagree: any reproducible spread between them is a property of the enzyme's
handling, not of the genomes.

Also calibrates the reported `std_err` against the observed spread across
replicates — if those disagree, `enzyme_chi2` has no meaningful null.

Usage:
    python3 enzyme_bias.py [reference.fasta] [true_ani] [n_replicates]
"""
import os
import statistics
import subprocess
import sys
import tempfile

import numpy as np

from simulate import mutate, read_fasta_single

BIN = os.environ.get("SYN2BANI", "../target/release/syn2bani")
PANEL = os.environ.get(
    "PANEL", "BcgI,AlfI,AloI,FalI,BplI,Bsp24I,PpiI,PsrI,BsaXI,CjeI,CjePI"
)


def main():
    ref = sys.argv[1] if len(sys.argv) > 1 else "mg1655.fasta"
    true_ani = float(sys.argv[2]) if len(sys.argv) > 2 else 0.95
    n_rep = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    outdir = "simrep"

    seq = "".join(read_fasta_single(ref)).upper()
    seq = "".join(c for c in seq if c in "ACGT")
    u8 = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    os.makedirs(outdir, exist_ok=True)

    def write(path, name, arr):
        with open(path, "w") as fh:
            fh.write(f">{name}\n")
            b = arr.tobytes().decode()
            for i in range(0, len(b), 80):
                fh.write(b[i : i + 80] + "\n")

    write(os.path.join(outdir, "ref.fasta"), "ref", u8)
    paths = []
    for r in range(n_rep):
        p = os.path.join(outdir, f"rep{r:02d}.fasta")
        paths.append(p)
        if os.path.exists(p):
            continue
        rng = np.random.default_rng(90000 + r)
        mut, _ = mutate(u8, true_ani, rng, inversion=False, indel_rate=0.0)
        write(p, f"rep{r:02d}", mut)

    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write("\n".join(paths) + "\n")
        ql = f.name
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(os.path.join(outdir, "ref.fasta") + "\n")
        rl = f.name
    out = subprocess.run(
        [BIN, "ani", "--ql", ql, "--rl", rl, "-p", "-t", "8", "-e", PANEL, "--verbose"],
        capture_output=True,
        text=True,
    ).stdout
    os.unlink(ql)
    os.unlink(rl)

    anis, ses, per = [], [], {}
    for i, line in enumerate(out.splitlines()):
        if i == 0:
            continue
        f = line.split("\t")
        anis.append(float(f[2]))
        ses.append(float(f[6]))
        for kv in f[13].split(","):
            if ":" in kv:
                name, v = kv.split(":")
                per.setdefault(name, []).append(float(v))

    if len(anis) < 2:
        print("not enough replicates produced a result", file=sys.stderr)
        return 1

    obs = statistics.stdev(anis)
    rep = statistics.mean(ses)
    print(f"{len(anis)} replicates at true ANI {true_ani*100:.2f}, panel: {PANEL}")
    print()
    print("std_err calibration (a ratio near 1 means enzyme_chi2 has a valid null)")
    print(f"  observed SD across replicates : {obs:.4f}")
    print(f"  mean reported std_err         : {rep:.4f}")
    print(f"  ratio                         : {obs/rep:.2f}x")
    print()
    print("per-enzyme systematic bias (uniform divergence, so truth is the same")
    print("for every enzyme; sigma is the bias divided by the SE of its mean)")
    print(f"  {'enzyme':<9}{'mean':>9}{'bias':>9}{'SD':>9}{'sigma':>8}")
    rows = []
    for name, v in per.items():
        if len(v) < 2:
            continue
        sd = statistics.stdev(v)
        bias = statistics.mean(v) - true_ani * 100
        sem = sd / (len(v) ** 0.5)
        rows.append((abs(bias / sem) if sem else 0.0, name, statistics.mean(v), bias, sd))
    for sigma, name, mean, bias, sd in sorted(rows, reverse=True):
        print(f"  {name:<9}{mean:>9.3f}{bias:>+9.3f}{sd:>9.4f}{sigma:>8.1f}")
    print()
    print("Enzymes above ~3 sigma are reproducibly biased, not noisy. Note that a")
    print("low bias here does NOT predict good real-genome behaviour: see")
    print("ALGORITHM_MLE.md 4.11.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
