#!/usr/bin/env python3
"""How does accuracy vary with genome GC content, and does the enzyme panel matter?

Type IIB recognition sites are mostly GC-rich — of the twelve enzymes whose tags
fit the 32-base packing, only FalI (33% site GC) and PsrI (43%) lean AT, the rest
sit at 57-80%. So tag density and, worse, the *enzyme composition* of the tag pool
both shift a lot with genome GC. Since the likelihood pools all enzymes under one
divergence and one shape, a compositionally different sample per genome pair is a
plausible route to genus-specific bias.

This tests it against exact ground truth: simulate a known ANI ladder on genomes
spanning the GC range, so any difference in MAE is attributable to the genome, not
to a reference tool's own error.

Usage:
    python3 gc_bench.py <genome.fasta> [<genome.fasta> ...]
"""
import os
import statistics
import subprocess
import sys
import tempfile

BIN = os.environ.get(
    "SYN2BANI",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "target", "release", "syn2bani"),
)
PANELS = [
    ("current-4", "BcgI,AlfI,AloI,FalI"),
    ("balanced-5", "FalI,PsrI,AloI,PpiI,BcgI"),
]


def gc_and_len(path):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    n = len(seq)
    return (seq.count("G") + seq.count("C")) / n * 100, n


def run_panel(simdir, panel):
    """Return (MAE, bias, median AF, median retention, n_below_detection)."""
    queries = sorted(
        os.path.join(simdir, f) for f in os.listdir(simdir) if f.startswith("q_")
    )
    ref = os.path.join(simdir, "ref.fasta")
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as qf:
        qf.write("\n".join(queries) + "\n")
        qlist = qf.name
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as rf:
        rf.write(ref + "\n")
        rlist = rf.name
    try:
        out = subprocess.run(
            [BIN, "ani", "--ql", qlist, "--rl", rlist, "-p", "-t", "8",
             "-e", panel, "--verbose"],
            capture_output=True, text=True,
        ).stdout
    finally:
        os.unlink(qlist)
        os.unlink(rlist)

    errs, afs, rets, below = [], [], [], 0
    for i, line in enumerate(out.splitlines()):
        if i == 0:
            continue
        f = line.split("\t")
        if len(f) < 15:
            continue
        # Query genome ids are "q_ani0.9500", so truth is in the name.
        truth = float(f[0].split("ani")[1]) * 100
        errs.append(float(f[2]) - truth)
        afs.append(float(f[4]))
        rets.append(float(f[8]))
        if f[14] == "BELOW_DETECTION":
            below += 1
    if not errs:
        return None
    return (
        sum(abs(e) for e in errs) / len(errs),
        sum(errs) / len(errs),
        statistics.median(afs),
        statistics.median(rets),
        below,
        len(errs),
    )


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    here = os.path.dirname(os.path.abspath(__file__))
    print(f"{'genome':<26}{'GC%':>6}{'Mb':>6}  {'panel':<12}"
          f"{'MAE':>8}{'bias':>9}{'medAF':>7}{'medRet':>8}{'below':>7}{'n':>4}")
    print("-" * 100)
    for path in sys.argv[1:]:
        name = os.path.basename(path).replace(".fasta", "")
        gc, n = gc_and_len(path)
        simdir = os.path.join(here, f"simgc_{name}")
        if not os.path.isdir(simdir):
            subprocess.run(
                [sys.executable, os.path.join(here, "simulate.py"), path, simdir],
                capture_output=True,
            )
        for label, panel in PANELS:
            r = run_panel(simdir, panel)
            if r is None:
                print(f"{name:<26}{gc:>6.1f}{n/1e6:>6.2f}  {label:<12}      no output")
                continue
            mae, bias, af, ret, below, cnt = r
            print(f"{name:<26}{gc:>6.1f}{n/1e6:>6.2f}  {label:<12}"
                  f"{mae:>8.4f}{bias:>+9.4f}{af:>7.3f}{ret:>8.3f}{below:>7}{cnt:>4}")
    print()
    print("Ground truth is exact (counted substitutions), so MAE differences between")
    print("rows are attributable to the genome and the panel, not to a reference tool.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
