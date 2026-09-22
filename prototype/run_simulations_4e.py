#!/usr/bin/env python3
"""Regenerate Task 0b simulation result TSVs with the current syn2bani binary.

Uses the default 4-enzyme panel (BcgI,AlfI,AloI,FalI). Writes *_4e.tsv files in
the prototype directory and preserves exact-truth labels from the manifests.
"""
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BIN = HERE.parent / "target" / "release" / "syn2bani"
REF = HERE / "mg1655.fasta"


def run(cmd, cwd=HERE, capture=True):
    print("  $ " + " ".join(str(c) for c in cmd))
    return subprocess.run(cmd, cwd=cwd, capture_output=capture, text=True,
                          check=True)


def syn2bani_tsv(queries, ref, threads=8):
    """Run syn2bani ani and return stdout TSV."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as qf:
        qf.write("\n".join(str(p) for p in queries) + "\n")
        qlist = qf.name
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as rf:
        rf.write(str(ref) + "\n")
        rlist = rf.name
    try:
        out = run([BIN, "ani", "--ql", qlist, "--rl", rlist,
                   "-p", "-t", str(threads), "--verbose"],
                  cwd=HERE).stdout
        return out
    finally:
        os.unlink(qlist)
        os.unlink(rlist)


def write_tsv(path, text):
    Path(path).write_text(text)
    print(f"  wrote {path}")


def simindel():
    print("\n=== simindel: ANI ladder with indels ===")
    simdir = HERE / "simindel"
    run([sys.executable, HERE / "simulate.py", REF, simdir, "1.0"])
    man = (simdir / "manifest.tsv").read_text().splitlines()
    queries = [simdir / line.split("\t")[3] for line in man[1:]]
    out = syn2bani_tsv(queries, simdir / "ref.fasta")
    write_tsv(HERE / "simindel_results_4e.tsv", out)


def simindel_sweep():
    print("\n=== simindel_sweep: indel-rate sweep at 95% ANI ===")
    simdir = HERE / "simindel_sweep"
    run([sys.executable, HERE / "simulate_indel_sweep.py", REF, simdir])
    man = (simdir / "manifest.tsv").read_text().splitlines()
    queries = [simdir / (line.split("\t")[0] + ".fasta") for line in man[1:]]
    out = syn2bani_tsv(queries, simdir / "ref.fasta")
    write_tsv(HERE / "simindel_sweep_results_4e.tsv", out)


def simfrag():
    print("\n=== simfrag: fragmentation of a 95% ANI draft ===")
    # Build a single 95% ANI query and name it q95.fasta so fragment.py uses the
    # expected naming convention (q95_c20, q95_c50, ...).
    tmpdir = HERE / "simfrag_tmp_src"
    run([sys.executable, HERE / "simulate.py", REF, tmpdir, "0.0"])
    src = tmpdir / "q95.fasta"
    (tmpdir / "q_ani0.9500.fasta").rename(src)
    simdir = HERE / "simfrag"
    run([sys.executable, HERE / "fragment.py", src, simdir])
    # The reference for the fragmented queries is the original MG1655 genome;
    # the source was a 95% ANI mutant, so true ANI remains 95% after fragmentation.
    import shutil
    shutil.copy(REF, simdir / "ref.fasta")
    # Use only the flipped/shuffled drafts (the forward-only controls are named *_fwd).
    queries = sorted(p for p in simdir.glob("q95_c*.fasta")
                     if "_fwd" not in p.name)
    out = syn2bani_tsv(queries, simdir / "ref.fasta")
    write_tsv(HERE / "simfrag_results_4e.tsv", out)


def simacc():
    print("\n=== simacc: accessory-fraction confound at 95% ANI ===")
    simdir = HERE / "simacc"
    # Regenerate so the result is tied to the current run, but FASTA is unchanged.
    run([sys.executable, HERE / "simulate_accessory.py", REF, simdir, "0.95"])
    queries = sorted(simdir.glob("acc*.fasta"))
    out = syn2bani_tsv(queries, simdir / "ref.fasta")
    write_tsv(HERE / "simacc_results_4e.tsv", out)


def simmosaic():
    print("\n=== simmosaic: gamma/bimodal rate heterogeneity ===")
    simdir = HERE / "simmosaic"
    run([sys.executable, HERE / "simulate_mosaic.py", REF, simdir, "5"])
    man = (simdir / "manifest.tsv").read_text().splitlines()
    queries = [simdir / (line.split("\t")[0] + ".fasta") for line in man[1:]]
    out = syn2bani_tsv(queries, simdir / "ref.fasta")
    write_tsv(HERE / "simmosaic_results_4e.tsv", out)


def main():
    if not BIN.exists():
        raise FileNotFoundError(f"syn2bani binary not found: {BIN}")
    if not REF.exists():
        raise FileNotFoundError(f"reference not found: {REF}")
    print(f"binary: {BIN}")
    print(f"reference: {REF}")
    simindel()
    simindel_sweep()
    simfrag()
    simacc()
    simmosaic()
    print("\nAll simulations complete.")


if __name__ == "__main__":
    main()
