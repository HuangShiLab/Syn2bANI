#!/usr/bin/env python3
"""Join syn2bani / skani / FastANI output for the real-genome benchmark.

There is no ground truth for real genome pairs, so this reports agreement.
Pairs flagged BELOW_DETECTION are excluded from the summary: that is the regime
where retention is too low for either rate model, and skani declines to report
those pairs at all, so including them would compare an extrapolation against
nothing.

Run from inside the benchmark directory created by realgenome_bench.sh.
"""
import csv
import math
import os

name_of_acc = {}
with open("accessions.tsv") as fh:
    for line in fh:
        acc, name = line.rstrip("\n").split("\t")
        name_of_acc[acc] = name
        name_of_acc[acc.split(".")[0]] = name


def acc_from_ena(s):
    """'ENA|U00096|U00096.3' -> 'U00096.3'"""
    parts = s.split("|")
    return parts[-1] if len(parts) > 1 else s


def basename(p):
    return os.path.basename(p).replace(".fasta", "")


syn = {}
with open("syn2bani.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        acc = acc_from_ena(row["query"])
        n = name_of_acc.get(acc) or name_of_acc.get(acc.split(".")[0])
        if n:
            syn[n] = row

skani = {}
if os.path.exists("skani.tsv"):
    with open("skani.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            skani[basename(row["Query_file"])] = float(row["ANI"])

fastani = {}
if os.path.exists("fastani.tsv"):
    with open("fastani.tsv") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            fastani[basename(f[0])] = float(f[2])

order = sorted(syn, key=lambda n: -float(syn[n]["ani_uniform"]))
print(
    f"{'genome':<30}{'ani':>7}{'unif':>8}{'skani':>7}{'fANI':>7}"
    f"{'shape':>7}{'ret':>6}{'AF':>6}  flag"
)
print("-" * 95)

reported = []
for n in order:
    s = syn[n]
    het = float(s["ani"])
    unif = float(s["ani_uniform"])
    shape = float(s.get("het_shape", "nan"))
    shape_str = "unif" if not math.isfinite(shape) or shape > 1e5 else f"{shape:.2f}"
    a_sk = skani.get(n)
    a_fa = fastani.get(n)
    ret = float(s.get("retention", "nan"))
    print(
        f"{n:<30}{het:>7.2f}{unif:>8.2f}"
        f"{(f'{a_sk:.2f}' if a_sk else '-'):>7}{(f'{a_fa:.2f}' if a_fa else '-'):>7}"
        f"{shape_str:>7}{ret:>6.2f}{float(s['af_query']):>6.2f}  {s.get('flag', '')}"
    )
    if s.get("flag") != "BELOW_DETECTION":
        reported.append((n, het, unif, a_sk, a_fa))

print()
for label, idx in [("skani", 3), ("fastANI", 4)]:
    sub = [(r[1] - r[idx], r[2] - r[idx]) for r in reported if r[idx] is not None]
    if not sub:
        continue
    hb = sum(a for a, _ in sub) / len(sub)
    hm = sum(abs(a) for a, _ in sub) / len(sub)
    ub = sum(b for _, b in sub) / len(sub)
    um = sum(abs(b) for _, b in sub) / len(sub)
    print(
        f"reported pairs only, vs {label:<8} n={len(sub):<3} "
        f"het: bias={hb:+.3f} MAE={hm:.3f}   uniform: bias={ub:+.3f} MAE={um:.3f}"
    )
