#!/usr/bin/env python3
"""Draw a stratified pair subset for independent (ANIm) truth.

Panel selection needs truth on a few hundred to a few thousand pairs, not on all
45,000 — panels are compared on the same pairs, so the comparison is paired and
most between-pair variance cancels. What it does need is *coverage*: the reported
bias varies by phylum (2.4–4.6%) and by ANI band, so a panel chosen on a pooled
sample optimises an average that may serve no regime well.

This samples evenly across ANI band x group cells, reports every under-filled
cell rather than quietly returning fewer pairs, and writes a pair list plus an
optional SLURM array driver for nucmer/dnadiff.

Two deliberate choices worth knowing about:

* **Banding uses an existing estimate, and that is fine.** The band only decides
  which pairs get truth computed, so it affects coverage, not the estimate being
  judged. Prefer a column produced by a method that is *not* the one under
  selection — skani rather than Syn2bANI — so that pairs are not chosen by the
  thing being evaluated. Never band on the truth column itself.
* **Symmetric pairs are collapsed.** ANIm is near-symmetric, so computing both
  (A,B) and (B,A) spends half the budget twice.

Usage:
    python3 stratified_sample.py PAIRS.tsv --out sample.tsv \\
        [--ani-col skani_ani] [--group-col phylum] \\
        [--per-cell 60] [--bands 80,85,90,95,99,100] \\
        [--genome-dir DIR] [--slurm anim_array.sh]
"""
import argparse
import os
import random
import sys
from collections import defaultdict


def sniff(header, candidates, what):
    """Pick the first header field containing any candidate substring, case-insensitively."""
    low = {h.lower(): h for h in header}
    for c in candidates:
        needle = c.lower()
        for h_lower, h in low.items():
            if needle in h_lower:
                return h
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pairs", help="TSV of candidate pairs, with a header")
    ap.add_argument("--out", required=True, help="sampled pair list (TSV)")
    ap.add_argument("--query-col", default=None)
    ap.add_argument("--ref-col", default=None)
    ap.add_argument("--ani-col", default=None,
                    help="column used only to assign ANI bands (default: prefer skani)")
    ap.add_argument("--group-col", default=None,
                    help="second stratification key, e.g. phylum (default: none)")
    ap.add_argument("--bands", default="80,85,90,95,99,100",
                    help="band edges in percent")
    ap.add_argument("--per-cell", type=int, default=60,
                    help="pairs to draw per band x group cell")
    ap.add_argument("--seed", type=int, default=20260728)
    ap.add_argument("--genome-dir", default=None,
                    help="prefix for genome filenames when writing paths")
    ap.add_argument("--suffix", default=".fna",
                    help="genome filename suffix under --genome-dir")
    ap.add_argument("--slurm", default=None,
                    help="also write a SLURM array driver running dnadiff")
    ap.add_argument("--chunk", type=int, default=50,
                    help="pairs per SLURM array task")
    args = ap.parse_args()

    with open(args.pairs) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    idx = {h: i for i, h in enumerate(header)}

    qc = args.query_col or sniff(header, ["query", "query_id", "genome_a", "q"], "query")
    rc = args.ref_col or sniff(header, ["reference", "ref", "ref_id", "genome_b", "r"], "reference")
    # Prefer an estimate from a method that is not the one being selected.
    ac = args.ani_col or sniff(
        header, ["skani_ani", "skani", "fastani_ani", "fastani", "ani_uniform", "ani"], "ani")
    gc = args.group_col or sniff(header, ["phylum", "p__", "group", "clade"], "group")

    missing = [n for n, v in [("query", qc), ("reference", rc), ("ani", ac)] if v is None]
    if missing:
        print(f"error: could not find column(s) for {', '.join(missing)}", file=sys.stderr)
        print(f"       header is: {header}", file=sys.stderr)
        return 1
    print(f"columns: query={qc}  reference={rc}  band-by={ac}  group={gc or '(none)'}")

    edges = [float(x) for x in args.bands.split(",")]
    def band_of(v):
        for i in range(len(edges) - 1):
            if edges[i] <= v < edges[i + 1]:
                return f"{edges[i]:g}-{edges[i+1]:g}"
        return None

    # Collapse symmetric duplicates: ANIm is near-symmetric, so (A,B) and (B,A)
    # would spend the same budget twice.
    cells = defaultdict(list)
    seen = set()
    skipped_nan = skipped_band = dup = 0
    for r in rows:
        try:
            q, ref, v = r[idx[qc]], r[idx[rc]], float(r[idx[ac]])
        except (IndexError, ValueError):
            skipped_nan += 1
            continue
        if q == ref:
            continue
        key = tuple(sorted((q, ref)))
        if key in seen:
            dup += 1
            continue
        seen.add(key)
        b = band_of(v)
        if b is None:
            skipped_band += 1
            continue
        g = r[idx[gc]] if gc and idx.get(gc) is not None and len(r) > idx[gc] else "all"
        cells[(b, g)].append((q, ref, v, g))

    rng = random.Random(args.seed)
    sampled, short = [], []
    for cell in sorted(cells):
        pool = cells[cell]
        take = min(args.per_cell, len(pool))
        if take < args.per_cell:
            short.append((cell, len(pool)))
        sampled.extend(rng.sample(pool, take))

    print(f"\n{len(rows)} rows -> {len(seen)} unique unordered pairs "
          f"({dup} symmetric duplicates dropped)")
    if skipped_nan or skipped_band:
        print(f"  skipped: {skipped_nan} unparseable, {skipped_band} outside the bands")
    print(f"\n{'band':<12}{'group':<28}{'available':>10}{'sampled':>9}")
    for cell in sorted(cells):
        n_av = len(cells[cell])
        n_s = min(args.per_cell, n_av)
        flag = "  <- under-filled" if n_av < args.per_cell else ""
        print(f"{cell[0]:<12}{cell[1]:<28}{n_av:>10}{n_s:>9}{flag}")
    print(f"\ntotal sampled: {len(sampled)} pairs across {len(cells)} cells")
    if short:
        print(f"WARNING: {len(short)} cell(s) could not supply {args.per_cell} pairs. "
              f"A stratified design that silently returns fewer pairs in exactly the "
              f"regimes that are hardest to sample is worse than no stratification — "
              f"either lower --per-cell or widen those bands.")

    def path(name):
        return os.path.join(args.genome_dir, name + args.suffix) if args.genome_dir else name

    with open(args.out, "w") as fh:
        fh.write("query\treference\tband\tgroup\tband_by_ani\tquery_path\tref_path\n")
        for q, r, v, g in sampled:
            fh.write(f"{q}\t{r}\t{band_of(v)}\t{g}\t{v:.4f}\t{path(q)}\t{path(r)}\n")
    print(f"wrote {args.out}")

    if args.slurm:
        n_tasks = (len(sampled) + args.chunk - 1) // args.chunk
        with open(args.slurm, "w") as fh:
            fh.write(f"""#!/bin/bash
#SBATCH --job-name=anim
#SBATCH --array=1-{n_tasks}
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=anim_%A_%a.out
# dnadiff over the sampled pairs. ANIm is the independent reference: selecting a
# panel against FastANI would optimise agreement with FastANI, including its own
# bias below 95% ANI.
set -euo pipefail
SAMPLE="{os.path.abspath(args.out)}"
CHUNK={args.chunk}
OUT="${{OUT:-anim_results}}"
mkdir -p "$OUT"
start=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 2 ))   # +2 skips the header
end=$(( start + CHUNK - 1 ))
sed -n "${{start}},${{end}}p" "$SAMPLE" | while IFS=$'\\t' read -r q r band grp v qp rp; do
    [ -z "${{q:-}}" ] && continue
    tag="${{q}}__${{r}}"
    pre="$OUT/$tag"
    [ -s "$pre.report" ] && continue
    dnadiff -p "$pre" "$rp" "$qp" >/dev/null 2>&1 || {{ echo "FAILED $tag" >&2; continue; }}
    # AvgIdentity from the 1-to-1 block is the ANIm value.
    ident=$(awk '/^AvgIdentity/ {{print $2; exit}}' "$pre.report")
    aln=$(awk '/^AlignedBases/ {{print $2; exit}}' "$pre.report")
    printf '%s\\t%s\\t%s\\t%s\\n' "$q" "$r" "${{ident:-NA}}" "${{aln:-NA}}" \\
        >> "$OUT/anim_${{SLURM_ARRAY_TASK_ID}}.tsv"
    rm -f "$pre".{{delta,1delta,mdelta,1coords,mcoords,qdiff,rdiff,snps,unref,unqry}}
done
""")
        os.chmod(args.slurm, 0o755)
        print(f"wrote {args.slurm}  ({n_tasks} array tasks x {args.chunk} pairs)")
        print("  after it finishes:  cat anim_results/anim_*.tsv > truth.tsv")
        print("  then:               syn2bani panel --strata strata.tsv --truth truth.tsv --greedy")
    return 0


if __name__ == "__main__":
    sys.exit(main())
