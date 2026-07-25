#!/usr/bin/env bash
# Real draft-assembly benchmark.
#
# Simulated partitioning (fragment.py) keeps every base and puts contig
# boundaries at random. Real drafts also lose sequence at repeat boundaries,
# carry contamination, and have contig lengths set by coverage and repeat
# structure. This downloads eight real contig-level E. coli assemblies from ENA,
# spanning 88 to 8025 contigs, and runs two tests:
#
# 1. **Reverse-complement self-control.** Every contig of a draft is reverse
#    complemented. Sequence content is unchanged, so ANI must be 100%. This is
#    the test that validates strand-canonical tag hashing on real assemblies —
#    real contigs come in arbitrary orientation, and without canonicalization the
#    tags in a reverse-oriented contig cannot match their homologs at all.
#
# 2. **Three-way comparison** against E. coli K-12 MG1655, with skani and
#    FastANI, to see how the estimate degrades with real fragmentation.
#
# Requires: skani, fastANI, python3. Downloads ~45 MB.
set -euo pipefail

OUT="${1:-draftbench}"
BIN="${SYN2BANI:-$(cd "$(dirname "$0")/.." && pwd)/target/release/syn2bani}"
mkdir -p "$OUT/drafts"
cd "$OUT"

# A spread of real fragmentation levels. The last one is 8025 contigs and
# 8.9 Mb — roughly twice an E. coli genome, so it is a contaminated or mixed
# assembly, i.e. a realistically bad MAG.
ACCESSIONS="GCA_001283865 GCA_001077875 GCA_001284645 GCA_001283245 GCA_001283605 GCA_001284145 GCA_001283205 GCA_001075925"

for g in $ACCESSIONS; do
    [ -s "drafts/$g.fasta" ] && continue
    curl -sSL "https://www.ebi.ac.uk/ena/browser/api/fasta/${g}?download=true" -o "drafts/$g.fasta"
done

[ -s ref.fasta ] || curl -sSL \
    "https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.3?download=true" -o ref.fasta

python3 - "$ACCESSIONS" <<'PY'
import sys, os
COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

def contigs(path):
    recs, name, cur = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(cur)))
            name, cur = line.strip(), []
        else:
            cur.append(line.strip())
    if name is not None:
        recs.append((name, "".join(cur)))
    return recs

def n50(lens):
    lens = sorted(lens, reverse=True)
    half, acc = sum(lens) / 2, 0
    for l in lens:
        acc += l
        if acc >= half:
            return l
    return 0

print(f"{'assembly':<16}{'contigs':>8}{'total bp':>12}{'N50':>10}")
for g in sys.argv[1].split():
    recs = contigs(f"drafts/{g}.fasta")
    lens = [len(s) for _, s in recs]
    print(f"{g:<16}{len(recs):>8}{sum(lens):>12,}{n50(lens):>10,}")
    out = f"drafts/{g}.rc.fasta"
    if not os.path.exists(out):
        with open(out, "w") as fh:
            for n, s in recs:
                fh.write(n + "_rc\n")
                d = s.encode().translate(COMP)[::-1].decode()
                for k in range(0, len(d), 80):
                    fh.write(d[k : k + 80] + "\n")
PY

echo
echo "== reverse-complement self-control (must be 100.00) =="
printf "%-16s %12s %8s %s\n" assembly ani AF flag
for g in $ACCESSIONS; do
    r=$("$BIN" ani "drafts/$g.rc.fasta" "drafts/$g.fasta" --verbose -p 2>/dev/null | tail -1)
    printf "%-16s %12s %8s %s\n" "$g" "$(echo "$r" | cut -f3)" \
        "$(echo "$r" | cut -f5)" "$(echo "$r" | cut -f15)"
done

echo
echo "== three-way vs E. coli K-12 MG1655 =="
ls drafts/GCA_*.fasta | grep -v '\.rc\.' > qlist.txt
echo ref.fasta > rlist.txt
"$BIN" ani $(cat qlist.txt) ref.fasta --verbose -p -o syn.tsv
skani dist --ql qlist.txt --rl rlist.txt -o skani.tsv -t 8 2>/dev/null
fastANI --ql qlist.txt -r ref.fasta -o fastani.tsv -t 8 2>/dev/null
echo "wrote syn.tsv skani.tsv fastani.tsv in $OUT"
