#!/usr/bin/env bash
# Wall time and peak RSS for `syn2bani ani` against skani and FastANI.
#
# Fairness notes, because the tools are not shaped the same way:
#
# - Thread counts are pinned identically. `syn2bani` needs `-p` to use more than
#   one thread at all, so `-t N -p` is the comparable form.
# - Each workload runs as ONE invocation of each tool, so per-genome setup
#   (digestion for syn2bani, sketching for skani) is amortised the same way.
# - skani is additionally measured in its two-phase form (`sketch` once, then
#   `dist` on the sketches). `syn2bani ani` has no equivalent reuse, so that row
#   is the honest picture of an all-vs-all or database workload, not a like-for-
#   like comparison.
# - Every measurement is the best of three, after a warm-up run, so the numbers
#   reflect warm page cache rather than first-read disk latency.
#
# Usage: bash perf_bench.sh [bench_dir]
# The bench dir needs genomes/ and optionally drafts/ as built by
# realgenome_bench.sh and draft_bench.sh.
set -euo pipefail

DIR="${1:-.}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${SYN2BANI:-$REPO/target/release/syn2bani}"
REPS=5

if [ ! -x "$BIN" ]; then
    echo "error: no binary at $BIN — run: (cd $REPO && cargo build --release)" >&2
    exit 1
fi
newer=$(find "$REPO/src" "$REPO/Cargo.toml" -newer "$BIN" -print -quit 2>/dev/null || true)
if [ -n "$newer" ]; then
    echo "error: $BIN is older than the source tree ($newer)." >&2
    echo "       run: (cd $REPO && cargo build --release)" >&2
    exit 1
fi

cd "$DIR"

# min/median/max wall time over REPS runs, plus the largest peak RSS seen.
#
# Report the spread, not a single number. On a machine doing anything else at
# all, every FASTA-reading workload here — including skani's — varies by more
# than an order of magnitude between its best and worst run, so a point estimate
# invites ratio claims the data cannot support. Only the sketch-input rows are
# reproducible to a few percent.
measure() {
    label="$1"; shift
    "$@" >/dev/null 2>&1 || true                     # warm-up, result discarded
    unset best_kb
    times=""
    i=0
    while [ "$i" -lt "$REPS" ]; do
        i=$((i + 1))
        out=$(/usr/bin/time -l "$@" 2>&1 >/dev/null || true)
        t=$(printf '%s\n' "$out" | awk '/ real /{print $1; exit}')
        b=$(printf '%s\n' "$out" | awk '/maximum resident set size/{print $1; exit}')
        [ -z "$t" ] && continue
        times="$times$t
"
        if [ -n "$b" ] && { [ -z "${best_kb:-}" ] || [ "$b" -gt "$best_kb" ]; }; then
            best_kb="$b"
        fi
    done
    sorted=$(printf '%s' "$times" | sort -n)
    n=$(printf '%s\n' "$sorted" | grep -c .)
    lo=$(printf '%s\n' "$sorted" | sed -n 1p)
    md=$(printf '%s\n' "$sorted" | sed -n "$(( (n + 1) / 2 ))p")
    hi=$(printf '%s\n' "$sorted" | sed -n "${n}p")
    printf "  %-24s %7s / %7s / %7s s  %6s MB\n" \
        "$label" "${lo:-n/a}" "${md:-n/a}" "${hi:-n/a}" "$(( ${best_kb:-0} / 1048576 ))"
}

REF=genomes/Ecoli_K12_MG1655.fasta
ls genomes/*.fasta | grep -v Ecoli_K12_MG1655 > perf_q13.txt
ls genomes/*.fasta > perf_all.txt
echo "$REF" > perf_ref.txt

for T in 1 8; do
    # Never leave this array empty: macOS ships bash 3.2, where "\${a[@]}" on an
    # empty array trips `set -u`. Without -p syn2bani stays single-threaded.
    if [ "$T" = 1 ]; then SYN_PAR=(-t 1); else SYN_PAR=(-t "$T" -p); fi
    echo
    echo "══ $T thread(s) ══"

    echo "-- W1: 1 pair, complete genomes --"
    measure "syn2bani ani"  "$BIN" ani genomes/Ecoli_O157H7_Sakai.fasta "$REF" "${SYN_PAR[@]}"
    measure "skani dist"    skani dist -q genomes/Ecoli_O157H7_Sakai.fasta -r "$REF" -t "$T"
    measure "fastANI"       fastANI -q genomes/Ecoli_O157H7_Sakai.fasta -r "$REF" -o /dev/null -t "$T"

    echo "-- W2: 13 pairs vs one reference --"
    measure "syn2bani ani"  "$BIN" ani $(cat perf_q13.txt) "$REF" "${SYN_PAR[@]}"
    measure "skani dist"    skani dist --ql perf_q13.txt --rl perf_ref.txt -t "$T"
    measure "fastANI"       fastANI --ql perf_q13.txt -r "$REF" -o /dev/null -t "$T"

    # --ql/--rl, not two positional lists: two greedy positional Vecs cannot be
    # split, so `ani a b c ... z` is (n-1) queries against one reference, not
    # all-vs-all. Measuring it that way understated the work by 7x here.
    echo "-- W3: 14x14 all-vs-all (196 pairs) --"
    measure "syn2bani FASTA"  "$BIN" ani --ql perf_all.txt --rl perf_all.txt "${SYN_PAR[@]}"
    measure "skani dist"      skani dist --ql perf_all.txt --rl perf_all.txt -t "$T"
    measure "skani triangle"  skani triangle -l perf_all.txt -t "$T"

    if [ -d drafts ]; then
        ls drafts/GCA_*.fasta | grep -v '\.rc\.' > perf_drafts.txt
        echo "-- W4: 8 real drafts vs one reference (incl. 8025-contig MAG) --"
        measure "syn2bani ani"  "$BIN" ani $(cat perf_drafts.txt) "$REF" "${SYN_PAR[@]}"
        measure "skani dist"    skani dist --ql perf_drafts.txt --rl perf_ref.txt -t "$T"
        measure "fastANI"       fastANI --ql perf_drafts.txt -r "$REF" -o /dev/null -t "$T"
    fi
done

echo
echo "══ sketch reuse (8 threads) ══"
PANEL="${PANEL:-BcgI,AlfI,AloI,FalI}"
rm -rf perf_s2ba perf_sketch && mkdir -p perf_s2ba perf_sketch
measure "syn2bani sketch (14)"  "$BIN" sketch $(cat perf_all.txt) -o perf_s2ba --enzymes "$PANEL" -t 8 -p
measure "skani sketch (14)"     skani sketch -l perf_all.txt -o perf_sketch -t 8
ls perf_s2ba/*.s2ba > perf_s2ba.txt 2>/dev/null || true
if [ -s perf_s2ba.txt ]; then
    du -sh perf_s2ba | sed 's/^/  syn2bani .s2ba dir: /'
    echo "-- 196 pairs from sketches --"
    measure "syn2bani sketches" "$BIN" ani --ql perf_s2ba.txt --rl perf_s2ba.txt -t 8 -p
fi
du -sh perf_sketch 2>/dev/null | sed 's/^/  skani sketch dir:   /'
