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
REPS=3

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

# Best-of-REPS wall time and the peak RSS seen. Arithmetic stays in the shell:
# nesting quotes into awk from inside a command substitution is how the format
# string gets eaten.
measure() {
    label="$1"; shift
    best_t=""; best_kb=""
    "$@" >/dev/null 2>&1 || true                     # warm-up, result discarded
    i=0
    while [ "$i" -lt "$REPS" ]; do
        i=$((i + 1))
        out=$(/usr/bin/time -l "$@" 2>&1 >/dev/null || true)
        t=$(printf '%s\n' "$out" | awk '/ real /{print $1; exit}')
        b=$(printf '%s\n' "$out" | awk '/maximum resident set size/{print $1; exit}')
        [ -z "$t" ] && continue
        # Compare times as integer milliseconds to stay out of awk.
        tms=$(printf '%s' "$t" | awk '{printf "%d", $1 * 1000}')
        if [ -z "$best_t" ] || [ "$tms" -lt "$best_tms" ]; then
            best_t="$t"; best_tms="$tms"
        fi
        if [ -n "$b" ] && { [ -z "$best_kb" ] || [ "$b" -gt "$best_kb" ]; }; then
            best_kb="$b"
        fi
    done
    mb=$(( ${best_kb:-0} / 1048576 ))
    printf "  %-26s %8s s %7s MB\n" "$label" "${best_t:-n/a}" "$mb"
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

    echo "-- W3: 14x14 all-vs-all (196 pairs) --"
    measure "syn2bani ani"  "$BIN" ani $(cat perf_all.txt) $(cat perf_all.txt) "${SYN_PAR[@]}"
    measure "skani dist"    skani dist --ql perf_all.txt --rl perf_all.txt -t "$T"
    measure "skani triangle" skani triangle -l perf_all.txt -t "$T"

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
echo "skani can sketch once and reuse; syn2bani ani re-digests every run."
rm -rf perf_sketch && mkdir -p perf_sketch
measure "skani sketch (14 genomes)" skani sketch -l perf_all.txt -o perf_sketch -t 8
ls perf_sketch/*.sketch >/dev/null 2>&1 && du -sh perf_sketch | sed 's/^/  sketch dir: /'
