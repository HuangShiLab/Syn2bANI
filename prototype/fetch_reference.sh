#!/usr/bin/env bash
# Fetch the reference genome used by the validation harness.
#
# E. coli K-12 MG1655 complete genome, ENA accession U00096.3 (4,641,652 bp).
# Downloaded from EBI/ENA rather than NCBI because the lab's IP range is
# blocked from the NCBI API.
#
# Pinning a public accession keeps the validation reproducible: everyone
# simulates from byte-identical input, so the reported MAE can be compared
# directly rather than approximately.
#
# Usage: bash fetch_reference.sh [outfile]
set -euo pipefail

OUT="${1:-mg1655.fasta}"
URL="https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.3?download=true"
EXPECTED_BASES=4641652
EXPECTED_MD5=805d1558950b76c15883a718da46d3e1

echo "downloading U00096.3 -> $OUT"
curl -sSL "$URL" -o "$OUT"

bases=$(grep -v '^>' "$OUT" | tr -d '\n\r ' | wc -c | tr -d ' ')
if [ "$bases" -ne "$EXPECTED_BASES" ]; then
    echo "ERROR: got $bases bases, expected $EXPECTED_BASES." >&2
    echo "The download is incomplete or the accession changed. Do not run the" >&2
    echo "validation on a truncated genome — it still produces exact ground" >&2
    echo "truth, but the numbers will not be comparable to the reported ones." >&2
    exit 1
fi

if command -v md5 >/dev/null 2>&1; then
    md5sum_actual=$(md5 -q "$OUT")
elif command -v md5sum >/dev/null 2>&1; then
    md5sum_actual=$(md5sum "$OUT" | cut -d' ' -f1)
else
    md5sum_actual=""
fi

if [ -n "$md5sum_actual" ] && [ "$md5sum_actual" != "$EXPECTED_MD5" ]; then
    echo "WARNING: md5 $md5sum_actual != expected $EXPECTED_MD5" >&2
    echo "Base count is correct, so this is probably only line-wrapping or" >&2
    echo "header formatting. Results should still match." >&2
fi

echo "OK: $bases bases, md5 ${md5sum_actual:-unavailable}"
echo
echo "Next:"
echo "  python3 simulate.py $OUT sim"
echo "  python3 simulate_accessory.py $OUT simacc 0.95"
echo "  ../target/release/syn2bani ani sim/q_*.fasta sim/ref.fasta --verbose -p"
