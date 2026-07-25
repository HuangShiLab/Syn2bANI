#!/usr/bin/env bash
# Three-way real-genome benchmark: syn2bani ani vs skani vs FastANI.
#
# Downloads 14 complete Enterobacteriaceae chromosomes from ENA spanning roughly
# 80-100% ANI against E. coli K-12 MG1655, then runs all three tools on the same
# pairs. Total download ~70 MB.
#
# There is no ground truth for real genome pairs, so this measures agreement,
# not accuracy — and FastANI's own reliability below ~92% ANI is contested. Read
# it as: does the estimate track the established tools where they are trusted,
# and does it decline to answer where they decline.
#
# Requires: skani, fastANI, python3.
set -euo pipefail

OUT="${1:-realbench}"
mkdir -p "$OUT/genomes"
cd "$OUT"

cat > accessions.tsv <<'EOF'
U00096.3	Ecoli_K12_MG1655
AP009048.1	Ecoli_K12_W3110
AM946981.2	Ecoli_BL21_DE3
BA000007.3	Ecoli_O157H7_Sakai
AE014075.1	Ecoli_CFT073
CP000243.1	Ecoli_UTI89
AE005674.2	Shigella_flexneri_301
CP000038.1	Shigella_sonnei_Ss046
CU928158.2	Escherichia_fergusonii
AE006468.2	Salmonella_Typhimurium_LT2
AL513382.1	Salmonella_Typhi_CT18
FN543502.1	Citrobacter_rodentium
CP000647.1	Klebsiella_pneumoniae_MGH78578
CP001918.1	Enterobacter_cloacae_13047
EOF

while IFS=$'\t' read -r acc name; do
    out="genomes/${name}.fasta"
    [ -s "$out" ] && continue
    curl -sSL "https://www.ebi.ac.uk/ena/browser/api/fasta/${acc}?download=true" -o "$out"
    b=$(grep -v '^>' "$out" | tr -d '\n\r ' | wc -c | tr -d ' ')
    printf "  %-32s %-12s %10s bp\n" "$name" "$acc" "$b"
done < accessions.tsv

REF=genomes/Ecoli_K12_MG1655.fasta
ls genomes/*.fasta | grep -v Ecoli_K12_MG1655 > qlist.txt
echo "$REF" > rlist.txt

BIN="${SYN2BANI:-../../target/release/syn2bani}"
echo "== syn2bani ani =="
"$BIN" ani $(cat qlist.txt) "$REF" --verbose -p -o syn2bani.tsv
echo "== skani =="
skani dist --ql qlist.txt --rl rlist.txt -o skani.tsv -t 8
echo "== fastANI =="
fastANI --ql qlist.txt -r "$REF" -o fastani.tsv -t 8 2>/dev/null

cp ../realgenome_compare.py compare.py
python3 compare.py
