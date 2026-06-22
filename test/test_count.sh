#!/usr/bin/env bash
set -euo pipefail
PROJECT="$1"; TMP="$2"
DATA="$PROJECT/test"

"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/a.kcounts" "$DATA/sample_a.fastq" 2>/dev/null
"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/b.kcounts" "$DATA/sample_b.fastq" 2>/dev/null

# Verify expected k-mer counts from binary header (magic + k + count)
a_count=$(od -A n -t u8 -j 8 -N 8 "$TMP/a.kcounts" | tr -d ' ')
b_count=$(od -A n -t u8 -j 8 -N 8 "$TMP/b.kcounts" | tr -d ' ')

[ "$a_count" -eq 12 ] || { echo "A count=$a_count, expected 12"; exit 1; }
[ "$b_count" -eq 20 ] || { echo "B count=$b_count, expected 20"; exit 1; }
