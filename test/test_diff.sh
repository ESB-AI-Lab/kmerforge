#!/usr/bin/env bash
set -euo pipefail
PROJECT="$1"; TMP="$2"
DATA="$PROJECT/test"

"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/a.kcounts" "$DATA/sample_a.fastq" 2>/dev/null
"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/b.kcounts" "$DATA/sample_b.fastq" 2>/dev/null
output=$("$PROJECT/kmer_diff" -a "$TMP/a.kcounts" -b "$TMP/b.kcounts" -o "$TMP/diff.kcounts" --only-b 2>&1)

gained=$(echo "$output" | grep "Only in B" | grep -oP '\d+')
shared=$(echo "$output" | grep "Shared (same count)" | grep -oP '\d+')

[ "$gained" -eq 19 ] || { echo "gained=$gained, expected 19"; exit 1; }
[ "$shared" -eq 1 ] || { echo "shared=$shared, expected 1"; exit 1; }

# Verify the output file has 19 entries
diff_count=$(od -A n -t u8 -j 8 -N 8 "$TMP/diff.kcounts" | tr -d ' ')
[ "$diff_count" -eq 19 ] || { echo "diff file has $diff_count entries, expected 19"; exit 1; }
