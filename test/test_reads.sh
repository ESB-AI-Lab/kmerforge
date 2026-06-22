#!/usr/bin/env bash
set -euo pipefail
PROJECT="$1"; TMP="$2"
DATA="$PROJECT/test"

"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/a.kcounts" "$DATA/sample_a.fastq" 2>/dev/null
"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/b.kcounts" "$DATA/sample_b.fastq" 2>/dev/null
"$PROJECT/kmer_diff" -a "$TMP/a.kcounts" -b "$TMP/b.kcounts" -o "$TMP/diff.kcounts" --only-b 2>/dev/null

output=$("$PROJECT/kmer_reads" -k 11 -d "$TMP/diff.kcounts" -i "$DATA/sample_b.fastq" -o "$TMP/unique.fastq.gz" 2>&1)

matched=$(echo "$output" | grep "Matched reads" | grep -oP ':\s*\K\d+')
[ "$matched" -eq 5 ] || { echo "matched=$matched, expected 5"; exit 1; }

# Decompress and count reads in output
read_count=$(zcat "$TMP/unique.fastq.gz" | awk 'NR%4==1' | wc -l)
[ "$read_count" -eq 5 ] || { echo "output has $read_count reads, expected 5"; exit 1; }
