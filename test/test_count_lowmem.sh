#!/usr/bin/env bash
set -euo pipefail
PROJECT="$1"; TMP="$2"
DATA="$PROJECT/test"

# Count with both modes
"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/a.kcounts" "$DATA/sample_a.fastq" 2>/dev/null
"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/b.kcounts" "$DATA/sample_b.fastq" 2>/dev/null
"$PROJECT/kmer_count" -k 11 -m 1 --lowmem -o "$TMP/a_low.kcounts" "$DATA/sample_a.fastq" 2>/dev/null
"$PROJECT/kmer_count" -k 11 -m 1 --lowmem -o "$TMP/b_low.kcounts" "$DATA/sample_b.fastq" 2>/dev/null

cmp "$TMP/a.kcounts" "$TMP/a_low.kcounts" || { echo "A: mem vs lowmem mismatch"; exit 1; }
cmp "$TMP/b.kcounts" "$TMP/b_low.kcounts" || { echo "B: mem vs lowmem mismatch"; exit 1; }
