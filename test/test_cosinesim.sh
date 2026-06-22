#!/usr/bin/env bash
set -euo pipefail
PROJECT="$1"; TMP="$2"
DATA="$PROJECT/test"

"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/a.kcounts" "$DATA/sample_a.fastq" 2>/dev/null
"$PROJECT/kmer_count" -k 11 -m 1 -o "$TMP/b.kcounts" "$DATA/sample_b.fastq" 2>/dev/null

output=$("$PROJECT/kmer_cosinesim" "$TMP/a.kcounts" "$TMP/b.kcounts" 2>&1)

# Extract cosine value from the pairwise line
cosine=$(echo "$output" | grep -E '^\s*a\s+b\s+' | awk '{print $3}')

# Verify cosine is approximately 0.19 (allow 0.18-0.20)
awk -v c="$cosine" 'BEGIN { if (c+0 < 0.18 || c+0 > 0.20) { printf "cosine=%s, expected ~0.19\n", c; exit 1 } }'

# Verify shared=1 from the pairwise output
shared=$(echo "$output" | grep -E '^\s*a\s+b\s+' | awk '{print $6}')
[ "$shared" -eq 1 ] || { echo "shared=$shared, expected 1"; exit 1; }
