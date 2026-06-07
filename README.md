# kmertools

Command-line tools for canonical k-mer counting, pairwise sample comparison, differential k-mer extraction, and read filtering on FASTQ/FASTA input (gzipped or plain).

Useful for cultivar identification, contamination screening, variant-enriched read extraction, and other k-mer-level comparisons between genomes or sequencing runs.

## Features

- **Canonical k-mers** - stores the lexicographically smaller of forward and reverse-complement, so strand is irrelevant
- **Gzip-transparent** - reads `.gz` files directly, no decompression step needed
- **Auto-detect format** - distinguishes FASTA (`>`) from FASTQ (`@`) automatically
- **Streaming merge-join** - `kmer_diff` and `kmer_cosinesim` run in O(n+m) time with O(1) extra memory
- **Binary `.kcounts` format** - sorted (hash, count) pairs for fast downstream operations

## Prerequisites

- C++17 compiler (GCC 7+, Clang 5+)
- zlib development headers (`zlib.h`)

On Debian/Ubuntu:
```bash
sudo apt install build-essential zlib1g-dev
```

On macOS:
```bash
xcode-select --install
# zlib is included with the Command Line Tools
```

## Installation

### Make (recommended)

```bash
git clone https://github.com/esolares/kmertools.git
cd kmertools
make -j$(nproc)
```

Binaries are built in the project root: `kmer_count`, `kmer_diff`, `kmer_cosinesim`, `kmer_reads`.

### CMake

```bash
git clone https://github.com/esolares/kmertools.git
cd kmertools
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Binaries are built in `build/`.

### Docker

```bash
docker build -t kmertools .
docker run --rm -v $(pwd):/data kmertools kmer_count -k 11 -m 1 -o /data/out.kcounts /data/input.fastq
```

The image uses a multi-stage build -- the final image contains only the runtime dependency (zlib) and the four binaries in `/usr/local/bin/`.

### Singularity / Apptainer

```bash
singularity build kmertools.sif singularity.def
singularity exec kmertools.sif kmer_count -k 11 -m 1 -o out.kcounts input.fastq
```

Singularity containers bind the current directory by default, so input/output files are accessible without extra flags.

## Quick start

Test data is included in `test/`. The examples below use it.

### 1. Count k-mers

Count canonical 11-mers in each sample, keeping everything (minimum count 1):

```bash
./kmer_count -k 11 -m 1 -o sample_a.kcounts test/sample_a.fastq
./kmer_count -k 11 -m 1 -o sample_b.kcounts test/sample_b.fastq
```

Output:
```
Counting 11-mers from test/sample_a.fastq
  6 reads processed, 12 unique k-mers
After min_count=1 filter: 12 k-mers
Wrote 12 k-mers to sample_a.kcounts (binary)
```

Use `-f tsv` for human-readable output:

```bash
./kmer_count -k 11 -m 1 -f tsv test/sample_a.fastq | head -5
```

### 2. Compare samples

Pairwise cosine similarity, Jaccard index, and angular distance:

```bash
./kmer_cosinesim sample_a.kcounts sample_b.kcounts
```

Output:
```
=== PAIRWISE SIMILARITY ===

Sample_A              Sample_B                  Cosine     Jaccard     Angular        Shared         Union
--------              --------                  ------     -------     -------        ------         -----
sample_a              sample_b                0.192726    0.032258     78.8881             1            31

=== INTERPRETATION ===

  Cosine ~1.000   : Clonal / identical
  Cosine 0.99+    : Very closely related (sport / bud mutation)
  Cosine 0.95-0.99: Close relatives (seedling, same race)
  Cosine < 0.95   : Different cultivars / races
```

Accepts any number of `.kcounts` files. Three or more produces a full pairwise matrix.

### 3. Find differential k-mers

K-mers present in sample B but not in sample A:

```bash
./kmer_diff -a sample_a.kcounts -b sample_b.kcounts -o diff.kcounts --only-b
```

Output:
```
=== K-MER DIFF SUMMARY ===
  A: sample_a.kcounts (12 k-mers)
  B: sample_b.kcounts (20 k-mers)
  Shared (same count) : 1
  Shared (diff count) : 0 (ratio >= 3.0)
  Only in A (lost)    : 11
  Only in B (gained)  : 19
  Total differential  : 30
  Written to          : diff.kcounts (19 entries)
```

Other modes: `--only-a` (lost k-mers), `--diff-ab` (shared but different abundance), `--all` (default).

Use `-t diff.tsv --seq` to write a TSV with decoded nucleotide sequences.

### 4. Extract differential reads

FASTQ reads from sample B that contain at least one differential k-mer:

```bash
./kmer_reads -k 11 -d diff.kcounts -i test/sample_b.fastq -o unique_reads.fastq.gz
```

Output:
```
=== READ EXTRACTION SUMMARY ===
  Input reads   : 6
  Matched reads : 5 (83.33%)
  Output reads  : 5 (differential reads)
  Output file   : unique_reads.fastq.gz
```

`--invert` extracts reads that do not match (shared reads):

```bash
./kmer_reads -k 11 -d diff.kcounts -i test/sample_b.fastq -o shared_reads.fastq.gz --invert
```

## Tool reference

| Tool | Description |
|------|-------------|
| `kmer_count` | Count canonical k-mers from FASTQ/FASTA (plain or gzipped). Output as `.kcounts` binary or TSV. |
| `kmer_diff` | Differential k-mers between two `.kcounts` files. Streaming merge-join. |
| `kmer_cosinesim` | Pairwise cosine similarity, Jaccard index, and angular distance across multiple samples. |
| `kmer_reads` | Extract (or exclude) FASTQ reads matching a k-mer set from `kmer_diff`. |

Run any tool with `-h` for full options.

## .kcounts format

Sorted binary k-mer count data:

| Field | Size | Description |
|-------|------|-------------|
| Magic | 4 bytes | `0x4B4D4552` ("KMER") |
| k | 4 bytes | k-mer size |
| n | 8 bytes | number of entries |
| entries | 12 bytes each | `uint64_t` k-mer hash + `uint32_t` count, sorted by hash |

K-mers use 2-bit encoding (A=0, C=1, G=2, T=3), canonical form (min of forward and reverse complement). Max k = 31. Sorted order is required -- `kmer_diff` and `kmer_cosinesim` depend on it for merge-join.

## License

Apache 2.0 -- see [LICENSE](LICENSE).
