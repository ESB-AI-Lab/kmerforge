# kmertools

Command-line tools for canonical k-mer counting, pairwise sample comparison, differential k-mer extraction, and read filtering on FASTQ/FASTA input (gzipped or plain).

Useful for cultivar identification, contamination screening, variant-enriched read extraction, and other k-mer-level comparisons between genomes or sequencing runs.

## Features

- **Canonical k-mers** - stores the lexicographically smaller of forward and reverse-complement, so strand is irrelevant
- **Gzip-transparent** - reads `.gz` files directly, no decompression step needed
- **Auto-detect format** - distinguishes FASTA (`>`) from FASTQ (`@`) automatically
- **Memory-efficient counting** - prefix-bucketed flat hash maps use ~9 bytes/slot vs ~48 bytes for `std::unordered_map` (~39 GB vs ~144 GB for a human genome at k=21)
- **Low-memory mode** - `--lowmem` disk-backed counting uses ~260 MB RAM regardless of genome size
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

For large genomes, use `--lowmem` to cap RAM at ~260 MB (uses temp files on disk):

```bash
./kmer_count -k 21 -m 2 --lowmem -o genome.kcounts reads.fastq.gz
./kmer_count -k 21 -m 2 --lowmem --tmpdir /scratch -o genome.kcounts reads.fastq.gz
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
| `kmer_count` | Count canonical k-mers from FASTQ/FASTA (plain or gzipped). Output as `.kcounts` binary or TSV. Supports `--lowmem` disk-backed mode for large genomes. |
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

## Algorithm

### K-mer encoding

Each nucleotide is encoded in 2 bits (A=0, C=1, G=2, T=3), so a k-mer of length k is packed into a single `uint64_t` (max k=31, using 62 bits). For each k-mer, the canonical form is computed as the lexicographic minimum of the forward and reverse-complement encodings, making counting strand-agnostic.

### Prefix-bucketed counting (k >= 9)

Rather than storing all k-mers in one large hash map, `kmer_count` decomposes each canonical k-mer into two parts:

```
k-mer (2k bits):  [  prefix (16 bits)  |  suffix (2k - 16 bits)  ]
                   top 8 bases           remaining bases
```

The 8-base prefix (16 bits) indexes into 65,536 independent flat hash map buckets. Each bucket stores only the suffix portion, which is smaller:

- **k <= 24**: suffix fits in `uint32_t` (up to 32 bits)
- **k 25-31**: suffix uses `uint64_t` (34-46 bits)

This decomposition has two benefits:

1. **Reduced key size** — for k <= 24, each key is 4 bytes instead of 8, saving ~30% of hash map memory
2. **Free global sort** — iterating buckets 0 through 65,535 and sorting within each bucket produces globally sorted output, which is required by the `.kcounts` format and the downstream merge-join algorithms in `kmer_diff` and `kmer_cosinesim`

#### FlatHashMap

Each bucket uses a struct-of-arrays (SoA) open-addressing hash map with linear probing:

| Array | Type | Description |
|-------|------|-------------|
| `keys[]` | `uint32_t` or `uint64_t` | Suffix values |
| `vals[]` | `uint32_t` | Counts |
| `state[]` | `uint8_t` | 0 = empty, 1 = occupied |

At 70% load factor, this costs ~9 bytes per slot (vs ~48 bytes per entry in `std::unordered_map`). Buckets use lazy initialization — empty buckets allocate nothing, so small inputs use minimal memory regardless of k.

**Estimated RAM for a human genome (~3 billion distinct 21-mers):**

| Method | Per entry | Total |
|--------|-----------|-------|
| `std::unordered_map` | ~48 bytes | ~144 GB |
| FlatHashMap (uint32_t suffix) | ~12.9 bytes | ~39 GB |

### Disk-backed mode (`--lowmem`)

When `--lowmem` is passed, suffixes are buffered in a flat array (1,024 entries per bucket, ~256 MB total for k <= 24) and flushed to per-bucket temp files on disk as the buffer fills. After scanning, each bucket file is loaded one at a time, sorted, counted by run-length, and written to the output. Only one bucket is in memory at a time during this phase.

```
Scan phase:     k-mer → prefix + suffix → buffer[prefix] → flush to disk when full
Count phase:    for each bucket file: load → sort → run-length count → write output → delete
```

Peak RAM is ~260 MB regardless of genome size. The tradeoff is disk I/O: temp files total ~4 bytes per k-mer occurrence (e.g., ~360 GB for 30x human coverage). Use `--tmpdir` to direct temp files to fast local storage.

### Small k (k <= 8)

For k <= 8, the entire k-mer space fits in at most 65,536 entries (4^8), so a simple `std::unordered_map` is used without prefix decomposition. The `--lowmem` flag is accepted but ignored with a note.

### Downstream tools

`kmer_diff` and `kmer_cosinesim` operate on the sorted `.kcounts` files using a merge-join (two-pointer) algorithm: both files are read in sorted order, matching k-mers are compared, and unmatched entries advance the smaller pointer. This runs in O(n+m) time with O(1) extra memory — no hash tables are built.

`kmer_reads` loads the differential k-mer set into a hash set, then streams the input FASTQ and emits reads containing (or lacking, with `--invert`) at least one matching k-mer.

## License

Apache 2.0 -- see [LICENSE](LICENSE).
