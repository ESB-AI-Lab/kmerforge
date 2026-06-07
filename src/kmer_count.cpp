/*
 * kmer_count.cpp — Fast canonical k-mer counter for FASTQ/FASTA (gzipped or plain)
 *
 * MIT License — Edwin Solares / ESB AI Lab
 *
 * Usage:
 *   ./kmer_count -k 21 -m 2 -o output.kcounts input.fastq.gz
 *
 * Options:
 *   -k INT      k-mer size (default: 21, max: 31)
 *   -t INT      threads for reading (default: 1)
 *   -m INT      minimum count to output (default: 2)
 *   -f tsv      output format: 'tsv' or 'bin' (default: bin)
 *   -o FILE     output file (default: stdout for tsv, required for bin)
 *   --lowmem    disk-backed mode (~260 MB RAM, uses temp files)
 *   --tmpdir DIR temp directory for --lowmem (default: $TMPDIR or /tmp)
 *
 * Compile:
 *   g++ -O3 -std=c++17 -o kmer_count kmer_count.cpp -lz -lpthread
 */

#include "kmer_common.h"
#include <getopt.h>
#include <chrono>
#include <functional>

static constexpr int PREFIX_BASES = 8;
static constexpr int PREFIX_BITS  = 2 * PREFIX_BASES;
static constexpr int NUM_BUCKETS  = 1 << PREFIX_BITS; // 65536

// --- Memory reporting (Linux) ---

static void report_memory() {
#ifdef __linux__
    FILE *f = fopen("/proc/self/status", "r");
    if (!f) return;
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (strncmp(line, "VmPeak:", 7) == 0 ||
            strncmp(line, "VmRSS:",  6) == 0)
            fprintf(stderr, "  %s", line);
    }
    fclose(f);
#endif
}

// --- MemBucketCounter: in-memory flat hash maps (default mode) ---

template<typename SuffixT>
struct MemBucketCounter {
    int k, suffix_bits;
    uint64_t suffix_mask;
    FlatHashMap<SuffixT> buckets[NUM_BUCKETS];

    explicit MemBucketCounter(int k_) : k(k_) {
        suffix_bits = 2 * k - PREFIX_BITS;
        suffix_mask = (1ULL << suffix_bits) - 1;
    }

    void operator()(uint64_t kmer) {
        uint32_t prefix = static_cast<uint32_t>(kmer >> suffix_bits);
        SuffixT suffix = static_cast<SuffixT>(kmer & suffix_mask);
        buckets[prefix].insert_or_increment(suffix);
    }

    uint64_t total_unique() const {
        uint64_t n = 0;
        for (int b = 0; b < NUM_BUCKETS; b++) n += buckets[b].size;
        return n;
    }

    void write_output(const char *path, uint32_t min_count, bool is_tsv) {
        FILE *fp;
        if (is_tsv) {
            fp = (path && path[0]) ? fopen(path, "w") : stdout;
        } else {
            fp = fopen(path, "wb");
        }
        if (!fp) { perror("fopen output"); exit(1); }

        if (!is_tsv) write_kcounts_header(fp, static_cast<uint32_t>(k), 0);

        uint64_t total = 0;
        std::vector<std::pair<SuffixT, uint32_t>> entries;

        for (int b = 0; b < NUM_BUCKETS; b++) {
            buckets[b].extract_sorted(min_count, entries);
            for (auto &[suffix, count] : entries) {
                uint64_t full = (static_cast<uint64_t>(b) << suffix_bits) | suffix;
                if (is_tsv)
                    fprintf(fp, "%lu\t%u\n", (unsigned long)full, count);
                else
                    write_kmer_entry(fp, full, count);
                total++;
            }
        }

        if (!is_tsv) {
            fseek(fp, 8, SEEK_SET);
            fwrite(&total, 8, 1, fp);
        }
        if (fp != stdout) fclose(fp);

        fprintf(stderr, "After min_count=%u filter: %lu k-mers\n",
                min_count, (unsigned long)total);
        if (is_tsv)
            fprintf(stderr, "Wrote %lu k-mers (TSV)\n", (unsigned long)total);
        else
            fprintf(stderr, "Wrote %lu k-mers to %s (binary)\n",
                    (unsigned long)total, path);
    }
};

// --- DiskBucketCounter: disk-backed mode (--lowmem) ---

template<typename SuffixT>
struct DiskBucketCounter {
    int k, suffix_bits;
    uint64_t suffix_mask;
    std::string tmpdir;

    static constexpr size_t BUF_ENTRIES = 1024;
    SuffixT *buffer_storage;
    size_t buf_pos[NUM_BUCKETS];

    DiskBucketCounter(int k_, const std::string &tmpdir_) : k(k_), tmpdir(tmpdir_) {
        suffix_bits = 2 * k - PREFIX_BITS;
        suffix_mask = (1ULL << suffix_bits) - 1;
        buffer_storage = static_cast<SuffixT *>(
            calloc(static_cast<size_t>(NUM_BUCKETS) * BUF_ENTRIES, sizeof(SuffixT)));
        if (!buffer_storage) {
            fprintf(stderr, "ERROR: DiskBucketCounter buffer alloc failed\n");
            exit(1);
        }
        memset(buf_pos, 0, sizeof(buf_pos));
    }

    ~DiskBucketCounter() {
        free(buffer_storage);
        for (int b = 0; b < NUM_BUCKETS; b++)
            remove(bucket_path(b).c_str());
        rmdir(tmpdir.c_str());
    }

    DiskBucketCounter(const DiskBucketCounter &) = delete;
    DiskBucketCounter &operator=(const DiskBucketCounter &) = delete;

    std::string bucket_path(int b) const {
        char buf[32];
        snprintf(buf, sizeof(buf), "/b%05d.tmp", b);
        return tmpdir + buf;
    }

    SuffixT *bucket_buf(int b) {
        return buffer_storage + static_cast<size_t>(b) * BUF_ENTRIES;
    }

    void flush_bucket(int b) {
        if (buf_pos[b] == 0) return;
        FILE *fp = fopen(bucket_path(b).c_str(), "ab");
        if (!fp) { perror("fopen bucket"); exit(1); }
        fwrite(bucket_buf(b), sizeof(SuffixT), buf_pos[b], fp);
        fclose(fp);
        buf_pos[b] = 0;
    }

    void operator()(uint64_t kmer) {
        uint32_t prefix = static_cast<uint32_t>(kmer >> suffix_bits);
        SuffixT suffix = static_cast<SuffixT>(kmer & suffix_mask);
        bucket_buf(prefix)[buf_pos[prefix]++] = suffix;
        if (buf_pos[prefix] == BUF_ENTRIES)
            flush_bucket(prefix);
    }

    void flush_all() {
        for (int b = 0; b < NUM_BUCKETS; b++)
            flush_bucket(b);
    }

    void write_output(const char *path, uint32_t min_count, bool is_tsv) {
        flush_all();

        FILE *fp;
        if (is_tsv) {
            fp = (path && path[0]) ? fopen(path, "w") : stdout;
        } else {
            fp = fopen(path, "wb");
        }
        if (!fp) { perror("fopen output"); exit(1); }

        if (!is_tsv) write_kcounts_header(fp, static_cast<uint32_t>(k), 0);

        uint64_t total = 0;
        std::vector<SuffixT> suffixes;

        for (int b = 0; b < NUM_BUCKETS; b++) {
            std::string bpath = bucket_path(b);
            FILE *bf = fopen(bpath.c_str(), "rb");
            if (!bf) continue;

            fseek(bf, 0, SEEK_END);
            long fsize = ftell(bf);
            fseek(bf, 0, SEEK_SET);
            size_t n = static_cast<size_t>(fsize) / sizeof(SuffixT);

            suffixes.resize(n);
            if (fread(suffixes.data(), sizeof(SuffixT), n, bf) != n) {
                fprintf(stderr, "ERROR: short read on bucket %d\n", b);
                fclose(bf);
                continue;
            }
            fclose(bf);

            std::sort(suffixes.begin(), suffixes.end());

            size_t i = 0;
            while (i < n) {
                SuffixT s = suffixes[i];
                uint32_t count = 1;
                while (i + count < n && suffixes[i + count] == s)
                    count++;
                if (count >= min_count) {
                    uint64_t full = (static_cast<uint64_t>(b) << suffix_bits) | s;
                    if (is_tsv)
                        fprintf(fp, "%lu\t%u\n", (unsigned long)full, count);
                    else
                        write_kmer_entry(fp, full, count);
                    total++;
                }
                i += count;
            }

            remove(bpath.c_str());

            if (b % 10000 == 0 && b > 0)
                fprintf(stderr, "\r  Processing bucket %d / %d", b, NUM_BUCKETS);
        }

        if (!is_tsv) {
            fseek(fp, 8, SEEK_SET);
            fwrite(&total, 8, 1, fp);
        }
        if (fp != stdout) fclose(fp);

        fprintf(stderr, "\rAfter min_count=%u filter: %lu k-mers\n",
                min_count, (unsigned long)total);
        if (is_tsv)
            fprintf(stderr, "Wrote %lu k-mers (TSV)\n", (unsigned long)total);
        else
            fprintf(stderr, "Wrote %lu k-mers to %s (binary)\n",
                    (unsigned long)total, path);
    }
};

// --- Temp directory creation ---

static std::string create_tmpdir(const std::string &base) {
    std::string tmpl = base + "/kmertools_XXXXXX";
    char *result = mkdtemp(&tmpl[0]);
    if (!result) {
        perror("mkdtemp");
        fprintf(stderr, "ERROR: Cannot create temp directory under %s\n", base.c_str());
        exit(1);
    }
    return std::string(result);
}

// --- Shared bucketed counting logic (avoids duplicating dispatch code) ---

template<typename SuffixT>
static void run_bucketed(const char *input_path, int k, uint32_t min_count,
                         bool is_tsv, const char *out_path,
                         bool lowmem, const std::string &tmpdir_base) {
    if (lowmem) {
        std::string base = tmpdir_base;
        if (base.empty()) {
            const char *env = getenv("TMPDIR");
            base = env ? env : "/tmp";
        }
        std::string tmpdir = create_tmpdir(base);
        fprintf(stderr, "  Temp directory: %s\n", tmpdir.c_str());

        DiskBucketCounter<SuffixT> counter(k, tmpdir);
        scan_kmers(input_path, k, std::ref(counter));
        counter.write_output(out_path, min_count, is_tsv);
    } else {
        MemBucketCounter<SuffixT> counter(k);
        scan_kmers(input_path, k, std::ref(counter));
        fprintf(stderr, "  %lu unique k-mers found\n",
                (unsigned long)counter.total_unique());
        counter.write_output(out_path, min_count, is_tsv);
    }
}

// --- Fallback counting for k <= 8 ---

static void run_fallback(const char *input_path, int k, int min_count,
                         bool is_tsv, const char *out_path) {
    auto counts = count_kmers(input_path, k);

    std::vector<std::pair<uint64_t, uint32_t>> filtered;
    filtered.reserve(counts.size());
    for (auto &[kmer, cnt] : counts) {
        if (static_cast<int>(cnt) >= min_count)
            filtered.push_back({kmer, cnt});
    }
    counts.clear();
    std::sort(filtered.begin(), filtered.end());

    fprintf(stderr, "After min_count=%d filter: %lu k-mers\n",
            min_count, (unsigned long)filtered.size());

    FILE *fp;
    if (is_tsv) {
        fp = (out_path && out_path[0]) ? fopen(out_path, "w") : stdout;
    } else {
        fp = fopen(out_path, "wb");
    }
    if (!fp) { perror("fopen"); exit(1); }

    if (!is_tsv) write_kcounts_header(fp, static_cast<uint32_t>(k), filtered.size());
    for (auto &[kmer, cnt] : filtered) {
        if (is_tsv)
            fprintf(fp, "%lu\t%u\n", (unsigned long)kmer, cnt);
        else
            write_kmer_entry(fp, kmer, cnt);
    }
    if (fp != stdout) fclose(fp);

    if (is_tsv)
        fprintf(stderr, "Wrote %lu k-mers (TSV)\n", (unsigned long)filtered.size());
    else
        fprintf(stderr, "Wrote %lu k-mers to %s (binary)\n",
                (unsigned long)filtered.size(), out_path);
}

// --- Usage ---

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] <input.fastq.gz|input.fasta>\n"
        "  -k INT      k-mer size (default: 21, max: 31)\n"
        "  -t INT      threads (default: 1, reserved)\n"
        "  -m INT      minimum count to output (default: 2)\n"
        "  -f STR      output format: 'tsv' or 'bin' (default: bin)\n"
        "  -o FILE     output file (required for bin, stdout for tsv)\n"
        "  --lowmem    disk-backed mode (~260 MB RAM, uses temp files)\n"
        "  --tmpdir D  temp directory for --lowmem (default: $TMPDIR or /tmp)\n"
        "  -h          show this help\n", prog);
}

// --- Main ---

int main(int argc, char **argv) {
    int k = 21;
    int min_count = 2;
    std::string out_format = "bin";
    std::string out_path;
    bool lowmem = false;
    std::string tmpdir_base;

    static struct option longopts[] = {
        {"lowmem", no_argument,       nullptr, 'L'},
        {"tmpdir", required_argument, nullptr, 'T'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "k:t:m:f:o:h", longopts, nullptr)) != -1) {
        switch (opt) {
            case 'k': k = atoi(optarg); break;
            case 't': /* threads — reserved */ break;
            case 'm': min_count = atoi(optarg); break;
            case 'f': out_format = optarg; break;
            case 'o': out_path = optarg; break;
            case 'L': lowmem = true; break;
            case 'T': tmpdir_base = optarg; break;
            case 'h': usage(argv[0]); return 0;
            default:  usage(argv[0]); return 1;
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "ERROR: No input file specified\n");
        usage(argv[0]);
        return 1;
    }
    if (k < 1 || k > 31) {
        fprintf(stderr, "ERROR: k must be between 1 and 31\n");
        return 1;
    }
    if (out_format == "bin" && out_path.empty()) {
        fprintf(stderr, "ERROR: -o required for binary output\n");
        return 1;
    }

    const char *input_path = argv[optind];
    bool is_tsv = (out_format != "bin");
    const char *out_cstr = is_tsv ? (out_path.empty() ? nullptr : out_path.c_str())
                                  : out_path.c_str();

    fprintf(stderr, "Counting %d-mers from %s", k, input_path);
    if (lowmem) fprintf(stderr, " [lowmem mode]");
    fprintf(stderr, "\n");

    auto t_start = std::chrono::steady_clock::now();

    if (k <= 8) {
        if (lowmem)
            fprintf(stderr, "  Note: --lowmem not needed for k=%d (k-mer space fits in memory), using hash map\n", k);
        run_fallback(input_path, k, min_count, is_tsv, out_cstr);
    } else if (k <= 24) {
        run_bucketed<uint32_t>(input_path, k, static_cast<uint32_t>(min_count),
                               is_tsv, out_cstr, lowmem, tmpdir_base);
    } else {
        run_bucketed<uint64_t>(input_path, k, static_cast<uint32_t>(min_count),
                               is_tsv, out_cstr, lowmem, tmpdir_base);
    }

    auto t_end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t_end - t_start).count();
    fprintf(stderr, "Wall time: %.3f sec\n", elapsed);
    report_memory();

    return 0;
}
