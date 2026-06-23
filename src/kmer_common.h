/*
 * kmer_common.h — Shared k-mer utilities for the kmertools suite.
 *
 * 2-bit encoding, canonical k-mer hashing, gzip-aware file reading,
 * and the .kcounts binary format.
 *
 * MIT License — Edwin Solares / ESB AI Lab
 */

#pragma once

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <unistd.h>
#include <zlib.h>

// --- 2-bit nucleotide encoding: A=0, C=1, G=2, T=3, else=4 ---

inline uint8_t ENCODE[256];

inline struct _InitEncode {
    _InitEncode() {
        memset(ENCODE, 4, 256);
        ENCODE['A'] = ENCODE['a'] = 0;
        ENCODE['C'] = ENCODE['c'] = 1;
        ENCODE['G'] = ENCODE['g'] = 2;
        ENCODE['T'] = ENCODE['t'] = 3;
    }
} _init_encode;

// --- Canonical k-mer hashing ---

inline uint64_t revcomp(uint64_t kmer, int k) {
    uint64_t rc = 0;
    for (int i = 0; i < k; i++) {
        rc = (rc << 2) | (3 - (kmer & 3));
        kmer >>= 2;
    }
    return rc;
}

inline uint64_t canonical(uint64_t fwd, int k) {
    uint64_t rc = revcomp(fwd, k);
    return fwd < rc ? fwd : rc;
}

// --- Gzip-aware line reader ---

struct GzReader {
    gzFile fp;
    char buf[1 << 16];
    std::string line;

    GzReader(const char *path) {
        fp = gzopen(path, "r");
        if (!fp) {
            fprintf(stderr, "ERROR: Cannot open %s\n", path);
            exit(1);
        }
        gzbuffer(fp, 1 << 18);
    }
    ~GzReader() { if (fp) gzclose(fp); }

    GzReader(const GzReader&) = delete;
    GzReader& operator=(const GzReader&) = delete;

    bool next_line() {
        line.clear();
        while (true) {
            if (gzgets(fp, buf, sizeof(buf)) == nullptr)
                return !line.empty();
            line.append(buf);
            if (!line.empty() && line.back() == '\n') {
                line.pop_back();
                if (!line.empty() && line.back() == '\r')
                    line.pop_back();
                return true;
            }
        }
    }
};

// --- .kcounts binary format ---

constexpr uint32_t KCOUNTS_MAGIC = 0x4B4D4552; // "KMER"

struct KcountsHeader {
    uint32_t magic;
    uint32_t k;
    uint64_t n;
};

struct KmerEntry {
    uint64_t kmer;
    uint32_t count;
};

inline bool read_kcounts_header(FILE *fp, KcountsHeader &hdr) {
    if (fread(&hdr.magic, 4, 1, fp) != 1) return false;
    if (fread(&hdr.k, 4, 1, fp) != 1) return false;
    if (fread(&hdr.n, 8, 1, fp) != 1) return false;
    return hdr.magic == KCOUNTS_MAGIC;
}

inline void write_kcounts_header(FILE *fp, uint32_t k, uint64_t n) {
    uint32_t magic = KCOUNTS_MAGIC;
    fwrite(&magic, 4, 1, fp);
    fwrite(&k, 4, 1, fp);
    fwrite(&n, 8, 1, fp);
}

inline bool read_kmer_entry(FILE *fp, KmerEntry &e) {
    if (fread(&e.kmer, 8, 1, fp) != 1) return false;
    if (fread(&e.count, 4, 1, fp) != 1) return false;
    return true;
}

inline void write_kmer_entry(FILE *fp, uint64_t kmer, uint32_t count) {
    fwrite(&kmer, 8, 1, fp);
    fwrite(&count, 4, 1, fp);
}

// --- Format detection ---

enum FileFormat { FMT_FASTA, FMT_FASTQ };

inline FileFormat detect_format(const char *path) {
    gzFile fp = gzopen(path, "r");
    if (!fp) return FMT_FASTQ;
    int c = gzgetc(fp);
    gzclose(fp);
    return (c == '>') ? FMT_FASTA : FMT_FASTQ;
}

// --- K-mer scanning (template-based, callback receives canonical k-mers) ---

template<typename EmitFn>
void scan_kmers_fastq(const char *path, int k, EmitFn &&emit) {
    GzReader reader(path);
    uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t n_reads = 0;

    while (reader.next_line()) {
        if (reader.line.empty() || reader.line[0] != '@') continue;
        if (!reader.next_line()) break;
        const std::string &seq = reader.line;

        uint64_t fwd = 0;
        int valid = 0;
        for (size_t i = 0; i < seq.size(); i++) {
            uint8_t code = ENCODE[static_cast<uint8_t>(seq[i])];
            if (code > 3) { valid = 0; fwd = 0; continue; }
            fwd = ((fwd << 2) | code) & mask;
            valid++;
            if (valid >= k)
                emit(canonical(fwd, k));
        }
        n_reads++;

        reader.next_line(); // +
        reader.next_line(); // quality

        if (n_reads % 1000000 == 0)
            fprintf(stderr, "\r  %lu M reads processed", n_reads / 1000000);
    }
    fprintf(stderr, "\r  %lu reads processed\n", n_reads);
}

// FASTA: fwd/valid must persist across continuation lines within a contig
template<typename EmitFn>
void scan_kmers_fasta(const char *path, int k, EmitFn &&emit) {
    GzReader reader(path);
    uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t fwd = 0;
    int valid = 0;
    uint64_t n_seqs = 0;

    while (reader.next_line()) {
        if (reader.line[0] == '>') {
            valid = 0; fwd = 0;
            n_seqs++;
            if (n_seqs % 1000 == 0)
                fprintf(stderr, "\r  %lu contigs processed", n_seqs);
            continue;
        }
        const std::string &seq = reader.line;
        for (size_t i = 0; i < seq.size(); i++) {
            uint8_t code = ENCODE[static_cast<uint8_t>(seq[i])];
            if (code > 3) { valid = 0; fwd = 0; continue; }
            fwd = ((fwd << 2) | code) & mask;
            valid++;
            if (valid >= k)
                emit(canonical(fwd, k));
        }
    }
    fprintf(stderr, "\r  %lu contigs processed\n", n_seqs);
}

template<typename EmitFn>
void scan_kmers(const char *path, int k, EmitFn &&emit) {
    if (detect_format(path) == FMT_FASTA)
        scan_kmers_fasta(path, k, std::forward<EmitFn>(emit));
    else
        scan_kmers_fastq(path, k, std::forward<EmitFn>(emit));
}

// Backward-compatible wrapper using unordered_map
inline std::unordered_map<uint64_t, uint32_t> count_kmers(const char *path, int k) {
    std::unordered_map<uint64_t, uint32_t> counts;
    counts.reserve(1 << 24);
    scan_kmers(path, k, [&](uint64_t kmer) { counts[kmer]++; });
    return counts;
}

// --- FlatHashMap: open-addressing linear-probing hash map (SoA layout) ---

template<typename KeyT>
struct FlatHashMap {
    KeyT     *keys;
    uint32_t *vals;
    uint8_t  *state;
    size_t    capacity;
    size_t    size;

    FlatHashMap() : keys(nullptr), vals(nullptr), state(nullptr),
                    capacity(0), size(0) {}

    ~FlatHashMap() { release(); }

    void release() {
        free(keys); keys = nullptr;
        free(vals); vals = nullptr;
        free(state); state = nullptr;
        capacity = 0; size = 0;
    }

    FlatHashMap(const FlatHashMap &) = delete;
    FlatHashMap &operator=(const FlatHashMap &) = delete;

    FlatHashMap(FlatHashMap &&o) noexcept
        : keys(o.keys), vals(o.vals), state(o.state),
          capacity(o.capacity), size(o.size) {
        o.keys = nullptr; o.vals = nullptr; o.state = nullptr;
        o.capacity = 0; o.size = 0;
    }

    void init(size_t cap) {
        capacity = cap;
        size = 0;
        keys  = static_cast<KeyT *>(calloc(cap, sizeof(KeyT)));
        vals  = static_cast<uint32_t *>(calloc(cap, sizeof(uint32_t)));
        state = static_cast<uint8_t *>(calloc(cap, sizeof(uint8_t)));
        if (!keys || !vals || !state) {
            fprintf(stderr, "ERROR: FlatHashMap alloc failed (cap=%zu)\n", cap);
            exit(1);
        }
    }

    void insert_or_increment(KeyT key) {
        if (capacity == 0) init(1024);
        if (size * 10 >= capacity * 7) grow();
        size_t h = hash_key(key) & (capacity - 1);
        while (true) {
            if (state[h] == 0) {
                keys[h] = key; vals[h] = 1; state[h] = 1;
                size++;
                return;
            }
            if (keys[h] == key) { vals[h]++; return; }
            h = (h + 1) & (capacity - 1);
        }
    }

    bool empty() const { return size == 0; }

    uint32_t lookup(KeyT key) const {
        if (capacity == 0) return 0;
        size_t h = hash_key(key) & (capacity - 1);
        while (true) {
            if (state[h] == 0) return 0;
            if (keys[h] == key) return vals[h];
            h = (h + 1) & (capacity - 1);
        }
    }

    void extract_sorted(uint32_t min_count,
                        std::vector<std::pair<KeyT, uint32_t>> &out) const {
        out.clear();
        out.reserve(size);
        for (size_t i = 0; i < capacity; i++) {
            if (state[i] && vals[i] >= min_count)
                out.push_back({keys[i], vals[i]});
        }
        std::sort(out.begin(), out.end());
    }

private:
    void grow() {
        size_t old_cap = capacity;
        KeyT *old_keys = keys;
        uint32_t *old_vals = vals;
        uint8_t *old_state = state;

        capacity *= 2;
        keys  = static_cast<KeyT *>(calloc(capacity, sizeof(KeyT)));
        vals  = static_cast<uint32_t *>(calloc(capacity, sizeof(uint32_t)));
        state = static_cast<uint8_t *>(calloc(capacity, sizeof(uint8_t)));
        if (!keys || !vals || !state) {
            fprintf(stderr, "ERROR: FlatHashMap grow failed (cap=%zu)\n", capacity);
            exit(1);
        }
        size = 0;
        for (size_t i = 0; i < old_cap; i++) {
            if (old_state[i]) {
                size_t h = hash_key(old_keys[i]) & (capacity - 1);
                while (state[h]) h = (h + 1) & (capacity - 1);
                keys[h] = old_keys[i]; vals[h] = old_vals[i]; state[h] = 1;
                size++;
            }
        }
        free(old_keys); free(old_vals); free(old_state);
    }

    static size_t hash_key(uint32_t k) {
        k ^= k >> 16; k *= 0x45d9f3b; k ^= k >> 16;
        k *= 0x45d9f3b; k ^= k >> 16;
        return k;
    }
    static size_t hash_key(uint64_t k) {
        k ^= k >> 33; k *= 0xff51afd7ed558ccdULL; k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53ULL; k ^= k >> 33;
        return static_cast<size_t>(k);
    }
};

// --- Compact quotient hash set ---
//
// Memory-efficient static set for uint64_t keys. Uses the Jellyfish-style
// technique: table index encodes L bits of the key, only the remaining
// (64-L) bits + 1 occupied bit are stored per slot. ~5-6 bytes/entry at
// 80% load vs ~60 bytes/entry for std::unordered_set.

class CompactHashSet {
    uint32_t L_;
    uint32_t R_;
    uint32_t slot_bytes_;
    uint64_t mask_;
    size_t n_elements_;
    std::vector<uint8_t> storage_;

    static uint64_t mix(uint64_t h) {
        h ^= h >> 30;
        h *= 0xbf58476d1ce4e5b9ULL;
        h ^= h >> 27;
        h *= 0x94d049bb133111ebULL;
        h ^= h >> 31;
        return h;
    }

    static uint64_t next_pow2(uint64_t v) {
        v--;
        v |= v >> 1;  v |= v >> 2;  v |= v >> 4;
        v |= v >> 8;  v |= v >> 16; v |= v >> 32;
        return v + 1;
    }

    uint64_t read_slot(size_t idx) const {
        uint64_t val = 0;
        memcpy(&val, &storage_[idx * slot_bytes_], slot_bytes_);
        return val;
    }

    void write_slot(size_t idx, uint64_t val) {
        memcpy(&storage_[idx * slot_bytes_], &val, slot_bytes_);
    }

public:
    CompactHashSet() : L_(0), R_(64), slot_bytes_(0), mask_(0), n_elements_(0) {}

    explicit CompactHashSet(const std::vector<uint64_t> &keys, double load_factor = 0.8) {
        n_elements_ = keys.size();
        if (n_elements_ == 0) {
            L_ = 0; R_ = 64; slot_bytes_ = 0; mask_ = 0;
            return;
        }

        uint64_t table_size = next_pow2((uint64_t)(n_elements_ / load_factor) + 1);
        if (table_size < 16) table_size = 16;

        L_ = 0;
        for (uint64_t v = table_size; v > 1; v >>= 1) L_++;
        R_ = 64 - L_;
        slot_bytes_ = (R_ + 1 + 7) / 8;
        mask_ = table_size - 1;

        storage_.assign((size_t)table_size * slot_bytes_, 0);

        for (uint64_t key : keys) {
            uint64_t h = mix(key);
            uint64_t idx = h & mask_;
            uint64_t rem = h >> L_;
            uint64_t packed = (rem << 1) | 1;

            while (true) {
                uint64_t val = read_slot(idx);
                if (!(val & 1)) {
                    write_slot(idx, packed);
                    break;
                }
                if (val == packed) break;
                idx = (idx + 1) & mask_;
            }
        }
    }

    bool contains(uint64_t key) const {
        if (storage_.empty()) return false;
        uint64_t h = mix(key);
        uint64_t idx = h & mask_;
        uint64_t packed = ((h >> L_) << 1) | 1;

        while (true) {
            uint64_t val = read_slot(idx);
            if (!(val & 1)) return false;
            if (val == packed) return true;
            idx = (idx + 1) & mask_;
        }
    }

    size_t size() const { return n_elements_; }
    size_t capacity() const { return mask_ ? mask_ + 1 : 0; }
    size_t memory_bytes() const { return storage_.size(); }

    CompactHashSet(const CompactHashSet&) = delete;
    CompactHashSet& operator=(const CompactHashSet&) = delete;
    CompactHashSet(CompactHashSet&&) = default;
    CompactHashSet& operator=(CompactHashSet&&) = default;
};

// --- Bucketed k-mer counting ---
//
// Efficient k-mer counting using prefix bucketing. Two modes:
//   MemBucketCounter  — in-memory FlatHashMap per bucket (~9-13 bytes/entry)
//   DiskBucketCounter — disk-backed temp files (~260 MB RAM total)

inline constexpr int KMER_PREFIX_BASES = 8;
inline constexpr int KMER_PREFIX_BITS  = 2 * KMER_PREFIX_BASES;
inline constexpr int KMER_NUM_BUCKETS  = 1 << KMER_PREFIX_BITS; // 65536

template<typename SuffixT>
struct MemBucketCounter {
    int k, suffix_bits;
    uint64_t suffix_mask;
    FlatHashMap<SuffixT> buckets[KMER_NUM_BUCKETS];

    explicit MemBucketCounter(int k_) : k(k_) {
        suffix_bits = 2 * k - KMER_PREFIX_BITS;
        suffix_mask = (1ULL << suffix_bits) - 1;
    }

    void operator()(uint64_t kmer) {
        uint32_t prefix = static_cast<uint32_t>(kmer >> suffix_bits);
        SuffixT suffix = static_cast<SuffixT>(kmer & suffix_mask);
        buckets[prefix].insert_or_increment(suffix);
    }

    uint64_t total_unique() const {
        uint64_t n = 0;
        for (int b = 0; b < KMER_NUM_BUCKETS; b++) n += buckets[b].size;
        return n;
    }

    void extract_kmers(uint32_t min_count, uint32_t max_count,
                       std::vector<uint64_t> &out) {
        out.clear();
        std::vector<std::pair<SuffixT, uint32_t>> entries;
        for (int b = 0; b < KMER_NUM_BUCKETS; b++) {
            buckets[b].extract_sorted(min_count, entries);
            buckets[b].release();
            for (auto &[suffix, count] : entries) {
                if (max_count > 0 && count > max_count) continue;
                out.push_back((static_cast<uint64_t>(b) << suffix_bits) | suffix);
            }
        }
    }

    void partition_against(uint32_t min_count, uint32_t max_count,
                           const CompactHashSet &ref_set,
                           std::vector<uint64_t> &out_match,
                           uint64_t &out_nonmatch_count) {
        out_match.clear();
        out_nonmatch_count = 0;
        std::vector<std::pair<SuffixT, uint32_t>> entries;
        for (int b = 0; b < KMER_NUM_BUCKETS; b++) {
            buckets[b].extract_sorted(min_count, entries);
            buckets[b].release();
            for (auto &[suffix, count] : entries) {
                if (max_count > 0 && count > max_count) continue;
                uint64_t kmer = (static_cast<uint64_t>(b) << suffix_bits) | suffix;
                if (ref_set.contains(kmer))
                    out_match.push_back(kmer);
                else
                    out_nonmatch_count++;
            }
            if (b % 10000 == 0 && b > 0)
                fprintf(stderr, "\r       Partitioning: bucket %d / %d", b, KMER_NUM_BUCKETS);
        }
        if (KMER_NUM_BUCKETS > 10000)
            fprintf(stderr, "\r       Partitioning: done                          \n");
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

        for (int b = 0; b < KMER_NUM_BUCKETS; b++) {
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

inline std::string create_kmer_tmpdir(const std::string &base) {
    std::string tmpl = base + "/kmertools_XXXXXX";
    char *result = mkdtemp(&tmpl[0]);
    if (!result) {
        perror("mkdtemp");
        fprintf(stderr, "ERROR: Cannot create temp directory under %s\n", base.c_str());
        exit(1);
    }
    return std::string(result);
}

template<typename SuffixT>
struct DiskBucketCounter {
    int k, suffix_bits;
    uint64_t suffix_mask;
    std::string tmpdir;

    static constexpr size_t BUF_ENTRIES = 1024;
    SuffixT *buffer_storage;
    size_t buf_pos[KMER_NUM_BUCKETS];

    DiskBucketCounter(int k_, const std::string &tmpdir_) : k(k_), tmpdir(tmpdir_) {
        suffix_bits = 2 * k - KMER_PREFIX_BITS;
        suffix_mask = (1ULL << suffix_bits) - 1;
        buffer_storage = static_cast<SuffixT *>(
            calloc(static_cast<size_t>(KMER_NUM_BUCKETS) * BUF_ENTRIES, sizeof(SuffixT)));
        if (!buffer_storage) {
            fprintf(stderr, "ERROR: DiskBucketCounter buffer alloc failed\n");
            exit(1);
        }
        memset(buf_pos, 0, sizeof(buf_pos));
    }

    ~DiskBucketCounter() {
        free(buffer_storage);
        for (int b = 0; b < KMER_NUM_BUCKETS; b++)
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
        for (int b = 0; b < KMER_NUM_BUCKETS; b++)
            flush_bucket(b);
    }

    void extract_kmers(uint32_t min_count, uint32_t max_count,
                       std::vector<uint64_t> &out) {
        flush_all();
        out.clear();
        std::vector<SuffixT> suffixes;
        for (int b = 0; b < KMER_NUM_BUCKETS; b++) {
            std::string bpath = bucket_path(b);
            FILE *bf = fopen(bpath.c_str(), "rb");
            if (!bf) continue;
            fseek(bf, 0, SEEK_END);
            long fsize = ftell(bf);
            fseek(bf, 0, SEEK_SET);
            size_t n = static_cast<size_t>(fsize) / sizeof(SuffixT);
            suffixes.resize(n);
            if (fread(suffixes.data(), sizeof(SuffixT), n, bf) != n) {
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
                if (count >= min_count && (max_count == 0 || count <= max_count))
                    out.push_back((static_cast<uint64_t>(b) << suffix_bits) | s);
                i += count;
            }
            remove(bpath.c_str());
            if (b % 10000 == 0 && b > 0)
                fprintf(stderr, "\r  Processing bucket %d / %d", b, KMER_NUM_BUCKETS);
        }
        fprintf(stderr, "\r  Extracted %lu k-mers                    \n",
                (unsigned long)out.size());
    }

    void partition_against(uint32_t min_count, uint32_t max_count,
                           const CompactHashSet &ref_set,
                           std::vector<uint64_t> &out_match,
                           uint64_t &out_nonmatch_count) {
        flush_all();
        out_match.clear();
        out_nonmatch_count = 0;
        std::vector<SuffixT> suffixes;
        for (int b = 0; b < KMER_NUM_BUCKETS; b++) {
            std::string bpath = bucket_path(b);
            FILE *bf = fopen(bpath.c_str(), "rb");
            if (!bf) continue;
            fseek(bf, 0, SEEK_END);
            long fsize = ftell(bf);
            fseek(bf, 0, SEEK_SET);
            size_t n = static_cast<size_t>(fsize) / sizeof(SuffixT);
            suffixes.resize(n);
            if (fread(suffixes.data(), sizeof(SuffixT), n, bf) != n) {
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
                if (count >= min_count && (max_count == 0 || count <= max_count)) {
                    uint64_t kmer = (static_cast<uint64_t>(b) << suffix_bits) | s;
                    if (ref_set.contains(kmer))
                        out_match.push_back(kmer);
                    else
                        out_nonmatch_count++;
                }
                i += count;
            }
            remove(bpath.c_str());
            if (b % 10000 == 0 && b > 0)
                fprintf(stderr, "\r       Partitioning: bucket %d / %d", b, KMER_NUM_BUCKETS);
        }
        if (KMER_NUM_BUCKETS > 10000)
            fprintf(stderr, "\r       Partitioning: done                          \n");
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

        for (int b = 0; b < KMER_NUM_BUCKETS; b++) {
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
                fprintf(stderr, "\r  Processing bucket %d / %d", b, KMER_NUM_BUCKETS);
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

// --- High-level k-mer counting helpers ---

namespace detail {
template<typename SuffixT>
inline void count_kmers_bucketed(const char *path, int k,
                                  uint32_t min_count, uint32_t max_count,
                                  std::vector<uint64_t> &out,
                                  bool lowmem,
                                  const std::string &tmpdir_base) {
    if (lowmem) {
        std::string base = tmpdir_base;
        if (base.empty()) {
            const char *env = getenv("TMPDIR");
            base = env ? env : "/tmp";
        }
        std::string tmpdir = create_kmer_tmpdir(base);
        fprintf(stderr, "  Temp directory: %s\n", tmpdir.c_str());
        DiskBucketCounter<SuffixT> counter(k, tmpdir);
        scan_kmers(path, k, std::ref(counter));
        counter.extract_kmers(min_count, max_count, out);
    } else {
        MemBucketCounter<SuffixT> counter(k);
        scan_kmers(path, k, std::ref(counter));
        fprintf(stderr, "  %lu unique k-mers found\n",
                (unsigned long)counter.total_unique());
        counter.extract_kmers(min_count, max_count, out);
        fprintf(stderr, "  %lu k-mers after filtering\n", (unsigned long)out.size());
    }
}
template<typename SuffixT>
inline void count_and_partition(const char *path, int k,
                                uint32_t min_count, uint32_t max_count,
                                const CompactHashSet &ref_set,
                                std::vector<uint64_t> &out_match,
                                uint64_t &out_nonmatch_count,
                                bool lowmem,
                                const std::string &tmpdir_base) {
    if (lowmem) {
        std::string base = tmpdir_base;
        if (base.empty()) {
            const char *env = getenv("TMPDIR");
            base = env ? env : "/tmp";
        }
        std::string tmpdir = create_kmer_tmpdir(base);
        fprintf(stderr, "  Temp directory: %s\n", tmpdir.c_str());
        DiskBucketCounter<SuffixT> counter(k, tmpdir);
        scan_kmers(path, k, std::ref(counter));
        counter.partition_against(min_count, max_count, ref_set,
                                  out_match, out_nonmatch_count);
    } else {
        MemBucketCounter<SuffixT> counter(k);
        scan_kmers(path, k, std::ref(counter));
        fprintf(stderr, "  %lu unique k-mers found\n",
                (unsigned long)counter.total_unique());
        counter.partition_against(min_count, max_count, ref_set,
                                  out_match, out_nonmatch_count);
    }
}
} // namespace detail

inline void count_kmers_filtered(const char *path, int k,
                                  uint32_t min_count, uint32_t max_count,
                                  std::vector<uint64_t> &out,
                                  bool lowmem = false,
                                  const std::string &tmpdir_base = "") {
    if (k <= 8) {
        auto counts = count_kmers(path, k);
        uint64_t total = counts.size();
        out.clear();
        out.reserve(total);
        for (auto &[kmer, cnt] : counts) {
            if (cnt >= min_count && (max_count == 0 || cnt <= max_count))
                out.push_back(kmer);
        }
        fprintf(stderr, "  %lu unique k-mers, %lu after filtering\n",
                (unsigned long)total, (unsigned long)out.size());
    } else if (k <= 24) {
        detail::count_kmers_bucketed<uint32_t>(path, k, min_count, max_count,
                                                out, lowmem, tmpdir_base);
    } else {
        detail::count_kmers_bucketed<uint64_t>(path, k, min_count, max_count,
                                                out, lowmem, tmpdir_base);
    }
}

inline void count_and_partition(const char *path, int k,
                                 uint32_t min_count, uint32_t max_count,
                                 const CompactHashSet &ref_set,
                                 std::vector<uint64_t> &out_match,
                                 uint64_t &out_nonmatch_count,
                                 bool lowmem = false,
                                 const std::string &tmpdir_base = "") {
    if (k <= 8) {
        auto counts = count_kmers(path, k);
        out_match.clear();
        out_nonmatch_count = 0;
        for (auto &[kmer, cnt] : counts) {
            if (cnt >= min_count && (max_count == 0 || cnt <= max_count)) {
                if (ref_set.contains(kmer))
                    out_match.push_back(kmer);
                else
                    out_nonmatch_count++;
            }
        }
    } else if (k <= 24) {
        detail::count_and_partition<uint32_t>(path, k, min_count, max_count,
                                              ref_set, out_match, out_nonmatch_count,
                                              lowmem, tmpdir_base);
    } else {
        detail::count_and_partition<uint64_t>(path, k, min_count, max_count,
                                              ref_set, out_match, out_nonmatch_count,
                                              lowmem, tmpdir_base);
    }
}

inline void load_kcounts_kmers(const char *path, int expected_k,
                                uint32_t min_count, uint32_t max_count,
                                std::vector<uint64_t> &out,
                                uint64_t *total_out = nullptr) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open .kcounts file: %s\n", path);
        exit(1);
    }
    KcountsHeader hdr;
    if (!read_kcounts_header(fp, hdr)) {
        fprintf(stderr, "ERROR: Invalid .kcounts file: %s\n", path);
        fclose(fp);
        exit(1);
    }
    if ((int)hdr.k != expected_k) {
        fprintf(stderr, "ERROR: .kcounts k=%u but expected k=%d in file: %s\n",
                hdr.k, expected_k, path);
        fclose(fp);
        exit(1);
    }
    out.clear();
    if (hdr.n > 0) out.reserve(hdr.n);
    uint64_t total = 0;
    KmerEntry e;
    while (read_kmer_entry(fp, e)) {
        total++;
        if (e.count >= min_count && (max_count == 0 || e.count <= max_count))
            out.push_back(e.kmer);
    }
    fclose(fp);
    if (total_out) *total_out = total;
    fprintf(stderr, "  Loaded %lu k-mers, %lu after filtering (k=%d)\n",
            (unsigned long)total, (unsigned long)out.size(), expected_k);
}

// --- K-mer sequence decoding ---

inline std::string decode_kmer(uint64_t val, int k) {
    static constexpr char bases[] = "ACGT";
    std::string s(k, 'N');
    for (int i = k - 1; i >= 0; i--) {
        s[i] = bases[val & 3];
        val >>= 2;
    }
    return s;
}
