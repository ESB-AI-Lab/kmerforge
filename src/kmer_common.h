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

// --- K-mer counting ---

inline void count_kmers_fastq(const char *path, int k,
                              std::unordered_map<uint64_t, uint32_t> &counts) {
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
                counts[canonical(fwd, k)]++;
        }
        n_reads++;

        reader.next_line(); // +
        reader.next_line(); // quality

        if (n_reads % 1000000 == 0)
            fprintf(stderr, "\r  %lu M reads processed", n_reads / 1000000);
    }
    fprintf(stderr, "\r  %lu reads processed, %lu unique k-mers\n",
            n_reads, (unsigned long)counts.size());
}

inline void count_kmers_fasta(const char *path, int k,
                              std::unordered_map<uint64_t, uint32_t> &counts) {
    GzReader reader(path);
    uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t fwd = 0;
    int valid = 0;
    uint64_t n_seqs = 0;

    while (reader.next_line()) {
        if (reader.line[0] == '>') {
            valid = 0;
            fwd = 0;
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
                counts[canonical(fwd, k)]++;
        }
    }
    fprintf(stderr, "\r  %lu contigs processed, %lu unique k-mers\n",
            n_seqs, (unsigned long)counts.size());
}

inline std::unordered_map<uint64_t, uint32_t> count_kmers(const char *path, int k) {
    std::unordered_map<uint64_t, uint32_t> counts;
    counts.reserve(1 << 24);
    if (detect_format(path) == FMT_FASTA)
        count_kmers_fasta(path, k, counts);
    else
        count_kmers_fastq(path, k, counts);
    return counts;
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
