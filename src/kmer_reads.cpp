/*
 * kmer_reads.cpp — Extract FASTQ reads containing differential k-mers.
 *
 * Loads a set of k-mer hashes (from kmer_diff output), then streams
 * through a FASTQ.gz and outputs reads containing >=1 matching k-mer.
 *
 * MIT License — Edwin Solares / ESB AI Lab
 *
 * Usage:
 *   ./kmer_reads -k 21 -d diff.kcounts -i sample.fastq.gz -o unique_reads.fastq.gz
 *   ./kmer_reads -k 21 -d diff.kcounts -i sample.fastq.gz --invert -o shared_reads.fastq.gz
 *
 * Compile:
 *   g++ -O3 -std=c++17 -o kmer_reads kmer_reads.cpp -lz
 */

#include "kmer_common.h"
#include <unordered_set>
#include <getopt.h>

static size_t load_kmer_set(const char *path, std::unordered_set<uint64_t> &kset) {
    FILE *fp = fopen(path, "rb");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); exit(1); }

    KcountsHeader hdr;
    if (!read_kcounts_header(fp, hdr)) {
        fprintf(stderr, "Invalid kcounts file: %s\n", path);
        exit(1);
    }

    kset.reserve(hdr.n);
    KmerEntry e;
    for (uint64_t i = 0; i < hdr.n; i++) {
        if (!read_kmer_entry(fp, e)) break;
        kset.insert(e.kmer);
    }
    fclose(fp);
    return static_cast<size_t>(hdr.n);
}

static bool read_has_kmer(const std::string &seq, int k,
                          const std::unordered_set<uint64_t> &kset) {
    uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t fwd = 0;
    int valid = 0;

    for (size_t i = 0; i < seq.size(); i++) {
        uint8_t code = ENCODE[static_cast<uint8_t>(seq[i])];
        if (code > 3) { valid = 0; fwd = 0; continue; }
        fwd = ((fwd << 2) | code) & mask;
        valid++;
        if (valid >= k) {
            if (kset.count(canonical(fwd, k)))
                return true;
        }
    }
    return false;
}

int main(int argc, char **argv) {
    int k = 21;
    const char *diff_path = nullptr;
    const char *in_path = nullptr;
    const char *out_path = nullptr;
    bool invert = false;

    static struct option longopts[] = {
        {"invert", no_argument, nullptr, 'v'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "k:d:i:o:vh", longopts, nullptr)) != -1) {
        switch (opt) {
            case 'k': k = atoi(optarg); break;
            case 'd': diff_path = optarg; break;
            case 'i': in_path = optarg; break;
            case 'o': out_path = optarg; break;
            case 'v': invert = true; break;
            case 'h':
            default:
                fprintf(stderr,
                    "Usage: %s -k 21 -d diff.kcounts -i reads.fastq.gz -o filtered.fastq.gz\n"
                    "  -k INT   k-mer size (must match diff file)\n"
                    "  -d FILE  differential k-mers (.kcounts from kmer_diff)\n"
                    "  -i FILE  input FASTQ (gzipped or plain)\n"
                    "  -o FILE  output FASTQ (gzipped)\n"
                    "  --invert output reads that do NOT match (shared reads)\n", argv[0]);
                return opt == 'h' ? 0 : 1;
        }
    }

    if (!diff_path || !in_path || !out_path) {
        fprintf(stderr, "ERROR: -d, -i, and -o required\n");
        return 1;
    }

    std::unordered_set<uint64_t> kset;
    fprintf(stderr, "Loading differential k-mers from %s...\n", diff_path);
    size_t nloaded = load_kmer_set(diff_path, kset);
    fprintf(stderr, "Loaded %lu differential k-mers\n", (unsigned long)nloaded);

    GzReader reader(in_path);
    gzFile out = gzopen(out_path, "w");
    if (!out) { perror("gzopen output"); return 1; }

    uint64_t total_reads = 0;
    uint64_t matched_reads = 0;
    uint64_t output_reads = 0;
    std::string header, seq, sep, qual;

    while (reader.next_line()) {
        header = reader.line;
        if (header.empty() || header[0] != '@') continue;
        if (!reader.next_line()) break;
        seq = reader.line;
        if (!reader.next_line()) break;
        sep = reader.line;
        if (!reader.next_line()) break;
        qual = reader.line;

        total_reads++;
        bool has_match = read_has_kmer(seq, k, kset);
        if (has_match) matched_reads++;

        bool should_output = invert ? !has_match : has_match;
        if (should_output) {
            gzprintf(out, "%s\n%s\n%s\n%s\n",
                     header.c_str(), seq.c_str(),
                     sep.c_str(), qual.c_str());
            output_reads++;
        }

        if (total_reads % 500000 == 0)
            fprintf(stderr, "\r  %luM reads, %lu matched, %lu output",
                    total_reads / 1000000, matched_reads, output_reads);
    }

    gzclose(out);

    fprintf(stderr, "\n");
    fprintf(stdout, "=== READ EXTRACTION SUMMARY ===\n");
    fprintf(stdout, "  Input reads   : %lu\n", (unsigned long)total_reads);
    fprintf(stdout, "  Matched reads : %lu (%.2f%%)\n",
            (unsigned long)matched_reads,
            100.0 * matched_reads / (total_reads ? total_reads : 1));
    fprintf(stdout, "  Output reads  : %lu (%s)\n",
            (unsigned long)output_reads,
            invert ? "inverted — shared reads" : "differential reads");
    fprintf(stdout, "  Output file   : %s\n", out_path);

    return 0;
}
