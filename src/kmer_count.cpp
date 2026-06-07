/*
 * kmer_count.cpp — Fast canonical k-mer counter for FASTQ/FASTA (gzipped or plain)
 *
 * MIT License — Edwin Solares / ESB AI Lab
 *
 * Usage:
 *   ./kmer_count -k 21 -t 1 -m 2 -o output.kcounts input.fastq.gz
 *
 * Options:
 *   -k INT   k-mer size (default: 21, max: 31)
 *   -t INT   threads for reading (default: 1)
 *   -m INT   minimum count to output (default: 2)
 *   -f tsv   output format: 'tsv' or 'bin' (default: bin)
 *   -o FILE  output file (default: stdout for tsv, required for bin)
 *
 * Compile:
 *   g++ -O3 -std=c++17 -o kmer_count kmer_count.cpp -lz -lpthread
 */

#include "kmer_common.h"
#include <getopt.h>

void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options] <input.fastq.gz|input.fasta>\n"
        "  -k INT   k-mer size (default: 21, max: 31)\n"
        "  -t INT   threads (default: 1)\n"
        "  -m INT   minimum count to output (default: 2)\n"
        "  -f STR   output format: 'tsv' or 'bin' (default: bin)\n"
        "  -o FILE  output file (required for bin, stdout for tsv)\n"
        "  -h       show this help\n", prog);
}

int main(int argc, char **argv) {
    int k = 21;
    int min_count = 2;
    std::string out_format = "bin";
    std::string out_path;

    int opt;
    while ((opt = getopt(argc, argv, "k:t:m:f:o:h")) != -1) {
        switch (opt) {
            case 'k': k = atoi(optarg); break;
            case 't': /* threads — reserved */ break;
            case 'm': min_count = atoi(optarg); break;
            case 'f': out_format = optarg; break;
            case 'o': out_path = optarg; break;
            case 'h': usage(argv[0]); return 0;
            default: usage(argv[0]); return 1;
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

    const char *input_path = argv[optind];
    fprintf(stderr, "Counting %d-mers from %s\n", k, input_path);

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

    if (out_format == "bin") {
        if (out_path.empty()) {
            fprintf(stderr, "ERROR: -o required for binary output\n");
            return 1;
        }
        FILE *fp = fopen(out_path.c_str(), "wb");
        if (!fp) { perror("fopen"); return 1; }
        write_kcounts_header(fp, static_cast<uint32_t>(k), filtered.size());
        for (auto &[kmer, cnt] : filtered)
            write_kmer_entry(fp, kmer, cnt);
        fclose(fp);
        fprintf(stderr, "Wrote %lu k-mers to %s (binary)\n",
                (unsigned long)filtered.size(), out_path.c_str());
    } else {
        FILE *fp = out_path.empty() ? stdout : fopen(out_path.c_str(), "w");
        if (!fp) { perror("fopen"); return 1; }
        for (auto &[kmer, cnt] : filtered)
            fprintf(fp, "%lu\t%u\n", (unsigned long)kmer, cnt);
        if (fp != stdout) fclose(fp);
        fprintf(stderr, "Wrote %lu k-mers (TSV)\n",
                (unsigned long)filtered.size());
    }

    return 0;
}
