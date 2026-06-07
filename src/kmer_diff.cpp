/*
 * kmer_diff.cpp — Extract differential k-mers between two sorted .kcounts files.
 *
 * Streaming merge-join: O(n+m) time, O(1) memory.
 * Outputs k-mers unique to sample A, unique to sample B, or shared with
 * significantly different abundance.
 *
 * MIT License — Edwin Solares / ESB AI Lab
 *
 * Usage:
 *   ./kmer_diff -a motherhass.kcounts -b mendez.kcounts -o mendez_unique.kcounts
 *
 * Modes:
 *   --only-a    : k-mers in A but not B (lost in B)
 *   --only-b    : k-mers in B but not A (gained in B)
 *   --diff-ab   : k-mers in both but count ratio > threshold
 *   --all       : all three categories (default)
 *
 * Output: binary .kcounts file (hash+count pairs), loadable by kmer_reads
 *
 * Compile:
 *   g++ -O3 -std=c++17 -o kmer_diff kmer_diff.cpp -lz
 */

#include "kmer_common.h"
#include <cmath>
#include <getopt.h>

int main(int argc, char **argv) {
    const char *path_a = nullptr;
    const char *path_b = nullptr;
    const char *out_path = nullptr;
    const char *out_tsv = nullptr;
    double ratio_thresh = 3.0;
    int mode = 0; // 0=all, 1=only-a, 2=only-b, 3=diff-ab
    bool show_seq = false;

    static struct option longopts[] = {
        {"only-a",  no_argument,       nullptr, 'A'},
        {"only-b",  no_argument,       nullptr, 'B'},
        {"diff-ab", no_argument,       nullptr, 'D'},
        {"all",     no_argument,       nullptr, 'L'},
        {"seq",     no_argument,       nullptr, 'S'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "a:b:o:t:r:h", longopts, nullptr)) != -1) {
        switch (opt) {
            case 'a': path_a = optarg; break;
            case 'b': path_b = optarg; break;
            case 'o': out_path = optarg; break;
            case 't': out_tsv = optarg; break;
            case 'r': ratio_thresh = atof(optarg); break;
            case 'A': mode = 1; break;
            case 'B': mode = 2; break;
            case 'D': mode = 3; break;
            case 'L': mode = 0; break;
            case 'S': show_seq = true; break;
            case 'h':
            default:
                fprintf(stderr,
                    "Usage: %s -a ref.kcounts -b query.kcounts -o diff.kcounts [-t diff.tsv]\n"
                    "  -a FILE       Reference sample (e.g. Mother Hass)\n"
                    "  -b FILE       Query sample (e.g. Carmen)\n"
                    "  -o FILE       Binary output (.kcounts format: hash+count pairs)\n"
                    "  -t FILE       TSV output with counts and categories\n"
                    "  -r FLOAT      Count ratio threshold for 'diff' category (default: 3.0)\n"
                    "  --only-a      K-mers in A but not B (lost)\n"
                    "  --only-b      K-mers in B but not A (gained)\n"
                    "  --diff-ab     K-mers in both but ratio > threshold\n"
                    "  --all         All categories (default)\n"
                    "  --seq         Include nucleotide sequence in TSV\n", argv[0]);
                return opt == 'h' ? 0 : 1;
        }
    }

    if (!path_a || !path_b) {
        fprintf(stderr, "ERROR: -a and -b required\n");
        return 1;
    }

    FILE *fa = fopen(path_a, "rb");
    FILE *fb = fopen(path_b, "rb");
    if (!fa) { fprintf(stderr, "Cannot open %s\n", path_a); return 1; }
    if (!fb) { fprintf(stderr, "Cannot open %s\n", path_b); return 1; }

    KcountsHeader ha, hb;
    if (!read_kcounts_header(fa, ha) || !read_kcounts_header(fb, hb)) {
        fprintf(stderr, "Invalid kcounts file\n");
        return 1;
    }
    if (ha.k != hb.k) {
        fprintf(stderr, "ERROR: k-mer sizes don't match (%u vs %u)\n", ha.k, hb.k);
        return 1;
    }
    int k = static_cast<int>(ha.k);

    FILE *fout_bin = nullptr;
    FILE *fout_tsv = nullptr;

    if (out_path) {
        fout_bin = fopen(out_path, "wb");
        if (!fout_bin) { perror("fopen output"); return 1; }
        write_kcounts_header(fout_bin, ha.k, 0); // placeholder count
    }

    if (out_tsv) {
        fout_tsv = fopen(out_tsv, "w");
        if (!fout_tsv) { perror("fopen tsv"); return 1; }
        if (show_seq)
            fprintf(fout_tsv, "kmer_hash\tsequence\tcount_a\tcount_b\tcategory\n");
        else
            fprintf(fout_tsv, "kmer_hash\tcount_a\tcount_b\tcategory\n");
    }

    KmerEntry ea, eb;
    bool has_a = read_kmer_entry(fa, ea);
    bool has_b = read_kmer_entry(fb, eb);

    uint64_t only_a_count = 0;
    uint64_t only_b_count = 0;
    uint64_t diff_count = 0;
    uint64_t shared_same = 0;
    uint64_t total_written = 0;

    auto emit = [&](uint64_t kmer, uint32_t cnt) {
        if (fout_bin) { write_kmer_entry(fout_bin, kmer, cnt); total_written++; }
    };
    auto emit_tsv = [&](uint64_t kmer, uint32_t ca, uint32_t cb, const char *cat) {
        if (!fout_tsv) return;
        if (show_seq)
            fprintf(fout_tsv, "%lu\t%s\t%u\t%u\t%s\n",
                    (unsigned long)kmer, decode_kmer(kmer, k).c_str(), ca, cb, cat);
        else
            fprintf(fout_tsv, "%lu\t%u\t%u\t%s\n", (unsigned long)kmer, ca, cb, cat);
    };

    while (has_a && has_b) {
        if (ea.kmer == eb.kmer) {
            double ratio = (ea.count > eb.count)
                ? static_cast<double>(ea.count) / eb.count
                : static_cast<double>(eb.count) / ea.count;
            if (ratio >= ratio_thresh && (mode == 0 || mode == 3)) {
                diff_count++;
                emit(ea.kmer, ea.count);
                emit_tsv(ea.kmer, ea.count, eb.count, "diff");
            } else {
                shared_same++;
            }
            has_a = read_kmer_entry(fa, ea);
            has_b = read_kmer_entry(fb, eb);
        } else if (ea.kmer < eb.kmer) {
            only_a_count++;
            if (mode == 0 || mode == 1) {
                emit(ea.kmer, ea.count);
                emit_tsv(ea.kmer, ea.count, 0, "only_a");
            }
            has_a = read_kmer_entry(fa, ea);
        } else {
            only_b_count++;
            if (mode == 0 || mode == 2) {
                emit(eb.kmer, eb.count);
                emit_tsv(eb.kmer, 0, eb.count, "only_b");
            }
            has_b = read_kmer_entry(fb, eb);
        }
    }

    while (has_a) {
        only_a_count++;
        if (mode == 0 || mode == 1) {
            emit(ea.kmer, ea.count);
            emit_tsv(ea.kmer, ea.count, 0, "only_a");
        }
        has_a = read_kmer_entry(fa, ea);
    }
    while (has_b) {
        only_b_count++;
        if (mode == 0 || mode == 2) {
            emit(eb.kmer, eb.count);
            emit_tsv(eb.kmer, 0, eb.count, "only_b");
        }
        has_b = read_kmer_entry(fb, eb);
    }

    fclose(fa);
    fclose(fb);

    if (fout_bin) {
        fseek(fout_bin, 8, SEEK_SET);
        fwrite(&total_written, 8, 1, fout_bin);
        fclose(fout_bin);
    }
    if (fout_tsv) fclose(fout_tsv);

    fprintf(stdout, "=== K-MER DIFF SUMMARY ===\n");
    fprintf(stdout, "  A: %s (%lu k-mers)\n", path_a, (unsigned long)ha.n);
    fprintf(stdout, "  B: %s (%lu k-mers)\n", path_b, (unsigned long)hb.n);
    fprintf(stdout, "  Shared (same count) : %lu\n", (unsigned long)shared_same);
    fprintf(stdout, "  Shared (diff count) : %lu (ratio >= %.1f)\n",
            (unsigned long)diff_count, ratio_thresh);
    fprintf(stdout, "  Only in A (lost)    : %lu\n", (unsigned long)only_a_count);
    fprintf(stdout, "  Only in B (gained)  : %lu\n", (unsigned long)only_b_count);
    fprintf(stdout, "  Total differential  : %lu\n",
            (unsigned long)(only_a_count + only_b_count + diff_count));
    if (out_path)
        fprintf(stdout, "  Written to          : %s (%lu entries)\n",
                out_path, (unsigned long)total_written);
    if (out_tsv)
        fprintf(stdout, "  TSV written to      : %s\n", out_tsv);

    return 0;
}
