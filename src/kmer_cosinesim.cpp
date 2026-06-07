/*
 * kmer_cosine.cpp — Pairwise cosine similarity from sorted .kcounts files.
 *
 * Uses merge-join on sorted k-mer arrays — O(n+m) time, O(1) extra memory.
 * No hash tables. Reads are streamed, not loaded entirely into RAM.
 *
 * MIT License — Edwin Solares / ESB AI Lab
 *
 * Usage:
 *   ./kmer_cosine file1.kcounts file2.kcounts [file3.kcounts ...]
 *
 * Compile:
 *   g++ -O3 -std=c++17 -o kmer_cosine kmer_cosinesim.cpp -lz -lm
 */

#include "kmer_common.h"
#include <cmath>

struct PairResult {
    double cosine = 0.0;
    double jaccard = 0.0;
    double angular = 0.0;
    uint64_t shared = 0;
    uint64_t union_size = 0;
    uint64_t n_a = 0;
    uint64_t n_b = 0;
};

static PairResult compute_pair(const char *path_a, const char *path_b) {
    PairResult res{};

    FILE *fa = fopen(path_a, "rb");
    FILE *fb = fopen(path_b, "rb");
    if (!fa) { fprintf(stderr, "Cannot open %s\n", path_a); return res; }
    if (!fb) { fprintf(stderr, "Cannot open %s\n", path_b); fclose(fa); return res; }

    KcountsHeader ha, hb;
    if (!read_kcounts_header(fa, ha) || !read_kcounts_header(fb, hb)) {
        fclose(fa); fclose(fb);
        return res;
    }

    res.n_a = ha.n;
    res.n_b = hb.n;

    double dot = 0.0;
    double norm_a_sq = 0.0;
    double norm_b_sq = 0.0;
    uint64_t shared = 0;
    uint64_t only_a = 0;
    uint64_t only_b = 0;

    KmerEntry ea, eb;
    bool has_a = read_kmer_entry(fa, ea);
    bool has_b = read_kmer_entry(fb, eb);

    while (has_a && has_b) {
        if (ea.kmer == eb.kmer) {
            dot += static_cast<double>(ea.count) * static_cast<double>(eb.count);
            norm_a_sq += static_cast<double>(ea.count) * static_cast<double>(ea.count);
            norm_b_sq += static_cast<double>(eb.count) * static_cast<double>(eb.count);
            shared++;
            has_a = read_kmer_entry(fa, ea);
            has_b = read_kmer_entry(fb, eb);
        } else if (ea.kmer < eb.kmer) {
            norm_a_sq += static_cast<double>(ea.count) * static_cast<double>(ea.count);
            only_a++;
            has_a = read_kmer_entry(fa, ea);
        } else {
            norm_b_sq += static_cast<double>(eb.count) * static_cast<double>(eb.count);
            only_b++;
            has_b = read_kmer_entry(fb, eb);
        }
    }

    while (has_a) {
        norm_a_sq += static_cast<double>(ea.count) * static_cast<double>(ea.count);
        only_a++;
        has_a = read_kmer_entry(fa, ea);
    }
    while (has_b) {
        norm_b_sq += static_cast<double>(eb.count) * static_cast<double>(eb.count);
        only_b++;
        has_b = read_kmer_entry(fb, eb);
    }

    fclose(fa);
    fclose(fb);

    double norm_a = sqrt(norm_a_sq);
    double norm_b = sqrt(norm_b_sq);

    res.cosine = (norm_a > 0 && norm_b > 0) ? dot / (norm_a * norm_b) : 0.0;
    res.shared = shared;
    res.union_size = shared + only_a + only_b;
    res.jaccard = res.union_size > 0 ? static_cast<double>(shared) / static_cast<double>(res.union_size) : 0.0;
    res.angular = acos(fmin(fmax(res.cosine, -1.0), 1.0)) * 180.0 / M_PI;

    return res;
}

static std::string basename_no_ext(const char *path) {
    std::string s(path);
    if (auto slash = s.rfind('/'); slash != std::string::npos)
        s = s.substr(slash + 1);
    if (auto dot = s.rfind(".kcounts"); dot != std::string::npos)
        s = s.substr(0, dot);
    return s;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <file1.kcounts> <file2.kcounts> [...]\n", argv[0]);
        return 1;
    }

    int n = argc - 1;
    std::vector<std::string> names(n);
    std::vector<const char*> paths(n);

    for (int i = 0; i < n; i++) {
        paths[i] = argv[i + 1];
        names[i] = basename_no_ext(paths[i]);
    }

    fprintf(stdout, "=== K-MER COUNTS ===\n");
    for (int i = 0; i < n; i++) {
        FILE *f = fopen(paths[i], "rb");
        if (!f) { fprintf(stderr, "Cannot open %s\n", paths[i]); continue; }
        KcountsHeader h;
        read_kcounts_header(f, h);
        fclose(f);
        fprintf(stdout, "  %-25s  k=%u  kmers=%lu\n",
                names[i].c_str(), h.k, (unsigned long)h.n);
    }

    fprintf(stdout, "\n=== PAIRWISE SIMILARITY ===\n\n");
    fprintf(stdout, "%-20s  %-20s  %10s  %10s  %10s  %12s  %12s\n",
            "Sample_A", "Sample_B", "Cosine", "Jaccard", "Angular", "Shared", "Union");
    fprintf(stdout, "%-20s  %-20s  %10s  %10s  %10s  %12s  %12s\n",
            "--------", "--------", "------", "-------", "-------", "------", "-----");

    std::vector<std::vector<double>> matrix(n, std::vector<double>(n, 1.0));

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            fprintf(stderr, "Comparing %s vs %s...\n",
                    names[i].c_str(), names[j].c_str());

            auto r = compute_pair(paths[i], paths[j]);

            fprintf(stdout, "%-20s  %-20s  %10.6f  %10.6f  %10.4f  %12lu  %12lu\n",
                    names[i].c_str(), names[j].c_str(),
                    r.cosine, r.jaccard, r.angular,
                    (unsigned long)r.shared, (unsigned long)r.union_size);

            matrix[i][j] = r.cosine;
            matrix[j][i] = r.cosine;
        }
    }

    fprintf(stdout, "\n=== COSINE SIMILARITY MATRIX ===\n\n");

    int maxlen = 8;
    for (int i = 0; i < n; i++)
        if (static_cast<int>(names[i].size()) > maxlen)
            maxlen = static_cast<int>(names[i].size());

    fprintf(stdout, "%*s", maxlen, "");
    for (int j = 0; j < n; j++)
        fprintf(stdout, "  %*s", maxlen, names[j].c_str());
    fprintf(stdout, "\n");

    for (int i = 0; i < n; i++) {
        fprintf(stdout, "%*s", maxlen, names[i].c_str());
        for (int j = 0; j < n; j++)
            fprintf(stdout, "  %*.6f", maxlen, matrix[i][j]);
        fprintf(stdout, "\n");
    }

    fprintf(stdout, "\n=== INTERPRETATION ===\n\n");
    fprintf(stdout, "  Cosine ~1.000   : Clonal / identical\n");
    fprintf(stdout, "  Cosine 0.99+    : Very closely related (sport / bud mutation)\n");
    fprintf(stdout, "  Cosine 0.95-0.99: Close relatives (seedling, same race)\n");
    fprintf(stdout, "  Cosine < 0.95   : Different cultivars / races\n");
    fprintf(stdout, "\n  Angular distance: 0 = identical, 90 = orthogonal\n");
    fprintf(stdout, "  (Better resolution than cosine at high similarity)\n\n");

    return 0;
}
