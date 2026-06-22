#include "src/kmer_common.h"
#include <cstdio>
#include <cstdlib>
#include <algorithm>

#define ASSERT(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); return 1; } } while(0)

int main() {
    const int k = 11;
    const uint32_t min_count = 1;
    const uint32_t max_count = 0;

    // Count A k-mers and build reference set
    std::vector<uint64_t> a_kmers;
    count_kmers_filtered("test/sample_a.fastq", k, min_count, max_count, a_kmers);
    ASSERT(a_kmers.size() == 12, "A should have 12 k-mers");
    CompactHashSet ref_set(a_kmers);

    // Partition B against A (in-memory)
    std::vector<uint64_t> match_mem;
    uint64_t nonmatch_mem = 0;
    count_and_partition("test/sample_b.fastq", k, min_count, max_count,
                        ref_set, match_mem, nonmatch_mem);
    ASSERT(match_mem.size() == 1, "should have 1 shared k-mer (mem)");
    ASSERT(nonmatch_mem == 19, "should have 19 non-matching k-mers (mem)");

    // Partition B against A (lowmem)
    std::vector<uint64_t> match_low;
    uint64_t nonmatch_low = 0;
    count_and_partition("test/sample_b.fastq", k, min_count, max_count,
                        ref_set, match_low, nonmatch_low, true);
    ASSERT(match_low.size() == 1, "should have 1 shared k-mer (lowmem)");
    ASSERT(nonmatch_low == 19, "should have 19 non-matching k-mers (lowmem)");

    // Verify mem and lowmem produce identical match sets
    std::sort(match_mem.begin(), match_mem.end());
    std::sort(match_low.begin(), match_low.end());
    ASSERT(match_mem == match_low, "mem and lowmem match sets should be identical");

    // Cross-check: count B independently and intersect manually
    std::vector<uint64_t> b_kmers;
    count_kmers_filtered("test/sample_b.fastq", k, min_count, max_count, b_kmers);
    std::sort(b_kmers.begin(), b_kmers.end());

    std::vector<uint64_t> expected_shared;
    std::set_intersection(a_kmers.begin(), a_kmers.end(),
                          b_kmers.begin(), b_kmers.end(),
                          std::back_inserter(expected_shared));
    ASSERT(match_mem == expected_shared, "partition matches should equal set_intersection result");
    ASSERT(nonmatch_mem == b_kmers.size() - expected_shared.size(),
           "non-match count should equal B minus shared");

    return 0;
}
