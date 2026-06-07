CXX      ?= g++
CXXFLAGS ?= -O3 -std=c++17
LDLIBS    = -lz

SRCDIR   = src
TOOLS    = kmer_count kmer_diff kmer_cosinesim kmer_reads

all: $(TOOLS)

kmer_count: $(SRCDIR)/kmer_count.cpp $(SRCDIR)/kmer_common.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

kmer_diff: $(SRCDIR)/kmer_diff.cpp $(SRCDIR)/kmer_common.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

kmer_cosinesim: $(SRCDIR)/kmer_cosinesim.cpp $(SRCDIR)/kmer_common.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS) -lm

kmer_reads: $(SRCDIR)/kmer_reads.cpp $(SRCDIR)/kmer_common.h
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(TOOLS)

.PHONY: all clean
