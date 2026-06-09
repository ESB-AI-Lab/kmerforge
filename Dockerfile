FROM debian:bookworm-slim AS build

RUN apt-get update && \
    apt-get install -y --no-install-recommends g++ make zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY src/ src/
COPY Makefile .

RUN make -j$(nproc)

FROM debian:bookworm-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends zlib1g && \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /build/kmer_count /build/kmer_diff \
     /build/kmer_cosinesim /build/kmer_reads /usr/local/bin/

CMD ["/bin/bash"]
