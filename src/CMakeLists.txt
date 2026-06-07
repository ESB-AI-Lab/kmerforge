cmake_minimum_required(VERSION 3.14)
project(kmertools LANGUAGES CXX)

find_package(ZLIB REQUIRED)

set(KMERTOOLS kmer_count kmer_diff kmer_cosinesim kmer_reads)

foreach(tool ${KMERTOOLS})
    add_executable(${tool} ${tool}.cpp)
    target_compile_features(${tool} PUBLIC cxx_std_17)
    target_link_libraries(${tool} PRIVATE ZLIB::ZLIB)
endforeach()

target_link_libraries(kmer_cosinesim PRIVATE m)
