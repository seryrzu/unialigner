project(minseq)

add_library(minseq_library STATIC minseq.cpp)

add_executable(minseq minseq.cpp min_interval.cpp lcp_interval.cpp cigar.cpp)

target_link_libraries(minseq minseq_library common sequence m)