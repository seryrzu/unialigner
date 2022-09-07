project(tandem_aligner)

add_library(tandem_aligner_library STATIC tandem_aligner.cpp)

add_executable(tandem_aligner tandem_aligner.cpp min_interval.cpp lcp_interval.cpp cigar.cpp)

target_link_libraries(tandem_aligner tandem_aligner_library common sequence m)