project(sequences)
set(CMAKE_CXX_STANDARD 14)

include_directories(.)
add_library(sequence STATIC contigs.cpp sequence.cpp)

find_package (ZLIB)
target_link_libraries (sequence ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} m)
