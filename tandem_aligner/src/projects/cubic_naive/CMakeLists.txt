project(cubic_naive)

add_library(cubic_naive_library STATIC cubic_naive.cpp cubic_naive.hpp)

add_executable(cubic_naive cubic_naive.cpp)

target_link_libraries(cubic_naive cubic_naive_library common sequence m)
