project(suffix_array)

add_library(suffix_array_library STATIC suffix_array.cpp)

add_executable(suffix_array suffix_array.cpp)

target_link_libraries(suffix_array suffix_array_library common sequence m)