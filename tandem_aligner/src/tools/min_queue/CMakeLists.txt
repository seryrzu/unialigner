project(min_queue)

add_library(min_queue_library STATIC min_queue.cpp)

add_executable(min_queue min_queue.cpp)

target_link_libraries(min_queue min_queue_library common sequence m)