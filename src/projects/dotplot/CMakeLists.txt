project(dotplot)

add_library(dotplot_library STATIC dotplot.cpp)

add_executable(dotplot dotplot.cpp)

target_link_libraries(dotplot dotplot_library common sequence m)