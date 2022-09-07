cmake_minimum_required(VERSION 3.15)
project(.)

set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-lstdc++fs -ggdb3" )
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_SHARED_LINKER_FLAGS "-static-libgcc static-libstdc++ -Wall -Wc++-compat -O2 -msse4.1 -DHAVE_KALLOC -DKSW_CPU_DISPATCH -D_FILE_OFFSET_BITS=64 -ltbb -fsigned-char -fsanitize=address")

# Define various dirs
set(MAIN_SRC_DIR src)
set(LIB_SRC_DIR ${MAIN_SRC_DIR}/lib)
set(PROJECTS_SRC_DIR ${MAIN_SRC_DIR}/projects)
set(TESTS_SRC_DIR ${MAIN_SRC_DIR}/tests)
set(TOOLS_SRC_DIR ${MAIN_SRC_DIR}/tools)

set(TA_SRC_DIR ${PROJECTS_SRC_DIR}/tandem_aligner)

include_directories(${TOOLS_SRC_DIR})
include_directories(${TESTS_SRC_DIR})
include_directories(${PROJECT_SRC_DIR})

add_subdirectory(src/projects)
add_subdirectory(src/tools)
