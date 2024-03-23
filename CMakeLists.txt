cmake_minimum_required(VERSION 3.20.0)

project(Graph_BQ)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(FLAG "-Wall -Wextra -Wunused -Wshadow -Wnon-virtual-dtor -Wconversion -pedantic -Wpedantic -Wnull-dereference \
          -g3 -O3 -march=native -mtune=native")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")

set(SRCS
    lib/tinyxml2/tinyxml2.cpp
    src/configuration.cpp
    src/graphs.cpp
    src/graph_algorithms.cpp
    src/graph_viz.hpp
    )

set(HEADERS
    lib/tinyxml2/tinyxml2.h
    src/plf_nanotimer.h
    src/configuration.hpp
    src/graphs.hpp
    src/graph_algorithms.hpp
    src/graph_viz.cpp
    )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lcgraph -lgvc ")

add_executable(graph ${SRCS} src/main.cpp ${HEADERS})

add_executable(test_graph ${SRCS} src/test_graph.cpp ${HEADERS})



