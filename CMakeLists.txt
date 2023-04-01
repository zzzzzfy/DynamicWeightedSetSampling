cmake_minimum_required(VERSION 3.22)
project(DynamicSetSampling)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
add_executable(DynamicSetSampling main.cpp ScapegoatTree.hpp BSTSampling.hpp ChunkScapegoatTree.hpp BF.hpp)
