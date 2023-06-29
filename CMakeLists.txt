cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
project(annoy_sratch VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

include_directories(./src)

add_executable(annoy_sratch main.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
