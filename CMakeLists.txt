cmake_minimum_required(VERSION 2.6)

list(APPEND CMAKE_CXX_FLAGS "-std=c++11")

add_subdirectory(libbiosim)
add_subdirectory(test)
add_subdirectory(src)
