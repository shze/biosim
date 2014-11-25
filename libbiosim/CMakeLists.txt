cmake_minimum_required(VERSION 2.6)
project(libbiosim)

# include current dir to avoid CMakeLists.txt in every subdir and a add_subdirectory(...) command for every subdir
include_directories(${CMAKE_CURRENT_LIST_DIR}) 

set(LIBBIOSIM_SOURCES
  che/cc.cpp
  che/cchb_dssp.cpp
  che/qs.cpp
  che/ss.cpp
  che/ts.cpp
  che/io/file_fasta.cpp
  che/io/file_psipredv.cpp
  che/io/file_qs.cpp
  che/io/file_sse_pool.cpp
  tools/file.cpp
  tools/incrementor.cpp
  tools/log.cpp
  tools/string.cpp
)

add_library(biosim STATIC ${LIBBIOSIM_SOURCES})

find_package(Boost REQUIRED COMPONENTS regex system filesystem)
target_link_libraries(biosim ${Boost_LIBRARIES})
