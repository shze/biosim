cmake_minimum_required(VERSION 2.6)
project(libbiosim)

# include current dir to avoid CMakeLists.txt in every subdir and a add_subdirectory(...) command for every subdir
include_directories(${CMAKE_CURRENT_LIST_DIR}) 

set(LIBBIOSIM_SOURCES
  che/assembly.cpp
  che/cc.cpp
  che/cchb_dssp.cpp
  che/molecule.cpp
  che/ss.cpp
  che/io/file_assembly.cpp
  che/io/file_fasta.cpp
  che/io/file_pdb.cpp
  che/io/file_psipredv.cpp
  che/io/file_sse_pool.cpp
  score/cm.cpp
  score/cm_assembly_ssq3.cpp
  score/cm_cc_blosum.cpp
  tools/file.cpp
  tools/incrementor.cpp
  tools/log.cpp
  tools/string.cpp
)

add_library(biosim STATIC ${LIBBIOSIM_SOURCES})

find_package(Boost REQUIRED COMPONENTS regex system filesystem)
target_link_libraries(biosim ${Boost_LIBRARIES})
