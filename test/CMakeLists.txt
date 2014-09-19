cmake_minimum_required(VERSION 2.6)
project(test)

include_directories(../libbiosim/)

set(SOURCES_UNIT_TESTER
  che/cc.cpp
  che/complex.cpp
  che/molecule.cpp
  che/io/file_fasta.cpp
  tools/enumerate.cpp
  tools/file.cpp
  tools/incrementor.cpp
  unit_tester.cpp
)

add_executable(unit_tester ${SOURCES_UNIT_TESTER})
target_link_libraries(unit_tester biosim)
