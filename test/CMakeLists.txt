cmake_minimum_required(VERSION 2.6)
project(test)

include_directories(../libbiosim/)

set(SOURCES_UNIT_TESTER
  che/assembly.cpp
  che/cc.cpp
  che/cchb_dssp.cpp
  che/molecule.cpp
  che/sequence_interval.cpp
  che/ss.cpp
  che/io/file_assembly.cpp
  che/io/file_fasta.cpp
  che/io/file_psipredv.cpp
  che/io/file_sse_pool.cpp
  che/score/cm.cpp
  che/score/cm_cc_blosum.cpp
  math/interval.cpp
  math/tensor.cpp
  math/algo/dp.cpp
  math/algo/interval_scheduler_maximize.cpp
  tools/enumerate.cpp
  tools/file.cpp
  tools/incrementor.cpp
  tools/string.cpp
  unit_tester.cpp
)

add_executable(unit_tester ${SOURCES_UNIT_TESTER})
target_link_libraries(unit_tester biosim)
