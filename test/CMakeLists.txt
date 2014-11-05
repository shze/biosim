cmake_minimum_required(VERSION 2.6)
project(test)

include_directories(../libbiosim/)

set(SOURCES_UNIT_TESTER
  che/cc.cpp
  che/cchb_dssp.cpp
  che/qs.cpp
  che/sequence_interval.cpp
  che/ss.cpp
  che/ts.cpp
  che/io/file_fasta.cpp
  che/io/file_psipredv.cpp
  che/io/file_qs.cpp
  che/io/file_sse_pool.cpp
  math/interval.cpp
  math/interval_scheduler_maximize.cpp
  tools/enumerate.cpp
  tools/file.cpp
  tools/incrementor.cpp
  unit_tester.cpp
)

add_executable(unit_tester ${SOURCES_UNIT_TESTER})
target_link_libraries(unit_tester biosim)
