cmake_minimum_required(VERSION 2.6)
project(test)

include_directories(../libbiosim/)

set(SOURCES_UNIT_TESTER
  che/assembly.cpp
  che/assembly_atom.cpp
  che/assembly_cc.cpp
  che/atom.cpp
  che/bond.cpp
  che/cc.cpp
  che/cchb_dssp.cpp
  che/sequence_interval.cpp
  che/ss.cpp
  che/structure.cpp
  che/structure_ensemble.cpp
  che/algo/aligner_dp.cpp
  che/io/file_assembly.cpp
  che/io/file_fasta.cpp
  che/io/file_pdb.cpp
  che/io/file_psipredv.cpp
  che/io/file_pssm.cpp
  che/io/file_sse_pool.cpp
  che/score/cm.cpp
  che/score/cm_cc_blosum.cpp
  che/score/cm_cc_identity.cpp
  che/score/ev_alignment.cpp
  math/floating_point.cpp
  math/interval.cpp
  math/point.cpp
  math/point_typed.cpp
  math/tensor.cpp
  math/algo/dp.cpp
  math/algo/factorial.cpp
  math/algo/interval_scheduler_maximize.cpp
  math/algo/point_distance.cpp
  tools/enumerate.cpp
  tools/file.cpp
  tools/mapper.cpp
  tools/string.cpp
  unit_tester.cpp
)

add_executable(unit_tester ${SOURCES_UNIT_TESTER})
target_link_libraries(unit_tester biosim)
