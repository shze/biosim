cmake_minimum_required(VERSION 2.6)
project(libbiosim)

# include current dir to avoid CMakeLists.txt in every subdir and a add_subdirectory(...) command for every subdir
include_directories(${CMAKE_CURRENT_LIST_DIR}) 

find_package(Boost 1.55 REQUIRED COMPONENTS regex system filesystem)  # 1.41 is too old
include_directories(${Boost_INCLUDE_DIRS}) 

set(LIBBIOSIM_SOURCES
  che/alignment.cpp
  che/assembly.cpp
  che/assembly_atom.cpp
  che/assembly_cc.cpp
  che/atom.cpp
  che/atom_contact.cpp
  che/bond.cpp
  che/cc.cpp
  che/cchb_dssp.cpp
  che/ss.cpp
  che/structure.cpp
  che/structure_ensemble.cpp
  che/algo/aligner_dp.cpp
  che/io/file_assembly.cpp
  che/io/file_fasta.cpp
  che/io/file_pdb.cpp
  che/io/file_pdb_data.cpp
  che/io/file_psipredv.cpp
  che/io/file_pssm.cpp
  che/io/file_sse_pool.cpp
  che/score/cm.cpp
  che/score/cm_assembly_ssq3.cpp
  che/score/cm_cc_blosum.cpp
  che/score/cm_cc_identity.cpp
  che/score/ev_alignment.cpp
  math/point.cpp
  math/point_typed.cpp
  math/algo/factorial.cpp
  tools/file.cpp
  tools/log.cpp
  tools/string.cpp
)

add_library(biosim STATIC ${LIBBIOSIM_SOURCES})
target_link_libraries(biosim ${Boost_LIBRARIES})
