cmake_minimum_required(VERSION 2.6)
project(src)
include_directories( ../libbiosim/)

set( SOURCES_SSSTAT
  ssstat.cpp
)

add_executable( ssstat ${SOURCES_SSSTAT})
target_link_libraries( ssstat biosim -lboost_program_options)


