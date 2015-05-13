cmake_minimum_required(VERSION 2.6)
project(src)
include_directories(../libbiosim/)

set(SOURCES_SSSTAT
  ssstat.cpp
)

add_executable(ssstat ${SOURCES_SSSTAT})

find_package(Boost REQUIRED COMPONENTS program_options)
target_link_libraries(ssstat biosim ${Boost_LIBRARIES})
