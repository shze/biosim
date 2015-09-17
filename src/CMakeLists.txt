cmake_minimum_required(VERSION 2.6)
project(src)
include_directories(../libbiosim/)

find_package(Boost 1.55 REQUIRED COMPONENTS program_options) # 1.41 is too old
include_directories(${Boost_INCLUDE_DIRS})

set(SOURCES_ALIGN align.cpp)
add_executable(align ${SOURCES_ALIGN})
target_link_libraries(align biosim ${Boost_LIBRARIES})

set(SOURCES_SSSTAT ssstat.cpp)
add_executable(ssstat ${SOURCES_SSSTAT})
target_link_libraries(ssstat biosim ${Boost_LIBRARIES})
