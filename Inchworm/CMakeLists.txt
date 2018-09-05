cmake_minimum_required (VERSION 3.1)
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH})
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH})

project (Inchworm)

enable_language(CXX)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "${Red}No build type specified so building in RELEASE mode.${ColourReset}")
    set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe -W -Wall -Wpedantic -fopenmp -m64")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm -ldl ")

MESSAGE(STATUS "system: ${CMAKE_SYSTEM_NAME}")
if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lrt")
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++")
endif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")


set(inchworm_SOURCES src/Fasta_entry.cpp
                     src/IRKE_run.cpp
                     src/sequenceUtil.cpp
                     src/IRKE.cpp
                     src/KmerCounter.cpp
                     src/string_util.cpp 
                     src/Fasta_reader.cpp
                     src/IRKE.hpp
                     src/stacktrace.cpp
                     src/argProcessor.cpp)
add_executable(inchworm ${inchworm_SOURCES})
install(TARGETS inchworm DESTINATION bin)


set(FastaToDeBruijn_SOURCES src/FastaToDeBruijn.cpp
                            src/argProcessor.cpp
                            src/Fasta_reader.cpp
                            src/Fasta_entry.cpp
                            src/sequenceUtil.cpp
                            src/string_util.cpp
                            src/stacktrace.cpp
                            src/DeBruijnGraph.cpp)
add_executable(FastaToDeBruijn ${FastaToDeBruijn_SOURCES})
install(TARGETS FastaToDeBruijn DESTINATION bin)



set(fastaToKmerCoverageStats_SOURCES src/fastaToKmerCoverageStats.cpp
                                     src/argProcessor.cpp
                                     src/Fasta_reader.cpp
                                     src/Fasta_entry.cpp
                                     src/sequenceUtil.cpp
                                     src/string_util.cpp
                                     src/stacktrace.cpp
                                     src/KmerCounter.cpp)
add_executable(fastaToKmerCoverageStats ${fastaToKmerCoverageStats_SOURCES})
install(TARGETS fastaToKmerCoverageStats DESTINATION bin)




