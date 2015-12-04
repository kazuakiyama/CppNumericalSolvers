#
# CMake for CppNumericalSolvers
#
cmake_minimum_required (VERSION 2.6.8)

# magic to find gfortran on Mac
find_program (FORTRAN_COMPILER_MAC "gfortran-mp-4.9" )
if (NOT FORTRAN_COMPILER_MAC)
find_program (FORTRAN_COMPILER_MAC "gfortran" )
endif (NOT FORTRAN_COMPILER_MAC)
if (FORTRAN_COMPILER_MAC)
set (CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER_MAC})
endif (FORTRAN_COMPILER_MAC)

# Default build mode
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE "RelWithDebInfo")
endif (NOT CMAKE_BUILD_TYPE)

#
# Project Name
#
project (CppNumericalSolvers CXX)
enable_language (C)

add_definitions(-DHAVE_ASA_CG)
include_directories (${CMAKE_SOURCE_DIR})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11 -Wall")

# local CMake module location
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

#
# Configuration options
#
option (ENABLE_OPENMP "Build with OpenMP support" ON)
if (ENABLE_OPENMP)
  find_package (OpenMP)
endif (ENABLE_OPENMP)
if (OPENMP_FOUND)
  add_definitions(-DHAVE_OPENMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else (OPENMP_FOUND)
  # no openmp: then dont warn about ignorning openmp pragmas
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif (OPENMP_FOUND)

#
# Packages
#
find_package (Eigen3 REQUIRED)
if (NOT EIGEN3_FOUND)
  message (FATAL_ERROR "Unable to find Eigen")
endif ()
include_directories (${EIGEN3_INCLUDE_DIRS})

find_package (IPOPT)
if (IPOPT_FOUND)
  include_directories (${IPOPT_INCLUDE_DIRS})
  add_definitions (-DHAVE_IPOPT)
endif (IPOPT_FOUND)

# gtest install
add_custom_command (OUTPUT "gtest-1.7.0/src/gtest-all.cc"
  COMMAND rm -f gtest-1.7.0.zip
  COMMAND rm -fR gtest-1.7.0
  COMMAND wget -O gtest-1.7.0.zip https://googletest.googlecode.com/files/gtest-1.7.0.zip
  COMMAND unzip gtest-1.7.0.zip
  COMMAND rm -f gtest-1.7.0.zip
  )
include_directories (${CMAKE_BINARY_DIR}/gtest-1.7.0/include)

# ASA_CG install
set (ASA_VER "ASA_CG")

# gtest
include_directories (${CMAKE_BINARY_DIR}/gtest-1.7.0)
add_library (gtest "gtest-1.7.0/src/gtest-all.cc")

#
# main
#
add_executable (main src/main.cpp ${ASA_VER}/asa_cg.c src/lbfgs.c)
target_link_libraries (main ${IPOPT_LIBRARIES})

#
# test
#
add_executable (test src/unittests.cpp ${ASA_VER}/asa_cg.c src/lbfgs.c)
target_link_libraries (test gtest)
target_link_libraries (test ${IPOPT_LIBRARIES})
