# Version number of the application
set(REALEVOL_VERSION "0.6")

# Start configuration
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(realevol CXX)
set(CMAKE_CXX_STANDARD 14)

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Load TRIQS and Cpp2Py
find_package(TRIQS 1.5 EXACT REQUIRED)
find_package(Cpp2Py 1.5 REQUIRED)

# Append TRIQS installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_ROOT}/share/cmake)

# Build the optimized version by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# Default install directory to TRIQS_ROOT if not given
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
endif()
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")

# Detect MPI
find_package(MPI)

# Detect Boost
find_package(Boost 1.53 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Detect ARPACK-NG
find_package(ARPACK REQUIRED)

# Get git hash
triqs_get_git_hash_of_source_dir(REALEVOL_GIT_HASH)
message(STATUS "Git hash: ${REALEVOL_GIT_HASH}")

option(Tests "Enable Tests" ON)

triqs_set_rpath_for_all_targets()

# Compile C++ code
add_subdirectory(c++)

include_directories(c++)

# Python interface
if(${TRIQS_WITH_PYTHON_SUPPORT})
  add_subdirectory(python)
endif()

# Testing
if(${Tests})
  enable_testing()
  add_subdirectory(test)
endif()

# Compile C++ benchmarks
#add_subdirectory(benchmark)
