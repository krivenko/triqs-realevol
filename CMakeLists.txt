# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
#
# realevol is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# realevol is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# realevol (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
# ##############################################################################

cmake_minimum_required(VERSION 3.3.2 FATAL_ERROR)
cmake_policy(VERSION 3.3.2)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()
if(POLICY CMP0077)
  cmake_policy(SET CMP0077 NEW)
endif()

# ############
# Define Project
project(triqs_realevol VERSION 0.6 LANGUAGES CXX)
get_directory_property(IS_SUBPROJECT PARENT_DIRECTORY)

# ############
# Load TRIQS and CPP2PY
find_package(TRIQS 3.0 REQUIRED)

# Get the git hash & print status
triqs_get_git_hash_of_source_dir(PROJECT_GIT_HASH)
message(STATUS "${PROJECT_NAME} version : ${PROJECT_VERSION}")
message(STATUS "${PROJECT_NAME} Git hash: ${PROJECT_GIT_HASH}")

# Default Install directory to TRIQS_ROOT if not given or invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
  set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)
endif()
if(NOT IS_SUBPROJECT)
  message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")
endif()
set(${PROJECT_NAME}_BINARY_DIR ${PROJECT_BINARY_DIR} CACHE
    STRING "Binary directory of the ${PROJECT_NAME} Project")

# ############
# Options

# Make additional Find Modules available
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/share/cmake/Modules)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
if(NOT IS_SUBPROJECT)
  message(STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} --------")
endif()

# Python Support
option(PythonSupport "Build with Python support" ON)
if(PythonSupport AND NOT TRIQS_WITH_PYTHON_SUPPORT)
  message(FATAL_ERROR "TRIQS was installed without Python support. Cannot build the Python Interface. Disable the build with -DPythonSupport=OFF")
endif()

option(Use_Antihermiticity "Use the antihermiticity of GFs/susceptibilities to skip computation of some (t,t') points" ON)
mark_as_advanced(Use_Antihermiticity)
if(Use_Antihermiticity)
  add_definitions(-DUSE_ANTIHERMITICITY)
endif(Use_Antihermiticity)

# Compile C++ benchmarks
option(Build_Benchmarks "Compile C++ benchmarks (requires Google benchmark)" OFF)

# Documentation
# TODO
# option(Build_Documentation "Build documentation" OFF)

# Testing
option(Build_Tests "Build tests" ON)
if(Build_Tests)
  enable_testing()
endif()

# Export the list of compile-commands into compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Global compiler options
option(BUILD_SHARED_LIBS "Enable compilation of shared libraries" OFF)
add_compile_options($<$<CONFIG:Debug>:-ggdb3>)

# Create an Interface target for compiler warnings
add_library(${PROJECT_NAME}_warnings INTERFACE)
target_compile_options(${PROJECT_NAME}_warnings
  INTERFACE
    -Wall
    -Wextra
    -Wpedantic
    -Wno-sign-compare
    $<$<CXX_COMPILER_ID:GNU>:-Wno-comma-subscript>
    $<$<CXX_COMPILER_ID:GNU>:-Wshadow=local>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-attributes>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-deprecated-comma-subscript>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-unknown-warning-option>
    $<$<CXX_COMPILER_ID:Clang>:-Wshadow>
    $<$<CXX_COMPILER_ID:Clang>:-Wno-gcc-compat>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wno-deprecated-comma-subscript>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wno-unknown-warning-option>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wshadow>
    $<$<CXX_COMPILER_ID:AppleClang>:-Wno-gcc-compat>
)

# #############
# Build Project

# Find / Build dependencies
add_subdirectory(deps)

# Build and install the library
add_subdirectory(c++/${PROJECT_NAME})

# Tests
if(Build_Tests)
  add_subdirectory(test)
endif()

# Python
if(PythonSupport)
  add_subdirectory(python/${PROJECT_NAME})
endif()

# Docs
# TODO
#if(Build_Documentation)
#  add_subdirectory(doc)
#endif()

# Benchmarks
if(Build_Benchmarks)
  add_subdirectory(benchmark)
endif()

# Additional configuration files
add_subdirectory(share)

# #############
# Debian Package

option(BUILD_DEBIAN_PACKAGE "Build a deb package" OFF)
if(BUILD_DEBIAN_PACKAGE AND NOT IS_SUBPROJECT)
  if(NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr")
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be /usr for packaging")
  endif()
  set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
  set(CPACK_GENERATOR "DEB")
  set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
  set(CPACK_PACKAGE_CONTACT "https://github.com/krivenko/triqs-realevol")
  execute_process(COMMAND dpkg --print-architecture
                  OUTPUT_VARIABLE CMAKE_DEBIAN_PACKAGE_ARCHITECTURE
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "triqs (>= 3.0)")
  set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
  include(CPack)
endif()

# TODO
# Detect MPI
#find_package(MPI)

# Detect Boost
#find_package(Boost 1.53 REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})

# Detect ARPACK-NG
#find_package(ARPACK REQUIRED)


