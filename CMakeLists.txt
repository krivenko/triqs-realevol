# start configuration
cmake_minimum_required(VERSION 2.8)

# Version number of the application
project(realevol CXX)
set(REALEVOL_VERSION "0.5")

# Build the optimized version by default
set(CMAKE_BUILD_TYPE Release)

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

if(NOT TRIQS_VERSION EQUAL 1.4)
  message(FATAL_ERROR "The application requires the TRIQS library version 1.4 (got ${TRIQS_VERSION})")
endif()

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Compile C++ code
add_subdirectory(c++)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)
 if(${Tests})
  enable_testing()
  add_subdirectory(test)
 endif(${Tests})
endif(${TRIQS_WITH_PYTHON_SUPPORT})

# Compile C++ benchmarks
add_subdirectory(benchmark)
