# start configuration
cmake_minimum_required(VERSION 2.8)

# Version number of the application
project(realevol CXX)
set(REALEVOL_VERSION "0.5")

# Build the optimized version by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

if(NOT TRIQS_VERSION EQUAL 1.4)
    message(FATAL_ERROR "The application requires the TRIQS library version 1.4 (got ${TRIQS_VERSION})")
endif()

# Detect Boost
find_package(Boost 1.53 REQUIRED COMPONENTS serialization)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

# Get hash
triqs_get_git_hash(${CMAKE_SOURCE_DIR} "REALEVOL")
if(${GIT_RESULT} EQUAL 0)
    message(STATUS "Hash: ${REALEVOL_GIT_HASH}")
endif(${GIT_RESULT} EQUAL 0)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

message(STATUS "TRIQS : Adding compilation flags detected by the library (C++11/14, libc++, etc...) ")
add_definitions(${TRIQS_CXX_DEFINITIONS})

option(Tests "Enable Tests" ON)

# Compile C++ code
add_subdirectory(c++)
include_directories(c++)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
    add_subdirectory(python)
    if (${Tests})
        enable_testing()
        add_subdirectory(test)
    endif()
# option(BUILD_DOC "Build documentation" OFF)
# if(${BUILD_DOC})
#   if(NOT ${TRIQS_WITH_DOCUMENTATION})
#     message(FATAL_ERROR "Error: TRIQS library has not been compiled with its documentation")
#   endif()
#   add_subdirectory(doc)
# endif(${BUILD_DOC})
endif()

# Compile C++ benchmarks
#add_subdirectory(benchmark)
