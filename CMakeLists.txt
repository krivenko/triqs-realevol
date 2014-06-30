# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# start configuration 
cmake_minimum_required(VERSION 2.8)
project(realevol CXX)
set(CMAKE_BUILD_TYPE Release)

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Detect muParser
find_package(PkgConfig)
pkg_check_modules(MUPARSER REQUIRED muparser)
add_definitions(${MUPARSER_CFLAGS})
link_libraries(${MUPARSER_LDFLAGS})
message(STATUS "muParser version: " ${MUPARSER_VERSION})

# Compile C++ code
add_subdirectory(c++)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
 if(${Tests})
  enable_testing()
  add_subdirectory(test)
 endif(${Tests})
endif(${TRIQS_WITH_PYTHON_SUPPORT})