include(external_dependency.cmake)

if(NOT DEFINED Build_Deps)
  set(Build_Deps "IfNotFound" CACHE STRING "Do we build dependencies from source? [Never/Always/IfNotFound]")
else()
  set(Build_Deps_Opts "Never" "Always" "IfNotFound")
  if(NOT ${Build_Deps} IN_LIST Build_Deps_Opts)
    message(FATAL_ERROR "Build_Deps option should be either 'Never', 'Always' or 'IfNotFound'")
  endif()
  set(Build_Deps ${Build_Deps} CACHE STRING "Do we build dependencies from source? [Never/Always/IfNotFound]")
  if(NOT IS_SUBPROJECT AND NOT Build_Deps STREQUAL "Always" AND (ASAN OR UBSAN))
    message(WARNING "For builds with llvm sanitizers (ASAN/UBSAN) it is recommended to use -DBuild_Deps=Always to avoid false positives.")
  endif()
endif()

# -- Cpp2Py --
if(PythonSupport OR Build_Documentation)
  external_dependency(Cpp2Py
    GIT_REPO https://github.com/TRIQS/cpp2py
    VERSION 2.0
    GIT_TAG master
    BUILD_ALWAYS
    EXCLUDE_FROM_ALL
  )
endif()

# -- GoogleTest --
external_dependency(GTest
  GIT_REPO https://github.com/google/googletest
  GIT_TAG release-1.10.0
  BUILD_ALWAYS
  EXCLUDE_FROM_ALL
)

# -- Boost --
find_package(Boost 1.53 REQUIRED)

# -- ezARPACK + ARPACK-NG --
set(Examples OFF CACHE BOOL "")
set(Tests OFF CACHE BOOL "")
external_dependency(ezARPACK
  VERSION 0.9
  GIT_REPO https://github.com/krivenko/ezARPACK
  GIT_TAG master
  BUILD_ALWAYS
  EXCLUDE_FROM_ALL
)

# -- Google Benchmark --
if(Build_Benchmarks)
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "")
  set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "")
  external_dependency(GBenchmark
    GIT_REPO https://github.com/google/benchmark
    GIT_TAG v1.5.0
    BUILD_ALWAYS
    EXCLUDE_FROM_ALL
  )
endif()
