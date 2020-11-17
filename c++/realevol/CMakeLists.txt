# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()


# ========= Dynamic Analyzer Checks ==========

option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")

if(ASAN)
  if(NOT TARGET asan)
    find_package(sanitizer REQUIRED "asan")
  endif()
endif()
if(UBSAN)
  if(NOT TARGET ubsan)
    find_package(sanitizer REQUIRED "ubsan")
  endif()
endif()

# ========= Build and install C++ libraries =========
macro(add_realevol_cpp_lib CPPLIBNAME)

  set(sources ${ARGN})

  # Add library target
  add_library(${CPPLIBNAME} ${sources})
  add_library(${PROJECT_NAME}::${CPPLIBNAME} ALIAS ${CPPLIBNAME})

  # Link against triqs and enable warnings
  target_link_libraries(${CPPLIBNAME} PUBLIC triqs
                        PRIVATE $<BUILD_INTERFACE:${PROJECT_NAME}_warnings>)

  # Configure target and compilation
  set_property(TARGET ${CPPLIBNAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
  target_include_directories(${CPPLIBNAME} PUBLIC
                             $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
  target_include_directories(${CPPLIBNAME} SYSTEM INTERFACE
                             $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
  target_compile_definitions(${CPPLIBNAME} PUBLIC
                             REALEVOL_GIT_HASH=${PROJECT_GIT_HASH}
                             TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
                             $<$<CONFIG:Debug>:REALEVOL_DEBUG>
                             $<$<CONFIG:Debug>:TRIQS_DEBUG>
                             $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
  )

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    set_target_properties(${CPPLIBNAME} PROPERTIES
                          CXX_CLANG_TIDY
                          "${CLANG_TIDY_EXECUTABLE};-header-filter=^.*/realevol/[a-z0-9_]+\.hpp$")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    add_custom_command(
      TARGET ${CPPLIBNAME}
      COMMAND ${CPPCHECK_EXECUTABLE}
        --enable=warning,style,performance,portability
        --std=c++17
        --template=gcc
        --verbose
        --force
        --quiet
        ${sources}
      WORKING_DIRECTORY
      ${CMAKE_CURRENT_SOURCE_DIR}
    )
  endif()

  # Runtime analyzers
  if(ASAN)
    target_link_libraries(${CPPLIBNAME} PUBLIC $<BUILD_INTERFACE:asan>)
  endif()
  if(UBSAN)
    target_link_libraries(${CPPLIBNAME} PUBLIC $<BUILD_INTERFACE:ubsan>)
  endif()

  # Install library
  install(TARGETS ${CPPLIBNAME} EXPORT ${PROJECT_NAME}-targets DESTINATION lib)
endmacro()

add_realevol_cpp_lib(realevol_time_expr_c time_expr.cpp)
add_realevol_cpp_lib(realevol_time_interp_c time_interp.cpp)

# Build main solver library
add_realevol_cpp_lib(realevol_c solver.cpp init_state.cpp make_init_state.cpp
                     propagator.cpp utility.cpp
                     operators/many_body_operator.cpp
                     hilbert_space/fundamental_operator_set.cpp)

target_link_libraries(realevol_c PUBLIC time_expr_c time_interp_c)
target_link_libraries(realevol_c PRIVATE ezarpack)
add_dependencies(realevol_c time_expr_c time_interp_c)

# Install headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")
