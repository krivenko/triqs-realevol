# Copy h5 files with input data
file(GLOB_RECURSE h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
foreach(f ${h5_files})
  configure_file(${f} ${f} COPYONLY)
endforeach()

# Add a test
macro(add_realevol_test test)
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_executable(${test_name} ${test})
  target_link_libraries(${test_name}
                        ${PROJECT_NAME}::realevol_time_expr_c
                        ${PROJECT_NAME}::realevol_time_interp_c
                        ${PROJECT_NAME}::realevol_c
                        ${PROJECT_NAME}_warnings
                        gtest_main)
  set_property(TARGET ${test_name} PROPERTY
               RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})

  if(${ARGC} GREATER 1)
    foreach(numproc ${ARGN})
      add_test(NAME ${test_name}_np${numproc}
               COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${numproc}
                       ${MPIEXEC_PREFLAGS} ${test_name} ${MPIEXEC_POSTFLAGS}
               WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
    endforeach()
  else()
    add_test(NAME ${test_name} COMMAND ${test_name}
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  endif()

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    set_target_properties(${test_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    add_custom_command(
      TARGET ${test_name}
      COMMAND ${CPPCHECK_EXECUTABLE}
      --enable=warning,style,performance,portability
      --std=c++17
      --template=gcc
      --verbose
      --force
      --quiet
      ${CMAKE_CURRENT_SOURCE_DIR}/${test}
    )
  endif()
endmacro()

# Simple tests
set(SIMPLE_TESTS time_expr time_interp many_body_operator hilbert_space space_partition
                 mesh_container lanczos_worker hs_structure sort2 interpolator
                 init_state_pure
                 trimer)
foreach(test ${SIMPLE_TESTS})
  add_realevol_test(${test})
endforeach()

# MPI tests
set(MPI_TESTS mpi_dispatcher init_state_equilibrium)
foreach(test ${MPI_TESTS})
  add_realevol_test(${test} 1 2 4)
endforeach()
