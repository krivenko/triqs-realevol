# ##############################################################################
#
# realevol - Real time evolution solver based on TRIQS
#
# Copyright (C) 2014-2023, I. Krivenko, M. Danilov, P. Kubiczek
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

# Copy h5 files with input data
file(GLOB_RECURSE h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
foreach(f ${h5_files})
  configure_file(${f} ${f} COPYONLY)
endforeach()

# Add a test
macro(add_realevol_cpp_test test)
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_executable(${test_name} ${test})
  target_link_libraries(${test_name}
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
    set_target_properties(${test_name} PROPERTIES
                          CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
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
set(SIMPLE_TESTS time_expr.cpp
                 time_interp.cpp
                 many_body_operator.cpp
                 hilbert_space.cpp
                 space_partition.cpp
                 lanczos_worker.cpp
                 hs_structure.cpp
                 sort2.cpp
                 interpolator.cpp
                 init_state_pure.cpp)
foreach(test ${SIMPLE_TESTS})
  add_realevol_cpp_test(${test})
endforeach()

# MPI tests
set(MPI_TESTS mpi_dispatcher.cpp init_state_equilibrium.cpp)
foreach(test ${MPI_TESTS})
  add_realevol_cpp_test(${test} 1 2 4)
endforeach()
add_realevol_cpp_test(worldlines.cpp 1)

