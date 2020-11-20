# Copy h5 files with input data
file(GLOB_RECURSE h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
foreach(f ${h5_files})
  configure_file(${f} ${f} COPYONLY)
endforeach()

# Add a test
macro(add_realevol_python_test test)
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)

  if(${ARGC} GREATER 1)
    foreach(numproc ${ARGN})
      add_test(NAME Py_${test_name}_np${numproc}
               COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${numproc}
                       ${MPIEXEC_PREFLAGS}
                       ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py
                       ${MPIEXEC_POSTFLAGS}
               WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
    endforeach()
  else()
    add_test(NAME Py_${test_name}
             COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  endif()

  set_property(TEST Py_${test_name} APPEND PROPERTY
               ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH}
               ${SANITIZER_RT_PRELOAD})
endmacro()

# Simple tests
set(SIMPLE_TESTS texpr tinterp)# operators_texpr operators_tinterp) #gf_retime
foreach(test ${SIMPLE_TESTS})
  add_realevol_python_test(${test})
endforeach()

# TODO
# MPI tests
#set(MPI_TESTS star_texpr star_tinterp)
#foreach(test ${MPI_TESTS})
#  add_realevol_python_test(${test} 1 2 4)
#endforeach()
