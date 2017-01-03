find_package(TriqsTest)

set(PYTHON_TESTS texpr operators gf_retime)
foreach(t ${PYTHON_TESTS})
  triqs_add_python_test(${t})
endforeach(t)

# Parallel tests
foreach(t star)

 foreach(np 1 2 4)
  set(TEST_MPI_NUMPROC ${np})
  triqs_add_python_test(${t})
 endforeach()
 unset(TEST_MPI_NUMPROC)

endforeach()

# Copy h5 files with input data
file(GLOB h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${h5_files} DESTINATION .)
