# load triqs helper to set up tests
find_package(TriqsTest)


include_directories(${CMAKE_SOURCE_DIR}/c++ ${TRIQS_INCLUDE_ALL})

# Simple tests
set(SIMPLE_TESTS time_expr many_body_operator hilbert_space space_partition
                 mesh_container lanczos_worker hs_structure sort2
                 init_state_pure
                 trimer) #mesh_container_cyclic runge_kutta lanczos)
foreach(t ${SIMPLE_TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    target_link_libraries(${t} realevol_c time_expr_c ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST})
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach()

# Parallel tests
foreach(t mpi_dispatcher init_state_equilibrium)
 add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
 target_link_libraries(${t} realevol_c ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST})
 triqs_set_rpath_for_target(${t})

 foreach(np 1 2 4)
  set(TEST_MPI_NUMPROC ${np})
  triqs_add_cpp_test(${t})
 endforeach()
 unset(TEST_MPI_NUMPROC)

endforeach()

# Copy h5 files with input data
file(GLOB h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${h5_files} DESTINATION .)

add_subdirectory(arpack)
