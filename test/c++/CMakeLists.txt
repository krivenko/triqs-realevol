# load triqs helper to set up tests
find_package(TriqsTest)

# Simple tests
set(SIMPLE_TESTS time_expr many_body_operator hilbert_space space_partition
                 mesh_container lanczos_worker hs_structure init_state trimer) #mesh_container_cyclic runge_kutta lanczos)
foreach(t ${SIMPLE_TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    target_link_libraries(${t} realevol_c time_expr_c ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST})
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach()

# Parallel tests
add_executable(mpi_dispatcher mpi_dispatcher.cpp)
target_link_libraries(mpi_dispatcher ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST})
triqs_set_rpath_for_target(mpi_dispatcher)

set(TEST_MPI_NUMPROC 1)
triqs_add_cpp_test(mpi_dispatcher)
set(TEST_MPI_NUMPROC 2)
triqs_add_cpp_test(mpi_dispatcher)
set(TEST_MPI_NUMPROC 4)
triqs_add_cpp_test(mpi_dispatcher)
unset(TEST_MPI_NUMPROC)

add_subdirectory(arpack)
