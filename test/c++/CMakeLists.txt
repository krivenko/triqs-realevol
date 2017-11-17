# Copy h5 files with input data
file(GLOB h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${h5_files} DESTINATION .)

# Simple tests
set(SIMPLE_TESTS time_expr many_body_operator hilbert_space space_partition
                 mesh_container lanczos_worker hs_structure sort2
                 init_state_pure
                 trimer)
foreach(t ${SIMPLE_TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    target_link_libraries(${t} realevol_c time_expr_c triqs gtest)
    add_test(${t} ${t})
endforeach()

# Parallel tests
foreach(t mpi_dispatcher init_state_equilibrium)
  add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
  target_link_libraries(${t} realevol_c triqs gtest)

  foreach(np 1 2 4)
    set(filename)
    add_test(NAME "${t}_${np}"
             COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${t})
  endforeach()
endforeach()
