# load triqs helper to set up tests
find_package(TriqsTest)

link_libraries(realevol_c ${TRIQS_LIBRARY_ALL} gtest)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../c++ ${TRIQS_INCLUDE_ALL})

# Simple tests
set(SIMPLE_TESTS time_expr many_body_operator hilbert_space space_partition lpt_scheduler
                 mesh_container lanczos_worker hs_structure init_state trimer) #mesh_container_cyclic runge_kutta lanczos)
foreach(t ${SIMPLE_TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach(t)
