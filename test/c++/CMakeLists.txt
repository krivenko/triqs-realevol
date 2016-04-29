# load triqs helper to set up tests
find_package(TriqsTest)

link_libraries(realevol_c ${TRIQS_LIBRARY_ALL} gtest)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../c++ ${TRIQS_INCLUDE_ALL})

# Simple tests
set(SIMPLE_TESTS time_expr many_body_operator hilbert_space space_partition)
#set(SIMPLE_TESTS time_expr mesh_container mesh_container_cyclic hilbert operator
#                 runge_kutta lanczos_worker lanczos lpt_scheduler)
foreach(t ${SIMPLE_TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach(t)

# HDF5 tests
#set(HDF5_TESTS trimer_lanczos trimer_rk)
#foreach(t ${HDF5_TESTS})
#    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
#    triqs_set_rpath_for_target(${t})
#    triqs_add_test_C_hdf(${t} " -d 2.e-4")
#endforeach(t)
