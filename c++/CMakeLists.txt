add_library(time_expr_c time_expr.cpp)
target_link_libraries(time_expr_c ${TRIQS_LIBRARY_ALL})

add_library(realevol_c solver.cpp init_state.cpp triqs/operators/many_body_operator.cpp)
target_link_libraries(realevol_c time_expr_c ${TRIQS_LIBRARY_ALL} ${ARPACK_LDFLAGS})
triqs_set_rpath_for_target(realevol_c)
add_dependencies(realevol_c time_expr_c)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TRIQS_INCLUDE_ALL})

install(TARGETS realevol_c time_expr_c DESTINATION lib)
