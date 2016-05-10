add_library(realevol_c solver.cpp init_state.cpp time_expr.cpp triqs/operators/many_body_operator.cpp)
target_link_libraries(realevol_c ${TRIQS_LIBRARY_ALL})
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TRIQS_INCLUDE_ALL})
triqs_set_rpath_for_target(realevol_c)

install(TARGETS realevol_c DESTINATION lib)
