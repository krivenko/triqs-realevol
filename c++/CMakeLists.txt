#add_library(realevol_c solver.cpp time_expr.cpp)
add_library(realevol_c time_expr.cpp)
target_link_libraries(realevol_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(realevol_c)

install(TARGETS realevol_c DESTINATION lib)
