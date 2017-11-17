# Build time_expr C++ library
add_library(time_expr_c time_expr.cpp)
target_link_libraries(time_expr_c PUBLIC triqs)

# Build main solver library
add_library(realevol_c solver.cpp init_state.cpp make_init_state.cpp
            propagator.cpp utility.cpp
            triqs/operators/many_body_operator.cpp
            triqs/hilbert_space/fundamental_operator_set.cpp)
target_link_libraries(realevol_c PUBLIC time_expr_c triqs ${ARPACK_LIB})
add_dependencies(realevol_c time_expr_c)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TRIQS_INCLUDE_ALL})

# Install library
install(TARGETS realevol_c time_expr_c DESTINATION lib EXPORT realevol-targets)
install(EXPORT realevol-targets DESTINATION share/cmake)

# Install headers
set(HEADERS_DIR include/realevol)
install(FILES common.hpp time_expr.hpp solver.hpp utility.hpp
              compute_2t_obs_parameters.hpp init_state.hpp
        DESTINATION ${HEADERS_DIR})
install(DIRECTORY exprtk DESTINATION ${HEADERS_DIR}/exprtk FILES_MATCHING PATTERN *.hpp)
install(DIRECTORY triqs/operators DESTINATION ${HEADERS_DIR}/triqs FILES_MATCHING PATTERN *.hpp)
install(DIRECTORY triqs/hilbert_space DESTINATION ${HEADERS_DIR}/triqs)
