# Build time_expr C++ library
add_library(time_expr_c time_expr.cpp)
target_link_libraries(time_expr_c ${TRIQS_LIBRARY_ALL})

# Build time_interp C++ library
add_library(time_interp_c time_interp.cpp)
target_link_libraries(time_interp_c ${TRIQS_LIBRARY_ALL})

# Build main solver library
add_library(realevol_c solver.cpp init_state.cpp make_init_state.cpp
            propagator.cpp utility.cpp
            triqs/operators/many_body_operator.cpp
            triqs/hilbert_space/fundamental_operator_set.cpp)
target_link_libraries(realevol_c time_expr_c time_interp_c ${TRIQS_LIBRARY_ALL} arpack)
triqs_set_rpath_for_target(realevol_c)
add_dependencies(realevol_c time_expr_c time_interp_c)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${TRIQS_INCLUDE_ALL})

# Install library
install(TARGETS realevol_c time_expr_c time_interp_c DESTINATION lib)

# Install headers
set(HEADERS include/realevol)
install(FILES common.hpp time_expr.hpp time_interp.hpp solver.hpp utility.hpp
              compute_2t_obs_parameters.hpp init_state.hpp
        DESTINATION ${HEADERS})
install(DIRECTORY exprtk DESTINATION ${HEADERS}/exprtk FILES_MATCHING PATTERN *.hpp)
install(DIRECTORY triqs/operators DESTINATION ${HEADERS}/triqs FILES_MATCHING PATTERN *.hpp)
install(DIRECTORY triqs/hilbert_space DESTINATION ${HEADERS}/triqs)
