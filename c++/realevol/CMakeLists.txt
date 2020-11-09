# Build C++ libraries

foreach(NAME time_expr)

  set(CPPLIBNAME ${PROJECT_NAME}_${NAME}_c)

  # Add library target
  add_library(${CPPLIBNAME} time_expr.cpp)
  add_library(${PROJECT_NAME}::${CPPLIBNAME} ALIAS ${CPPLIBNAME})

  # Link against triqs and enable warnings
  target_link_libraries(${CPPLIBNAME} PUBLIC triqs
                        PRIVATE $<BUILD_INTERFACE:${PROJECT_NAME}_warnings>)

  # Configure target and compilation
  set_property(TARGET ${CPPLIBNAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
  target_include_directories(${CPPLIBNAME} PUBLIC
                             $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
  target_include_directories(${CPPLIBNAME} SYSTEM INTERFACE
                             $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
  target_compile_definitions(${CPPLIBNAME} PUBLIC
                             REALEVOL_GIT_HASH=${PROJECT_GIT_HASH}
                             TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
                             $<$<CONFIG:Debug>:REALEVOL_DEBUG>
                             $<$<CONFIG:Debug>:TRIQS_DEBUG>
                             $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
  )

  # Install library
  install(TARGETS ${CPPLIBNAME} EXPORT ${PROJECT_NAME}-targets DESTINATION lib)

endforeach()

# Install headers
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp")

# TODO

## Build time_interp C++ library
#add_library(time_interp_c time_interp.cpp)
#target_link_libraries(time_interp_c ${TRIQS_LIBRARY_ALL})

## Build main solver library
#add_library(realevol_c solver.cpp init_state.cpp make_init_state.cpp
            #propagator.cpp utility.cpp
            #triqs/operators/many_body_operator.cpp
            #triqs/hilbert_space/fundamental_operator_set.cpp)
#target_link_libraries(realevol_c time_expr_c time_interp_c ${TRIQS_LIBRARY_ALL} ${ARPACK_LIB})
#triqs_set_rpath_for_target(realevol_c)
#add_dependencies(realevol_c time_expr_c time_interp_c)
