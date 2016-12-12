# make a local pytriqs copy
triqs_prepare_local_pytriqs(${CMAKE_PROJECT_NAME})

# Build the python modules
triqs_python_extension_and_install(texpr)
target_link_libraries(texpr realevol_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(texpr)

triqs_python_extension_and_install(operators)
target_link_libraries(operators realevol_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(operators)

triqs_python_extension_and_install(gf_retime)
target_link_libraries(gf_retime realevol_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(gf_retime)

triqs_python_extension_and_install(init_state)
target_link_libraries(init_state realevol_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(init_state)

triqs_python_extension_and_install(realevol)
target_link_libraries(realevol realevol_c ${TRIQS_LIBRARY_ALL})
triqs_set_rpath_for_target(realevol)

include_directories(${CMAKE_BINARY_DIR}/include ${CMAKE_SOURCE_DIR}/c++ ${CMAKE_CURRENT_SOURCE_DIR} ${TRIQS_INCLUDE_ALL})

file(COPY _gf_retime.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
install(FILES _gf_retime.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${CMAKE_PROJECT_NAME})

# add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${CMAKE_PROJECT_NAME})
