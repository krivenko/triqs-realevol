# make a local pytriqs copy
triqs_prepare_local_pytriqs(${CMAKE_PROJECT_NAME})

# Build the python modules
triqs_python_extension_and_install(texpr)
target_link_libraries(texpr realevol_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(texpr)

triqs_python_extension_and_install(operators)
target_link_libraries(operators realevol_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(operators)

# add version file
configure_file(version.py.in version.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${CMAKE_PROJECT_NAME})
