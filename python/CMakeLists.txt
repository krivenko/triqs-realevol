set(PYTHON_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
)

# where will the python end up in triqs?
set(python_destination applications/impurity_solvers/realevol)
# make a local pytriqs copy
triqs_prepare_local_pytriqs_merged_with_my_python(${python_destination})

include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_SOURCE_DIR})

set(all_modules texpr)
foreach(module ${all_modules})
  triqs_python_extension(${module} ${python_destination})
  target_link_libraries(${module} realevol_c)
  triqs_set_rpath_for_target(${module})
  execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/${module}.so ${CMAKE_BINARY_DIR}/pytriqs/${python_destination})
endforeach()

configure_file(version.py.in version.py)
install (FILES ${CMAKE_BINARY_DIR}/pytriqs/${python_destination}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST}/${python_destination})
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/${python_destination})
