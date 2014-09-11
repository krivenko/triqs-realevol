set(PYTHON_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py
)

# where will the python end up in triqs?
set(python_destination applications/realevol)
# make a local pytriqs copy
triqs_prepare_local_pytriqs_merged_with_my_python(${python_destination})

foreach(ext texpr ctexpr meshes)
    triqs_python_extension(${ext} ${python_destination})
    target_link_libraries(${ext} realevol_c)
    include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_SOURCE_DIR})
    triqs_set_rpath_for_target(${ext})
endforeach()

add_dependencies(ctexpr texpr)

#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.py.in ${CMAKE_BINARY_DIR}/pytriqs/${python_destination}/version.py)
#install (FILES ${CMAKE_BINARY_DIR}/pytriqs/${python_destination}/version.py DESTINATION ${TRIQS_PYTHON_LIB_DEST}/${python_destination})
