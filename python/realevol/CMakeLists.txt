# Configure the version
configure_file(version.py.in version.py)

# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
# FIXME
#file(GLOB_RECURSE wrap_generators RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *_desc.py)
file(GLOB_RECURSE wrap_generators RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
     texpr_desc.py tinterp_desc.py operators_texpr_desc.py operators_tinterp_desc.py init_state_desc.py)
list(REMOVE_ITEM python_sources "${wrap_generators}")
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install Python files to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}
        FILES_MATCHING PATTERN "*.py" PATTERN "*_desc.py" EXCLUDE)

# Build and install any Python modules
foreach(gen ${wrap_generators})
  string(REPLACE "_desc.py" "" gen ${gen})
  get_filename_component(module_name ${gen} NAME_WE)
  get_filename_component(module_dir ${gen} DIRECTORY)
  add_cpp2py_module(NAME ${module_name} DIRECTORY ${module_dir})
  target_link_libraries(${module_name} triqs_py)
  install(TARGETS ${module_name} DESTINATION ${PYTHON_LIB_DEST}/${module_dir})
endforeach()

# Link compiled Python modules to C++ libraries
target_link_libraries(texpr realevol_time_expr_c)
target_link_libraries(tinterp realevol_time_interp_c)
target_link_libraries(operators_texpr realevol_c)
target_link_libraries(operators_tinterp realevol_c)
target_link_libraries(init_state realevol_c)

#triqs_python_extension_and_install(gf_retime)
#target_link_libraries(gf_retime realevol_c ${TRIQS_LIBRARY_ALL})
#triqs_set_rpath_for_target(gf_retime)

#triqs_python_extension_and_install(realevol)
#target_link_libraries(realevol realevol_c ${TRIQS_LIBRARY_ALL})
#triqs_set_rpath_for_target(realevol)

#file(COPY _gf_retime.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#install(FILES _gf_retime.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${CMAKE_PROJECT_NAME})
