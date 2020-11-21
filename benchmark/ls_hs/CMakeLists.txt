add_executable(ls_hs ls_hs.cpp)
target_link_libraries(ls_hs realevol_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_SOURCE_DIR}/c++ ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(ls_hs)
