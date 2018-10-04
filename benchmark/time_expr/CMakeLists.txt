add_executable(time_expr_benchmark time_expr_benchmark.cpp)
target_link_libraries(time_expr_benchmark time_expr_c ${TRIQS_LIBRARY_ALL} benchmark)
include_directories(${CMAKE_SOURCE_DIR}/c++ ${TRIQS_INCLUDE_ALL})
