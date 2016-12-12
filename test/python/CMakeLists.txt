find_package(TriqsTest)

set(PYTHON_TESTS texpr operators gf_retime)
foreach(t ${PYTHON_TESTS})
  triqs_add_python_test(${t})
endforeach(t)
