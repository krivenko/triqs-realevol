find_package(TriqsTest)

set(PYTHON_TESTS texpr)
foreach(t ${PYTHON_TESTS})
  triqs_add_python_test(${t})
endforeach(t)
