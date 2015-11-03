find_package(TriqsTest)

set(PYTHON_TESTS texpr)
foreach(t ${PYTHON_TESTS})
  triqs_add_test_txt(${t})
endforeach(t)
