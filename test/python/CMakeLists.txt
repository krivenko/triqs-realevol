find_package(TriqsTest)

set(PYTHON_TESTS texpr ctexpr meshes operator)
foreach(t ${PYTHON_TESTS})
    triqs_add_test_txt(py_${t})
endforeach(t)