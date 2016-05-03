// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/solver.hpp -p -m realevol -o realevol --appname realevol --moduledoc "The Real-time evolution solver"


// --- C++ Python converter for solve_parameters_t

namespace triqs { namespace py_tools {

template <> struct py_converter<solve_parameters_t> {
 static PyObject *c2py(solve_parameters_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "h"        , convert_to_python(x.h));
  PyDict_SetItemString( d, "verbosity", convert_to_python(x.verbosity));
  PyDict_SetItemString( d, "hbar"     , convert_to_python(x.hbar));
  PyDict_SetItemString( d, "method"   , convert_to_python(x.method));
  return d;
 }

 template <typename T, typename U> static void _get_optional(PyObject *dic, const char *name, T &r, U const &init_default) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = init_default;
 }

 template <typename T> static void _get_optional(PyObject *dic, const char *name, T &r) {
  if (PyDict_Contains(dic, pyref::string(name)))
   r = convert_from_python<T>(PyDict_GetItemString(dic, name));
  else
   r = T{};
 }

 static solve_parameters_t py2c(PyObject *dic) {
  solve_parameters_t res;
  res.h = convert_from_python<operator_t>(PyDict_GetItemString(dic, "h"));
  _get_optional(dic, "verbosity", res.verbosity   ,((triqs::mpi::communicator().rank()==0)?3:0));
  _get_optional(dic, "hbar"     , res.hbar        ,1.0);
  _get_optional(dic, "method"   , res.method      ,Lanczos);
  return res;
 }

 template <typename T>
 static void _check(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!convertible_from_python<T>(PyDict_GetItemString(dic, name), false))
   fs << "\n" << ++err << " The parameter " << name << " does not have the right type : expecting " << tname
      << " in C++, but got '" << PyDict_GetItemString(dic, name)->ob_type->tp_name << "' in Python.";
 }

 template <typename T>
 static void _check_mandatory(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (!PyDict_Contains(dic, pyref::string(name)))
   fs << "\n" << ++err << " Mandatory parameter " << name << " is missing.";
  else _check<T>(dic,fs,err,name,tname);
 }

 template <typename T>
 static void _check_optional(PyObject *dic, std::stringstream &fs, int &err, const char *name, const char *tname) {
  if (PyDict_Contains(dic, pyref::string(name))) _check<T>(dic, fs, err, name, tname);
 }

 static bool is_convertible(PyObject *dic, bool raise_exception) {
  if (!PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "Not a python dict");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"h","verbosity","hbar","method"};
  pyref keys = PyDict_Keys(dic);
  if (!convertible_from_python<std::vector<std::string>>(keys, true)) {
   fs << "\nThe dict keys are not strings";
   goto _error;
  }
  ks = convert_from_python<std::vector<std::string>>(keys);
  for (auto & k : ks)
   if (std::find(all_keys.begin(), all_keys.end(), k) == all_keys.end())
    fs << "\n"<< ++err << " The parameter '" << k << "' is not recognized.";
#endif

  _check_mandatory<operator_t                >(dic, fs, err, "h"        , "operator_t");
  _check_optional <int                       >(dic, fs, err, "verbosity", "int");
  _check_optional <double                    >(dic, fs, err, "hbar"     , "double");
  _check_optional <realevol::ode_solve_method>(dic, fs, err, "method"   , "realevol::ode_solve_method");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class solve_parameters_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}
