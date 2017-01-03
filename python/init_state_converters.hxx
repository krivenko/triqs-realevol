// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/init_state.hpp --only_converters -I../c++ --only wrap_eq_solver_parameters -N realevol


// --- C++ Python converter for eq_solver_parameters_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<eq_solver_parameters_t> {
 static PyObject *c2py(eq_solver_parameters_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "verbosity"             , convert_to_python(x.verbosity));
  PyDict_SetItemString( d, "min_rel_weight"        , convert_to_python(x.min_rel_weight));
  PyDict_SetItemString( d, "arpack_min_matrix_size", convert_to_python(x.arpack_min_matrix_size));
  PyDict_SetItemString( d, "arpack_tolerance"      , convert_to_python(x.arpack_tolerance));
  PyDict_SetItemString( d, "arpack_ncv"            , convert_to_python(x.arpack_ncv));
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

 static eq_solver_parameters_t py2c(PyObject *dic) {
  eq_solver_parameters_t res;
  _get_optional(dic, "verbosity"             , res.verbosity                ,0);
  _get_optional(dic, "min_rel_weight"        , res.min_rel_weight           ,std::numeric_limits<double>::epsilon());
  _get_optional(dic, "arpack_min_matrix_size", res.arpack_min_matrix_size   ,101);
  _get_optional(dic, "arpack_tolerance"      , res.arpack_tolerance         ,0);
  _get_optional(dic, "arpack_ncv"            , res.arpack_ncv               ,std::map<long,int>({}));
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
  if (dic == nullptr or !PyDict_Check(dic)) {
   if (raise_exception) { PyErr_SetString(PyExc_TypeError, "The function must be called with named arguments");}
   return false;
  }
  std::stringstream fs, fs2; int err=0;

#ifndef TRIQS_ALLOW_UNUSED_PARAMETERS
  std::vector<std::string> ks, all_keys = {"verbosity","min_rel_weight","arpack_min_matrix_size","arpack_tolerance","arpack_ncv"};
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

  _check_optional <int                >(dic, fs, err, "verbosity"             , "int");
  _check_optional <double             >(dic, fs, err, "min_rel_weight"        , "double");
  _check_optional <int                >(dic, fs, err, "arpack_min_matrix_size", "int");
  _check_optional <double             >(dic, fs, err, "arpack_tolerance"      , "double");
  _check_optional <std::map<long, int>>(dic, fs, err, "arpack_ncv"            , "std::map<long, int>");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class eq_solver_parameters_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}
