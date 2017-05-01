// DO NOT EDIT
// Generated automatically using libclang using the command :
// c++2py.py ../c++/solver.hpp -I../c++ --only_converters -p -m realevol -o realevol --appname realevol --moduledoc "The Real-time evolution solver"


// --- C++ Python converter for compute_2t_obs_parameters_t
#include <triqs/python_tools/converters/vector.hpp>
#include <triqs/python_tools/converters/string.hpp>
#include <algorithm>

namespace triqs { namespace py_tools {

template <> struct py_converter<compute_2t_obs_parameters_t> {
 static PyObject *c2py(compute_2t_obs_parameters_t const & x) {
  PyObject * d = PyDict_New();
  PyDict_SetItemString( d, "h"                      , convert_to_python(x.h));
  PyDict_SetItemString( d, "verbosity"              , convert_to_python(x.verbosity));
  PyDict_SetItemString( d, "hbar"                   , convert_to_python(x.hbar));
  PyDict_SetItemString( d, "hamiltonian_interpol"   , convert_to_python(x.hamiltonian_interpol));
  PyDict_SetItemString( d, "lanczos_min_matrix_size", convert_to_python(x.lanczos_min_matrix_size));
  PyDict_SetItemString( d, "lanczos_gs_energy_tol"  , convert_to_python(x.lanczos_gs_energy_tol));
  PyDict_SetItemString( d, "lanczos_max_krylov_dim" , convert_to_python(x.lanczos_max_krylov_dim));
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

 static compute_2t_obs_parameters_t py2c(PyObject *dic) {
  compute_2t_obs_parameters_t res;
  res.h = convert_from_python<operator_t>(PyDict_GetItemString(dic, "h"));
  _get_optional(dic, "verbosity"              , res.verbosity                 ,((triqs::mpi::communicator().rank()==0)?3:0));
  _get_optional(dic, "hbar"                   , res.hbar                      ,1.0);
  _get_optional(dic, "hamiltonian_interpol"   , res.hamiltonian_interpol      ,Rectangle);
  _get_optional(dic, "lanczos_min_matrix_size", res.lanczos_min_matrix_size   ,11);
  _get_optional(dic, "lanczos_gs_energy_tol"  , res.lanczos_gs_energy_tol     ,std::map<long,double>({}));
  _get_optional(dic, "lanczos_max_krylov_dim" , res.lanczos_max_krylov_dim    ,std::map<long,int>({}));
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
  std::vector<std::string> ks, all_keys = {"h","verbosity","hbar","hamiltonian_interpol","lanczos_min_matrix_size","lanczos_gs_energy_tol","lanczos_max_krylov_dim"};
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

  _check_mandatory<operator_t               >(dic, fs, err, "h"                      , "operator_t");
  _check_optional <int                      >(dic, fs, err, "verbosity"              , "int");
  _check_optional <double                   >(dic, fs, err, "hbar"                   , "double");
  _check_optional <realevol::h_interpolation>(dic, fs, err, "hamiltonian_interpol"   , "realevol::h_interpolation");
  _check_optional <int                      >(dic, fs, err, "lanczos_min_matrix_size", "int");
  _check_optional <std::map<long, double>   >(dic, fs, err, "lanczos_gs_energy_tol"  , "std::map<long, double>");
  _check_optional <std::map<long, int>      >(dic, fs, err, "lanczos_max_krylov_dim" , "std::map<long, int>");
  if (err) goto _error;
  return true;

 _error:
   fs2 << "\n---- There " << (err > 1 ? "are " : "is ") << err<< " error"<<(err >1 ?"s" : "")<< " in Python -> C++ transcription for the class compute_2t_obs_parameters_t\n" <<fs.str();
   if (raise_exception) PyErr_SetString(PyExc_TypeError, fs2.str().c_str());
  return false;
 }
};

}}