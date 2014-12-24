#include <triqs/utility/first_include.hpp>

#include <string>
#include <array>
#include <fstream>

#include <triqs/gfs/meshes/segment.hpp>
#include <triqs/h5.hpp>

#include "../c++/time_expr_r.hpp"
#include "../c++/solver.hpp"

using namespace realevol;

using triqs::utility::many_body_operator;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;
using triqs::gfs::segment_mesh;

namespace h5 = triqs::h5;

double hbar = 1.0;
double U = 1.0;
double mu = 0.5*U;
auto V = "0.3*(1 - exp(-4*t))"_te;
ode_solve_method method = method_runge_kutta;

std::string output_filename("trimer_rk.output.h5");

#include "./trimer.hpp"