from wrap_generator import *

module = module_(full_name = "tinterp", app_name = "realevol", doc = "Linear interpolator on time mesh")

module.use_module('gf', 'triqs')

module.add_include("time_interp.hpp")
module.add_include("<triqs/python_tools/converters/vector.hpp>")
module.add_include("<triqs/python_tools/converters/arrays.hpp>")
module.add_include("<triqs/python_tools/converters/gf.hpp>")

module.add_using("namespace realevol")
module.add_using("triqs::utility::is_zero")
module.add_using("triqs::utility::conj")
module.add_using("namespace triqs::gfs")

c = class_(
    py_type = "TInterp",
    c_type = "time_interp",
    c_type_absolute = "realevol::time_interp",
    is_printable = True,
    arithmetic = ("algebra","with_unit","with_unary_minus","double","std::complex<double>"),
    comparisons = "==",
    doc = "Linear interpolator on a time mesh"
    )

c.add_constructor(signature="()", doc="Create zero interpolator")
c.add_constructor(signature="(gf_mesh<retime> m)", doc="Create zero interpolator")
c.add_constructor(signature="(double r)", doc="Create real constant interpolator")
c.add_constructor(signature="(gf_mesh<retime> m, double r)", doc="Create real constant interpolator")
c.add_constructor(signature="(std::complex<double> z)", doc="Create complex constant interpolator")
c.add_constructor(signature="(gf_mesh<retime> m, std::complex<double> z)", doc="Create complex constant interpolator")
c.add_constructor(signature="(double r, double i)", doc="Create complex constant interpolator out of two real numbers")
c.add_constructor(signature="(gf_mesh<retime> m, double r, double i)", doc="Create complex constant interpolator out of two real numbers")
c.add_constructor(signature="(gf_mesh<retime> m, triqs::arrays::array<double, 1> r)", doc="Create real interpolator")
c.add_constructor(signature="(gf_mesh<retime> m, triqs::arrays::array<double, 1> r, double i)",
                  doc="Create complex interpolator with a constant imaginary part")
c.add_constructor(signature="(gf_mesh<retime> m, double r, triqs::arrays::array<double, 1> i)",
                  doc="Create complex interpolator with a constant real part")
c.add_constructor(signature="(gf_mesh<retime> m, triqs::arrays::array<double, 1> r, triqs::arrays::array<double, 1> i)",
                  doc="Create complex interpolator")
c.add_constructor(signature="(gf_mesh<retime> m, triqs::arrays::array<std::complex<double>, 1> z)",
                  doc="Create complex interpolator")

c.add_call(signature="std::complex<double>(double t)", doc="Interpolation result at time point t")
module.add_class(c)

module.add_function(signature="bool is_zero(time_interp ti)", doc="Boolean : is zero interpolator?")
module.add_function(signature="bool is_constant(time_interp ti)", doc="Boolean : is constant interpolator?")
module.add_function(signature="time_interp conj(time_interp ti)", doc="Complex conjugate of the interpolator")

module.generate_code()
