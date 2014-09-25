#pragma once

#include <ostream>
#include <complex>
#include <tuple>
#include <boost/operators.hpp>

#include "time_expr.hpp"

namespace realevol {

// Complex time-dependent expressions
class c_time_expr : public boost::operators<c_time_expr>
{
    time_expr re;
    time_expr im;

public:

    c_time_expr(time_expr const& re = time_expr(), time_expr const& im = time_expr()) : re(re), im(im) {};
    c_time_expr(std::complex<double> z) : re(z.real()), im(z.imag()) {}
    c_time_expr(double x) : re(x), im(.0) {}
    c_time_expr(std::string const& str) : re(str), im(.0) {}
    c_time_expr(const char* str) : re(str), im(.0) {}
    c_time_expr(c_time_expr const&) = default;
    c_time_expr(c_time_expr &&) = default;

    c_time_expr operator-() const { return {-re,-im}; };

    c_time_expr & operator=(c_time_expr const&) = default;
    c_time_expr & operator+=(c_time_expr const& cte) {
        re += cte.re; im += cte.im;
        return *this;
    }
    c_time_expr & operator-=(c_time_expr const& cte) {
        re -= cte.re; im -= cte.im;
        return *this;
    }
    c_time_expr & operator*=(c_time_expr const& cte) {
        std::tie(re,im) = std::make_tuple(re*cte.re - im*cte.im, re*cte.im + im*cte.re);
        return *this;
    }
    c_time_expr & operator/=(c_time_expr const& cte) {
        time_expr abs2 = cte.re*cte.re + cte.im*cte.im;
        std::tie(re,im) = std::make_tuple((re*cte.re + im*cte.im)/abs2,(im*cte.re - re*cte.im)/abs2);
        return *this;
    }
    bool operator==(c_time_expr const& cte) const { return (re==cte.re) && (im==cte.im); }

    std::complex<double> operator()(double t) const { return {re(t),im(t)}; }

    friend bool is_constant(c_time_expr const& cte) { return is_constant(cte.re) && is_constant(cte.im); }
    bool is_zero() const { return re.is_zero() && im.is_zero(); }

    time_expr real() const { return re; }
    time_expr imag() const { return im; }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, c_time_expr const& cte) {
        return (os << "(" << cte.re << "," << cte.im << ")");
    }

private:

    // Methods for boost::serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & TRIQS_MAKE_NVP("r",re) & TRIQS_MAKE_NVP("im",im);
    }
};

inline c_time_expr operator ""_cte(long double r){ return c_time_expr(r,.0); }
inline c_time_expr operator ""_cte(const char* expr, std::size_t) { return c_time_expr(expr); };

}

namespace triqs { namespace utility {

inline bool is_zero(realevol::c_time_expr const& cte, double = 0 /* neglected */) { return cte.is_zero(); }
inline realevol::c_time_expr _conj(realevol::c_time_expr const& cte) { return {cte.real(),-cte.imag()}; }

}}
