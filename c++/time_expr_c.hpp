#pragma once

#include <ostream>
#include <complex>
#include <tuple>
#include <boost/operators.hpp>

#include "time_expr_r.hpp"

namespace realevol {

// Complex time-dependent expressions
class time_expr_c : public boost::operators<time_expr_c>
{
    time_expr_r re;
    time_expr_r im;

public:

    time_expr_c(time_expr_r const& re = time_expr_r(), time_expr_r const& im = time_expr_r()) : re(re), im(im) {};
    time_expr_c(std::complex<double> z) : re(z.real()), im(z.imag()) {}
    time_expr_c(double x) : re(x), im(.0) {}
    time_expr_c(std::string const& str) : re(str), im(.0) {}
    time_expr_c(const char* str) : re(str), im(.0) {}
    time_expr_c(time_expr_c const&) = default;
    time_expr_c(time_expr_c &&) = default;

    time_expr_c operator-() const { return {-re,-im}; };

    time_expr_c & operator=(time_expr_c const&) = default;
    time_expr_c & operator+=(time_expr_c const& cte) {
        re += cte.re; im += cte.im;
        return *this;
    }
    time_expr_c & operator-=(time_expr_c const& cte) {
        re -= cte.re; im -= cte.im;
        return *this;
    }
    time_expr_c & operator*=(time_expr_c const& cte) {
        std::tie(re,im) = std::make_tuple(re*cte.re - im*cte.im, re*cte.im + im*cte.re);
        return *this;
    }
    time_expr_c & operator/=(time_expr_c const& cte) {
        time_expr_r abs2 = cte.re*cte.re + cte.im*cte.im;
        std::tie(re,im) = std::make_tuple((re*cte.re + im*cte.im)/abs2,(im*cte.re - re*cte.im)/abs2);
        return *this;
    }
    bool operator==(time_expr_c const& cte) const { return (re==cte.re) && (im==cte.im); }

    std::complex<double> operator()(double t) const { return {re(t),im(t)}; }

    friend bool is_constant(time_expr_c const& cte) { return is_constant(cte.re) && is_constant(cte.im); }
    bool is_zero() const { return re.is_zero() && im.is_zero(); }

    time_expr_r real() const { return re; }
    time_expr_r imag() const { return im; }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, time_expr_c const& cte) {
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

inline time_expr_c operator ""_cte(long double r){ return time_expr_c(r,.0); }
inline time_expr_c operator ""_cte(const char* expr, std::size_t) { return time_expr_c(expr); };

}

namespace triqs { namespace utility {

inline bool is_zero(realevol::time_expr_c const& cte, double = 0 /* neglected */) { return cte.is_zero(); }
inline realevol::time_expr_c _conj(realevol::time_expr_c const& cte) { return {cte.real(),-cte.imag()}; }

}}
