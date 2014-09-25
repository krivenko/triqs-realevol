#pragma once

#include <string>
#include <ostream>
#include <boost/operators.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>

#include <triqs/utility/draft/numeric_ops.hpp>
#include <triqs/utility/mini_vector.hpp>

// Disable some features of ExprTk to speed up compilation
#define exprtk_disable_comments
#define exprtk_disable_break_continue
#define exprtk_disable_string_capabilities
#include "exprtk/exprtk.hpp"

namespace realevol {

// Time-dependent expressions;
class time_expr : public boost::operators<time_expr>
{
    std::string str;
    exprtk::expression<double> expr;
    mutable double arg;

public:

    time_expr(std::string const& str);
    time_expr(const char* str);
    time_expr(double r = 0);
    time_expr(time_expr const& te);

    time_expr operator-() const;

    time_expr & operator=(time_expr const& te);
    time_expr & operator=(std::string const& str);
    time_expr & operator=(const char* str);
    time_expr & operator=(double r);

    time_expr & operator+=(time_expr const& te);
    time_expr & operator-=(time_expr const& te);
    time_expr & operator*=(time_expr const& te);
    time_expr & operator/=(time_expr const& te);
    bool operator==(time_expr const& te) const;

    double operator()(double t) const;

    friend bool is_constant(time_expr const& te)
    {
        return exprtk::expression_helper<double>::is_constant(te.expr);
    }
    bool is_zero() const {
        return is_constant(*this) && triqs::utility::is_zero(expr.value());
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, time_expr const& te);

private:

    // Methods for boost::serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar & TRIQS_MAKE_NVP("str",str);
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        std::string str;
        ar & str;
        *this = time_expr(str);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

inline time_expr operator ""_te (long double r){ return time_expr(r); }
inline time_expr operator ""_te(const char* expr, std::size_t) { return time_expr(expr); };

// Replace the expression with a constant if it takes equal values at all mesh points
template<class Mesh, class Expr>
bool try_reduce_to_constant(Expr& te, Mesh const& m)
{
    auto it = m.begin();
    auto value = te(*it);

    for(it++; it != m.end(); it++)
        if(!triqs::utility::is_zero(te(*it) - value)) return false;

    te = Expr(value);
    return true;
}

}

namespace triqs { namespace utility {

inline bool is_zero(realevol::time_expr const& te, double = 0 /* neglected */) { return te.is_zero(); }
inline realevol::time_expr _conj(realevol::time_expr const& te) { return te; }

}}
