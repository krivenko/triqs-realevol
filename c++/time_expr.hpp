#pragma once

#include <string>
#include <ostream>
#include <limits>
#include <type_traits>
#include <muParser.h>
#include <boost/operators.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>

#include <triqs/utility/draft/numeric_ops.hpp>

#include "mesh_base.hpp"
#include "mesh_iterator.hpp"

using boost::lexical_cast;

namespace realevol {

// Time-dependent expressions;

class time_expr : public boost::operators<time_expr>
{

public:

    time_expr(std::string const& expr);
    time_expr(const char* expr);
    time_expr(double r = 0);
    time_expr(time_expr const& te);
    ~time_expr();

    time_expr operator-() const;

    time_expr operator=(time_expr const& te);
    time_expr operator=(std::string const& expr);
    time_expr operator=(const char* expr);
    time_expr operator=(double r);

    time_expr operator+=(time_expr const& te);
    time_expr operator-=(time_expr const& te);
    time_expr operator*=(time_expr const& te);
    time_expr operator/=(time_expr const& te);
    bool operator==(time_expr const& te) const;

    double operator()(double t) const;

    friend bool is_constant(time_expr const& te)
    {
        return te.parser==nullptr;
    }
    bool is_zero() const {
        return parser==nullptr && triqs::utility::is_zero(arg);
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, time_expr const& te);

private:

    mutable mu::Parser* parser;
    mutable double arg;

    // Methods for boost::serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar & (is_constant(*this) ? lexical_cast<std::string>(arg) : parser->GetExpr());
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        std::string expr;
        ar & expr;
        *this = time_expr(expr);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

time_expr operator""_te (long double r){ return time_expr(r); }
time_expr operator ""_te(const char* expr, std::size_t) { return time_expr(expr); };

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

bool is_zero(realevol::time_expr const& te) { return te.is_zero(); }
realevol::time_expr _conj(realevol::time_expr const& te) { return te; }

}}
