#pragma once

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
#include "callable_complex.hpp"

using boost::lexical_cast;
using triqs::utility::numeric_ops;

namespace realevol {

// Time-dependent expressions;

class time_expr : public boost::operators<time_expr>
{

public:

    using result_type = double;

    time_expr();
    time_expr(std::string const& expr);
    time_expr(const char* expr);
    time_expr(result_type r);
    time_expr(time_expr const& te);
    ~time_expr();
    result_type operator()(double t) const;

    friend bool is_constant(time_expr const& te)
    {
        return te.parser==nullptr;
    }
    bool is_zero() const {
        return parser==nullptr && numeric_ops<double>::is_zero(arg);
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, time_expr const& te);

    time_expr operator-() const;

    time_expr operator=(time_expr const& te);
    time_expr operator=(std::string const& expr);
    time_expr operator=(const char* expr);
    time_expr operator=(result_type r);

    time_expr operator+=(time_expr const& te);
    time_expr operator-=(time_expr const& te);
    time_expr operator*=(time_expr const& te);
    time_expr operator/=(time_expr const& te);
    bool operator==(time_expr const& te) const;

private:

    mutable mu::Parser* parser;
    mutable result_type arg;

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

using complex_time_expr = callable_complex<time_expr>;

// Replace the expression with a constant if it takes equal values at all mesh points
template<class Mesh, class Expr>
bool try_reduce_to_constant(Expr& te, Mesh const& m)
{
    auto it = m.begin();
    auto value = te(*it);

    for(it++; it != m.end(); it++)
        if(!numeric_ops<decltype(value)>::is_zero(te(*it) - value)) return false;

    te = Expr(value);
    return true;
}

}

namespace triqs { namespace utility {
    template<> struct numeric_ops<realevol::time_expr> {
        static bool is_zero(realevol::time_expr const& te) { return te.is_zero(); }
        static realevol::time_expr conj(realevol::time_expr const& te) { return te; }
    };
}}
