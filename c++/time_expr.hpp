#pragma once

#pragma once

#include <complex>
#include <string>
#include <ostream>
#include <limits>
#include <type_traits>
#include <muParser.h>
#include <boost/operators.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>

#include "numerics.hpp"
#include "mesh_base.hpp"
#include "mesh_iterator.hpp"

using boost::lexical_cast;

namespace realevol {
    
// Time-dependent expressions
template<bool ComplexValued = false> class time_expr;
    
template<>
class time_expr<false> : public boost::operators<time_expr<false>>
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
    bool is_constant() const;

    friend bool is_zero(time_expr const& te)
    {
        return te.is_constant() && is_zero(te.arg);
    }
    
    // Stream output
    friend std::ostream& operator<<(std::ostream& os, time_expr const& te);
    
    time_expr operator-() const;
    
    time_expr operator=(time_expr const& te);
    time_expr operator+=(time_expr const& te);
    time_expr operator-=(time_expr const& te);
    time_expr operator*=(time_expr const& te);
    time_expr operator/=(time_expr const& te);
    bool operator==(time_expr const& te) const;
    
    friend time_expr _conj(time_expr const& te)
    {
        return te;
    }
    
private:
    
    mutable mu::Parser* parser;
    mutable result_type arg;
    
    // Methods for boost::serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar & (is_constant() ? lexical_cast<std::string>(arg) : parser->GetExpr());
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

// Time-dependent expressions (complex)
template<>
class time_expr<true> : public boost::operators<time_expr<true>>
{

public:
    
    using result_type = std::complex<double>;
    
    time_expr();
    time_expr(std::string const& expr_re, std::string const& expr_im = "0");
    time_expr(const char* expr_re, const char* expr_im = "0");
    time_expr(result_type c);
    time_expr(result_type::value_type r);
    time_expr(time_expr const& te);
    result_type operator()(double t) const;
    bool is_constant() const;

    friend bool is_zero(time_expr const& te)
    {
        return te.is_constant() && is_zero(te.precomputed_value);
    }
    
    // Stream output
    friend std::ostream& operator<<(std::ostream& os, time_expr<true> const& te);

    time_expr operator-() const;
    
    time_expr operator=(time_expr const& te);
    time_expr operator+=(time_expr const& te);
    time_expr operator-=(time_expr const& te);
    time_expr operator*=(time_expr const& te);
    time_expr operator/=(time_expr const& te);
    bool operator==(time_expr const& te) const;
    
    friend time_expr _conj(time_expr const& te)
    {
        if(te.is_constant())
            return _conj(te.precomputed_value);
        else
            return time_expr(te.parser_re.GetExpr(),"-(" + te.parser_im.GetExpr() + ")");
    }
    
private:
    
    mu::Parser parser_re, parser_im;
    
    mutable double t;
    
    mutable bool _is_constant;
    mutable result_type precomputed_value;
    
    // Methods for boost::serialization
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar & parser_re.GetExpr();
        ar & parser_im.GetExpr();
        ar & _is_constant;
        if(_is_constant) ar & precomputed_value;
    }
    
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        std::string expr;
        ar & expr;
        parser_re.SetExpr(expr);
        ar & expr;
        parser_im.SetExpr(expr);
        ar & _is_constant;
        if(_is_constant) ar & precomputed_value;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

// Replace the expression with a constant if it takes equal values at all mesh points
template<class Mesh, bool ComplexValued>
bool try_reduce_to_constant(time_expr<ComplexValued>& te, Mesh const& m)
{
    auto it = m.begin();
    auto value = te(*it);
    
    for(it++; it != m.end(); it++)
        if(!is_zero(te(*it) - value)) return false;
        
    te = time_expr<ComplexValued>(value);
    return true;
}

}
