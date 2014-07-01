#include "time_expr.hpp"

#include <limits>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

namespace realevol {

time_expr<true>::time_expr(std::string const& expr_re, std::string const& expr_im)
{
    parser_re.DefineVar("t", &t);
    parser_im.DefineVar("t", &t);
    
    parser_re.SetExpr(expr_re);
    parser_im.SetExpr(expr_im);
    
    _is_constant = (!parser_re.GetUsedVar().size()) && (!parser_im.GetUsedVar().size());
    if(_is_constant){
        precomputed_value = std::complex<double>(parser_re.Eval(), parser_im.Eval());
    }
}

time_expr<true>::time_expr(const char* expr_re, const char* expr_im)
{
    parser_re.DefineVar("t", &t);
    parser_im.DefineVar("t", &t);
    
    parser_re.SetExpr(expr_re);
    parser_im.SetExpr(expr_im);
    
    _is_constant = (!parser_re.GetUsedVar().size()) && (!parser_im.GetUsedVar().size());
    if(_is_constant){
        precomputed_value = std::complex<double>(parser_re.Eval(), parser_im.Eval());
    }
}

time_expr<true>::time_expr(result_type c) :
    _is_constant(true),
    precomputed_value(c)
{
    parser_re.DefineVar("t", &t);
    parser_im.DefineVar("t", &t);
    
    parser_re.SetExpr(boost::lexical_cast<std::string>(real(c)));
    parser_im.SetExpr(boost::lexical_cast<std::string>(imag(c)));
}

time_expr<true>::time_expr(result_type::value_type r) :
    _is_constant(true),
    precomputed_value(r,0.0)
{
    parser_re.DefineVar("t", &t);
    parser_im.DefineVar("t", &t);
    
    parser_re.SetExpr(boost::lexical_cast<std::string>(r));
    parser_im.SetExpr("0");
}

time_expr<true>::time_expr(void) :
    _is_constant(true),
    precomputed_value(0.0)
{
    parser_re.DefineVar("t", &t);
    parser_im.DefineVar("t", &t);
    
    parser_re.SetExpr("0");
    parser_im.SetExpr("0");
}

time_expr<true>::time_expr(time_expr const& te) :
    _is_constant(te._is_constant),
    precomputed_value(te.precomputed_value)
{
    parser_re.DefineVar("t",&t);
    parser_im.DefineVar("t",&t);
    
    parser_re.SetExpr(te.parser_re.GetExpr());
    parser_im.SetExpr(te.parser_im.GetExpr()); 
}

auto time_expr<true>::operator()(double t) const -> result_type
{
    if(_is_constant)
        return precomputed_value;
    else {
        this->t = t;
        return std::complex<double>(parser_re.Eval(), parser_im.Eval());
    }
}

bool time_expr<true>::is_constant() const
{
    return _is_constant;
}

bool time_expr<true>::is_zero() const
{
    return _is_constant && std::abs(precomputed_value) < comparison_tolerance;
}

std::ostream& operator<<(std::ostream& os, time_expr<true> const& te)
{
    if(te._is_constant)
        os << te.precomputed_value;
    else{
        os << "(" << boost::trim_copy(te.parser_re.GetExpr());
        os << "," << boost::trim_copy(te.parser_im.GetExpr()) << ")";
    }
    return os;
}

auto time_expr<true>::operator-(void) -> time_expr
{
    if(_is_constant)
        return time_expr(boost::lexical_cast<std::string>(-precomputed_value));
    else
        return time_expr("-(" + parser_re.GetExpr() + ")", "-(" + parser_im.GetExpr() + ")");
}

auto time_expr<true>::operator=(const time_expr& te) -> time_expr
{
    parser_re.SetExpr(te.parser_re.GetExpr());
    parser_im.SetExpr(te.parser_im.GetExpr());
    _is_constant = te._is_constant;
    if(_is_constant)
        precomputed_value = te.precomputed_value;
    return *this;
}

auto time_expr<true>::operator+=(const time_expr& te) -> time_expr
{
    _is_constant = _is_constant && te._is_constant;
    if(_is_constant){
        precomputed_value += te.precomputed_value;
        parser_re.SetExpr(boost::lexical_cast<std::string>(real(precomputed_value)));
        parser_im.SetExpr(boost::lexical_cast<std::string>(imag(precomputed_value)));
    } else {
        parser_re.SetExpr("(" + parser_re.GetExpr() + ")+(" + te.parser_re.GetExpr() + ")");
        parser_im.SetExpr("(" + parser_im.GetExpr() + ")+(" + te.parser_im.GetExpr() + ")");
    }
    return *this;
}

auto time_expr<true>::operator-=(const time_expr& te) -> time_expr
{
    _is_constant = _is_constant && te._is_constant;
    if(_is_constant){
        precomputed_value -= te.precomputed_value;
        parser_re.SetExpr(boost::lexical_cast<std::string>(real(precomputed_value)));
        parser_im.SetExpr(boost::lexical_cast<std::string>(imag(precomputed_value)));
    } else {
        parser_re.SetExpr("(" + parser_re.GetExpr() + ")-(" + te.parser_re.GetExpr() + ")");
        parser_im.SetExpr("(" + parser_im.GetExpr() + ")-(" + te.parser_im.GetExpr() + ")");
    }
    return *this;
}

auto time_expr<true>::operator*=(const time_expr& te) -> time_expr
{
    _is_constant = _is_constant && te._is_constant;
    if(_is_constant){
        precomputed_value *= te.precomputed_value;
        parser_re.SetExpr(boost::lexical_cast<std::string>(real(precomputed_value)));
        parser_im.SetExpr(boost::lexical_cast<std::string>(imag(precomputed_value)));
    } else {
        std::string this_re = parser_re.GetExpr(), te_re = te.parser_re.GetExpr();
        std::string this_im = parser_im.GetExpr(), te_im = te.parser_im.GetExpr();
        parser_re.SetExpr("(" + this_re + ")*(" + te_re + ") - (" + this_im + ")*(" + te_im + ")");
        parser_im.SetExpr("(" + this_re + ")*(" + te_im + ") + (" + this_im + ")*(" + te_re + ")");
    }
    return *this;
}

auto time_expr<true>::operator/=(const time_expr& te) -> time_expr
{
    _is_constant = _is_constant && te._is_constant;
    if(_is_constant){
        precomputed_value /= te.precomputed_value;
        parser_re.SetExpr(boost::lexical_cast<std::string>(real(precomputed_value)));
        parser_im.SetExpr(boost::lexical_cast<std::string>(imag(precomputed_value)));
    } else {
        std::string this_re = parser_re.GetExpr(), te_re = te.parser_re.GetExpr();
        std::string this_im = parser_im.GetExpr(), te_im = te.parser_im.GetExpr();
        std::string denom = "(" + te_re + ")*(" + te_re + ") + (" + te_im + ")*(" + te_im + ")";
        parser_re.SetExpr("((" + this_re + ")*(" + te_re + ") + (" + this_im + ")*(" + te_im + ")) / (" + denom + ")");
        parser_im.SetExpr("((" + this_im + ")*(" + te_re + ") - (" + this_re + ")*(" + te_im + ")) / (" + denom + ")");
    }
    return *this;
}

bool time_expr<true>::operator==(const time_expr& te) const
{
    return ((_is_constant == te._is_constant) &&
                (std::abs(precomputed_value - te.precomputed_value) < time_expr::comparison_tolerance))
        || ((parser_re.GetExpr() == te.parser_re.GetExpr()) && (parser_im.GetExpr() == te.parser_im.GetExpr()));
}
    
}