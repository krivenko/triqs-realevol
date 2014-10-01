#include "time_expr_r.hpp"

#include <boost/lexical_cast.hpp>

namespace realevol {

// Global exprtk::parser object
exprtk::parser<double> parser;

time_expr_r::time_expr_r(std::string const& str) :
    str(str)
{
    exprtk::symbol_table<double> symt;
    symt.add_variable("t",arg);
    symt.add_constants();
    expr.register_symbol_table(symt);
    parser.compile(str,expr);
    if(exprtk::expression_helper<double>::is_constant(expr))
        this->str = boost::lexical_cast<std::string>(expr());
}

time_expr_r::time_expr_r(const char* expr) : time_expr_r(std::string(expr))
{}

time_expr_r::time_expr_r(double r) : time_expr_r(boost::lexical_cast<std::string>(r))
{}

time_expr_r::time_expr_r(time_expr_r const& te) :
    str(te.str)
{
    exprtk::symbol_table<double> symt;
    symt.add_variable("t",arg);
    symt.add_constants();
    expr.register_symbol_table(symt);
    parser.compile(str,expr);
}

auto time_expr_r::operator()(double t) const -> double
{
    arg = t;
    return expr.value();
}

std::ostream& operator<<(std::ostream& os, time_expr_r const& te)
{
    return (os << te.str);
}

auto time_expr_r::operator-() const -> time_expr_r
{
    return time_expr_r("-(" + str + ")");
}

auto time_expr_r::operator=(const time_expr_r& te) -> time_expr_r &
{
    this->str = te.str;
    parser.compile(str,expr);
    return *this;
}

auto time_expr_r::operator=(std::string const& str) -> time_expr_r &
{
    this->str = str;
    parser.compile(str,expr);
    return *this;
}

auto time_expr_r::operator=(const char* expr) -> time_expr_r &
{
    *this = std::string(expr);
    return *this;
}

auto time_expr_r::operator=(double r) -> time_expr_r &
{
    *this = boost::lexical_cast<std::string>(r);
    return *this;
}

auto time_expr_r::operator+=(const time_expr_r& te) -> time_expr_r &
{
    str = "(" + str + ")+(" + te.str + ")";
    parser.compile(str,expr);
    if(exprtk::expression_helper<double>::is_constant(expr))
        str = boost::lexical_cast<std::string>(expr());
    return *this;
}

auto time_expr_r::operator-=(const time_expr_r& te) -> time_expr_r &
{
    str = "(" + str + ")-(" + te.str + ")";
    parser.compile(str,expr);
    if(exprtk::expression_helper<double>::is_constant(expr))
        str = boost::lexical_cast<std::string>(expr());
    return *this;
}

auto time_expr_r::operator*=(const time_expr_r& te) -> time_expr_r &
{
    str = "(" + str + ")*(" + te.str + ")";
    parser.compile(str,expr);
    if(exprtk::expression_helper<double>::is_constant(expr))
        str = boost::lexical_cast<std::string>(expr());
    return *this;
}

auto time_expr_r::operator/=(const time_expr_r& te) -> time_expr_r &
{
    str = "(" + str + ")/(" + te.str + ")";
    parser.compile(str,expr);
    if(exprtk::expression_helper<double>::is_constant(expr))
        str = boost::lexical_cast<std::string>(expr());
    return *this;
}

bool time_expr_r::operator==(const time_expr_r& te) const
{
    return str == te.str;
}

}