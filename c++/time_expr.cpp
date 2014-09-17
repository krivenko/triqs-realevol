#include "time_expr.hpp"

#include <boost/lexical_cast.hpp>

namespace realevol {

// Global exprtk::parser object
exprtk::parser<double> parser;

time_expr::time_expr(std::string const& str) :
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

time_expr::time_expr(const char* expr) : time_expr(std::string(expr))
{}

time_expr::time_expr(double r) : time_expr(boost::lexical_cast<std::string>(r))
{}

time_expr::time_expr(time_expr const& te) :
    str(te.str)
{
    exprtk::symbol_table<double> symt;
    symt.add_variable("t",arg);
    symt.add_constants();
    expr.register_symbol_table(symt);
    parser.compile(str,expr);
}

auto time_expr::operator()(double t) const -> double
{
    arg = t;
    return expr.value();
}

std::ostream& operator<<(std::ostream& os, time_expr const& te)
{
    return (os << te.str);
}

auto time_expr::operator-() const -> time_expr
{
    return time_expr("-(" + str + ")");
}

auto time_expr::operator=(const time_expr& te) -> time_expr
{
    this->str = te.str;
    parser.compile(str,expr);
    return *this;
}

auto time_expr::operator=(std::string const& str) -> time_expr
{
    this->str = str;
    parser.compile(str,expr);
    return *this;
}

auto time_expr::operator=(const char* expr) -> time_expr
{
    *this = std::string(expr);
    return *this;
}

auto time_expr::operator=(double r) -> time_expr
{
    *this = boost::lexical_cast<std::string>(r);
    return *this;
}

auto time_expr::operator+=(const time_expr& te) -> time_expr
{
    str = "(" + str + ")+(" + te.str + ")";
    parser.compile(str,expr);
    if(exprtk::expression_helper<double>::is_constant(expr))
        str = boost::lexical_cast<std::string>(expr());
    return *this;
}

auto time_expr::operator-=(const time_expr& te) -> time_expr
{
    str = "(" + str + ")-(" + te.str + ")";
    parser.compile(str,expr);
    if(exprtk::expression_helper<double>::is_constant(expr))
        str = boost::lexical_cast<std::string>(expr());
    return *this;
}

auto time_expr::operator*=(const time_expr& te) -> time_expr
{
    str = "(" + str + ")*(" + te.str + ")";
    parser.compile(str,expr);
    if(exprtk::expression_helper<double>::is_constant(expr))
        str = boost::lexical_cast<std::string>(expr());
    return *this;
}

auto time_expr::operator/=(const time_expr& te) -> time_expr
{
    str = "(" + str + ")/(" + te.str + ")";
    parser.compile(str,expr);
    if(exprtk::expression_helper<double>::is_constant(expr))
        str = boost::lexical_cast<std::string>(expr());
    return *this;
}

bool time_expr::operator==(const time_expr& te) const
{
    return str == te.str;
}

}