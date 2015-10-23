/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2015 I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#include "time_expr.hpp"

namespace realevol {

// Global exprtk::parser object
exprtk::parser<double> parser;

// Constructors
time_expr::time_expr() : time_expr(.0) {}

time_expr::time_expr(double r) : time_expr(std::to_string(r)) {}

time_expr::time_expr(double r, double i) : time_expr(std::to_string(r),std::to_string(i)) {}

time_expr::time_expr(std::complex<double> const& z) : time_expr(z.real(),z.imag()) {}

time_expr::time_expr(const char* str) : time_expr(std::string(str)) {}

time_expr::time_expr(const char* re_str, const char* im_str) :
  time_expr(std::string(re_str),std::string(im_str)) {}

time_expr::time_expr(std::string const& re_str, double i) : time_expr(re_str,std::to_string(i)) {}
time_expr::time_expr(double r, std::string const& im_str) : time_expr(std::to_string(r),im_str) {}
time_expr::time_expr(const char* re_str, double i) : time_expr(std::string(re_str),std::to_string(i)) {}
time_expr::time_expr(double r, const char* im_str) : time_expr(std::to_string(r),std::string(im_str)) {}

time_expr::time_expr(std::string const& str) :
  _is_real(true), re_str(str) {
  exprtk::symbol_table<double> symt;
  symt.add_variable("t",arg);
  symt.add_constants();
  re.register_symbol_table(symt);
  parser.compile(str,re);
  if(exprtk::expression_helper<double>::is_constant(re))
    this->re_str = std::to_string(re());
}

time_expr::time_expr(std::string const& re_str, std::string const& im_str) :
  _is_real(false), re_str(re_str), im_str(im_str) {
  exprtk::symbol_table<double> symt;
  symt.add_variable("t",arg);
  symt.add_constants();
  re.register_symbol_table(symt);
  im.register_symbol_table(symt);
  parser.compile(re_str,re);
  parser.compile(im_str,im);
  if(exprtk::expression_helper<double>::is_constant(re))
    this->re_str = std::to_string(re());
  if(exprtk::expression_helper<double>::is_constant(im))
    this->im_str = std::to_string(im());
}

/*
time_expr::time_expr(time_expr const& te) :
  _is_real(te._is_real), re_str(te.re_str), im_str(te.im_str), re(te.re), im(te.im), arg(te.arg) {
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
*/
}
