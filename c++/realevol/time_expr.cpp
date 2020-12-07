/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
 *
 * realevol is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * realevol is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * realevol. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "time_expr.hpp"

#include <sstream>
#include <utility>

#include <triqs/utility/exceptions.hpp>

namespace realevol {

// Analog of std::to_string with controllable formating
template<typename T> std::string to_string(T x) {
 static std::stringstream ss;
 ss.str("");
 ss << x;
 return ss.str();
}

// Helper methods
exprtk::symbol_table<double> time_expr::create_sym_table() const {
  exprtk::symbol_table<double> symt;
  symt.add_variable("t",arg);
  symt.add_constants();
  return symt;
}

void time_expr::recompile_re_expr() {
  bool res = parser.compile(re_str,re);
  if(!res)
    TRIQS_RUNTIME_ERROR << "ExprTk error in the real part '" << re_str
                        << "': " << parser.error();
  if(exprtk::expression_helper<double>::is_constant(re)) {
    re_str = to_string(re());
    if(re_str == "-0") re_str = "0";
  }
}
void time_expr::recompile_im_expr() {
  bool res = parser.compile(im_str,im);
  if(!res)
    TRIQS_RUNTIME_ERROR << "ExprTk error in the imaginary part '" << im_str
                        << "': " << parser.error();
  if(exprtk::expression_helper<double>::is_constant(im)) {
    im_str = to_string(im());
    if(im_str == "-0") im_str = "0";
  }
}

void time_expr::init_re_expr(exprtk::symbol_table<double> & symt) {
  re.register_symbol_table(symt);
  recompile_re_expr();
}
void time_expr::init_im_expr(exprtk::symbol_table<double> & symt) {
  im.register_symbol_table(symt);
  recompile_im_expr();
}

// Constructors
time_expr::time_expr() : time_expr(.0) {}

time_expr::time_expr(double r) : time_expr(to_string(r)) {}

time_expr::time_expr(double r, double i) : time_expr(to_string(r),to_string(i)) {}

time_expr::time_expr(std::complex<double> const& z) : time_expr(z.real(),z.imag()) {}

time_expr::time_expr(const char* str) : time_expr(std::string(str)) {}

time_expr::time_expr(const char* re_str, const char* im_str) :
  time_expr(std::string(re_str),std::string(im_str)) {}

time_expr::time_expr(std::string const& re_str, double i) : time_expr(re_str,to_string(i)) {}
time_expr::time_expr(double r, std::string const& im_str) : time_expr(to_string(r),im_str) {}
time_expr::time_expr(const char* re_str, double i) : time_expr(std::string(re_str),to_string(i)) {}
time_expr::time_expr(double r, const char* im_str) : time_expr(to_string(r),std::string(im_str)) {}

time_expr::time_expr(std::string str) :
  re_str(std::move(str)), _is_real(true) {
  auto symt = create_sym_table();
  init_re_expr(symt);
}

time_expr::time_expr(std::string re_str, std::string im_str) :
  re_str(std::move(re_str)), im_str(std::move(im_str)), _is_real(false) {
  auto symt = create_sym_table();
  init_re_expr(symt);
  init_im_expr(symt);
}

time_expr::time_expr(time_expr const& te) :
  re_str(te.re_str), im_str(te.im_str), _is_real(te._is_real) {
  auto symt = create_sym_table();
  init_re_expr(symt);
  if(!_is_real) init_im_expr(symt);
}

time_expr::time_expr(time_expr && te) noexcept :
  re_str(std::move(te.re_str)),
  im_str(std::move(te.im_str)),
  _is_real(std::move(te._is_real)) {
  auto symt = create_sym_table();
  init_re_expr(symt);
  if(!_is_real) init_im_expr(symt);
}

std::complex<double> time_expr::operator()(double t) const {
  arg = t;
  return _is_real ? re.value() : std::complex<double>(re.value(),im.value());
}

std::ostream& operator<<(std::ostream& os, time_expr const& te) {
  return te._is_real ? (os << te.re_str) :
                       (os << "(" << te.re_str << "," << te.im_str << ")");
}

time_expr time_expr::operator-() const {
  return _is_real ? time_expr("-(" + re_str + ")") :
                    time_expr("-(" + re_str + ")","-(" + im_str + ")");
}

time_expr & time_expr::operator=(const time_expr& te) {
  re_str = te.re_str;
  recompile_re_expr();
  if(!te._is_real) {
    im_str = te.im_str;
    if(_is_real) {
      auto symt = create_sym_table();
      init_im_expr(symt);
    } else
      recompile_im_expr();
  } else {
    im_str = "";
    im.release();
  }
  _is_real = te._is_real;
  return *this;
}

time_expr & time_expr::operator=(std::string const& str) {
  if(!_is_real) {
    im_str = "";
    im.release();
  }
  _is_real = true;
  re_str = str;
  recompile_re_expr();
  return *this;
}

time_expr & time_expr::operator=(const char* expr) {
  *this = std::string(expr);
  return *this;
}

time_expr & time_expr::operator=(double r) {
  *this = to_string(r);
  return *this;
}

inline std::string sum_str(std::string const& a, std::string const& b) {
 if(a == "0") return b;
 if(b == "0") return a;
 return "(" + a + ")+(" + b + ")";
}

time_expr & time_expr::operator+=(const time_expr& te) {
  re_str = sum_str(re_str, te.re_str);
  recompile_re_expr();
  if(!te._is_real) {
    if(_is_real) {
      im_str = te.im_str;
      auto symt = create_sym_table();
      init_im_expr(symt);
      _is_real = false;
    } else {
      im_str = sum_str(im_str, te.im_str);
      recompile_im_expr();
    }
  }
  return *this;
}

inline std::string unary_minus_str(std::string const& a) {
 if(a == "0") return "0";
 return "-(" + a + ")";
}

inline std::string sub_str(std::string const& a, std::string const& b) {
 if(a == "0" && b == "0") return "0";
 if(a == "0") return "-(" + b + ")";
 if(b == "0") return a;
 return "(" + a + ")-(" + b + ")";
}

time_expr & time_expr::operator-=(const time_expr& te) {
  re_str = sub_str(re_str, te.re_str);
  recompile_re_expr();
  if(!te._is_real) {
    if(_is_real) {
      im_str = unary_minus_str(te.im_str);
      auto symt = create_sym_table();
      init_im_expr(symt);
      _is_real = false;
    } else {
      im_str = sub_str(im_str, te.im_str);
      recompile_im_expr();
    }
  }
  return *this;
}

inline std::string prod_str(std::string const& a, std::string const& b) {
 if(a == "0" || b == "0") return "0";
 if(a == "1") return b;
 if(b == "1") return a;
 return "(" + a + ")*(" + b + ")";
}

time_expr & time_expr::operator*=(const time_expr& te) {
  auto new_re_str = prod_str(re_str, te.re_str);
  if(!_is_real || !te._is_real) {
    if(!_is_real && !te._is_real) {
      new_re_str = sub_str(new_re_str, prod_str(im_str, te.im_str));
      im_str = sum_str(prod_str(re_str, te.im_str), prod_str(im_str, te.re_str));
      recompile_im_expr();
    } else if(te._is_real) {
      im_str = prod_str(im_str, te.re_str);
      recompile_im_expr();
    } else { // _is_real
      im_str = prod_str(re_str, te.im_str);
      auto symt = create_sym_table();
      init_im_expr(symt);
      _is_real = false;
    }
  }
  re_str = new_re_str;
  recompile_re_expr();
  return *this;
}

inline std::string sqr_str(std::string const& a) {
 if(a == "0") return "0";
 if(a == "1") return "1";
 return "(" + a + ")^2";
}

inline std::string div_str(std::string const& a, std::string const& b) {
 if(a == "0" && b != "0") return "0";
 if(b == "1") return a;
 return "(" + a + ")/(" + b + ")";
}

time_expr & time_expr::operator/=(const time_expr& te) {
  auto denom = sqr_str(te.re_str);
  auto new_re_str = prod_str(re_str, te.re_str);
  if(!_is_real || !te._is_real) {
    if(!_is_real && !te._is_real) {
      denom = sum_str(denom, sqr_str(te.im_str));
      new_re_str = sum_str(new_re_str, prod_str(im_str, te.im_str));
      im_str = div_str(sub_str(prod_str(im_str, te.re_str), prod_str(re_str, te.im_str)), denom);
      recompile_im_expr();
    } else if(te._is_real) {
      im_str = div_str(prod_str(im_str, te.re_str), denom);
      recompile_im_expr();
    } else { // _is_real
      denom = sum_str(denom, sqr_str(te.im_str));
      im_str = div_str(unary_minus_str(prod_str(re_str, te.im_str)), denom);
      auto symt = create_sym_table();
      init_im_expr(symt);
      _is_real = false;
    }
  }
  re_str = div_str(new_re_str, denom);
  recompile_re_expr();
  return *this;
}

bool time_expr::operator==(const time_expr& te) const {
  if(_is_real != te._is_real) return false;
  return re_str == te.re_str &&
        (_is_real ? true : (im_str == te.im_str));
}

}
