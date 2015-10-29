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
#pragma once

#include <complex>
#include <string>
#include <ostream>
#include <boost/operators.hpp>

#include <boost/serialization/string.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/concept_check.hpp>

#include <triqs/utility/draft/numeric_ops.hpp>
#include <triqs/utility/mini_vector.hpp>

// Disable some features of ExprTk to speed up compilation
#define exprtk_disable_comments
#define exprtk_disable_break_continue
#define exprtk_disable_string_capabilities
#include "exprtk/exprtk.hpp"

namespace realevol {

using dcomplex = std::complex<double>;

// Time-dependent expressions;
class time_expr : public boost::operators<time_expr>
{
  std::string re_str, im_str;
  exprtk::expression<double> re, im;
  bool _is_real;

  mutable double arg;

public:

  // Constructors
  time_expr();
  time_expr(double r);
  time_expr(double r, double i);
  time_expr(dcomplex const& z);
  time_expr(std::string const& re_str);
  time_expr(const char* str);

  time_expr(std::string const& re_str, std::string const& im_str);
  time_expr(const char* re_str, const char* im_str);
  time_expr(std::string const& re_str, double i);
  time_expr(double r, std::string const& im_str);
  time_expr(const char* re_str, double i);
  time_expr(double r, const char* im_str);

  time_expr(time_expr const&);
  time_expr(time_expr &&) = default;

  // Assignments
  time_expr & operator=(double r);
  time_expr & operator=(std::string const& str);
  time_expr & operator=(const char* str);
  time_expr & operator=(time_expr const& te);

  // Arithmetics
  time_expr operator-() const;

  time_expr & operator+=(time_expr const& te);
  time_expr & operator-=(time_expr const& te);
  time_expr & operator*=(time_expr const& te);
  time_expr & operator/=(time_expr const& te);
  bool operator==(time_expr const& te) const;

  // Is this expression real?
  bool is_real() const { return _is_real; }

  // Evaluation of the expression
  dcomplex operator()(double t) const;

  friend bool is_constant(time_expr const& te) {
    return exprtk::expression_helper<double>::is_constant(te.re) &&
           (te._is_real ? true : exprtk::expression_helper<double>::is_constant(te.im));
  }
  bool is_zero() const {
    return is_constant(*this) &&
           triqs::utility::is_zero(re.value()) &&
           (_is_real ? true : triqs::utility::is_zero(im.value()));
  }

  // Stream output
  friend std::ostream& operator<<(std::ostream& os, time_expr const& te);

private:

  inline exprtk::symbol_table<double> create_sym_table() const;
  inline void init_re_expr(std::string const& str, exprtk::symbol_table<double> & symt);
  inline void init_im_expr(std::string const& str, exprtk::symbol_table<double> & symt);

  // Methods for boost::serialization
  friend class boost::serialization::access;
/*
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const {
    ar & _is_real;
    ar & TRIQS_MAKE_NVP("re_str",re_str);
    if(!_is_real) ar & TRIQS_MAKE_NVP("im_str",im_str);
  }

  template<class Archive>
  void load(Archive & ar, const unsigned int version) {
    // TODO
    ar & _is_real;
    std::string str;

    if(_is_real) {
      ar & str;
      *this = time_expr(str);
    } else {
      ar & str;
    }
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()*/
};

inline time_expr operator ""_te (long double r){ return time_expr(r); }
inline time_expr operator ""_te(const char* expr, std::size_t) { return time_expr(expr); };

/*
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
}*/

}

namespace triqs { namespace utility {

inline bool is_zero(realevol::time_expr const& te, double = 0 /* neglected */) { return te.is_zero(); }
inline realevol::time_expr conj(realevol::time_expr const& te) { return te; }

}}
