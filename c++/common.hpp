#pragma once

#include <complex>
#include <type_traits>

#include "time_expr_r.hpp"
#include "time_expr_c.hpp"

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>
#include <triqs/draft/hilbert_space_tools/hilbert_space.hpp>
#include <triqs/draft/hilbert_space_tools/imperative_operator.hpp>
#include <triqs/draft/hilbert_space_tools/state.hpp>
#include <triqs/utility/draft/numeric_ops.hpp>

#include "mesh_container.hpp"
#include "space_partition.hpp"

namespace realevol {

using dcomplex = std::complex<double>;

using triqs::utility::many_body_operator;
using triqs::utility::imperative_operator;
using triqs::utility::state;
using triqs::utility::fundamental_operator_set;
using triqs::utility::hilbert_space;
using triqs::utility::sub_hilbert_space;

template<bool ComplexOp> using operator_coeff_t = typename std::conditional<ComplexOp,time_expr_c,time_expr_r>::type;
template<bool ComplexOp> using operator_t = many_body_operator<operator_coeff_t<ComplexOp>>;


template<bool ComplexOp> using op_on_space_t = imperative_operator<hilbert_space,operator_coeff_t<ComplexOp>,false>;
template<bool ComplexOp> using op_on_subspace_t = imperative_operator<sub_hilbert_space,operator_coeff_t<ComplexOp>,false>;
using state_on_space_t = state<hilbert_space,dcomplex,false>;
using state_on_subspace_t = state<sub_hilbert_space,dcomplex,false>;

template<bool ComplexOp> using space_partition_t = space_partition<state<hilbert_space,operator_coeff_t<ComplexOp>,true>,op_on_space_t<ComplexOp>>;

template<typename Mesh> using results_t = std::map<std::string,mesh_container<double,Mesh>>;

}