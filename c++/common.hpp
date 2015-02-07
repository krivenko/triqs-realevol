#pragma once

#include <complex>
#include <type_traits>

#include "time_expr_r.hpp"
#include "time_expr_c.hpp"

#include <triqs/operators/many_body_operator.hpp>
#include <triqs/hilbert_space/fundamental_operator_set.hpp>
#include <triqs/hilbert_space/hilbert_space.hpp>
#include <triqs/hilbert_space/imperative_operator.hpp>
#include <triqs/hilbert_space/state.hpp>
#include <triqs/hilbert_space/space_partition.hpp>
#include <triqs/utility/draft/numeric_ops.hpp>

#include "mesh_container.hpp"

namespace realevol {

using dcomplex = std::complex<double>;

using triqs::utility::many_body_operator;
using triqs::hilbert_space::imperative_operator;
using triqs::hilbert_space::state;
using triqs::hilbert_space::fundamental_operator_set;
using triqs::hilbert_space::hilbert_space;
using triqs::hilbert_space::sub_hilbert_space;
using triqs::hilbert_space::space_partition;

template<bool ComplexOp> using operator_coeff_t = typename std::conditional<ComplexOp,time_expr_c,time_expr_r>::type;
template<bool ComplexOp> using operator_t = many_body_operator<operator_coeff_t<ComplexOp>>;


template<bool ComplexOp> using op_on_space_t = imperative_operator<hilbert_space,operator_coeff_t<ComplexOp>,false>;
template<bool ComplexOp> using op_on_subspace_t = imperative_operator<sub_hilbert_space,operator_coeff_t<ComplexOp>,false>;
using state_on_space_t = state<hilbert_space,dcomplex,false>;
using state_on_subspace_t = state<sub_hilbert_space,dcomplex,false>;

template<bool ComplexOp> using space_partition_t = space_partition<state<hilbert_space,operator_coeff_t<ComplexOp>,true>,op_on_space_t<ComplexOp>>;

template<typename Mesh> using results_t = std::map<std::string,mesh_container<double,Mesh>>;

}