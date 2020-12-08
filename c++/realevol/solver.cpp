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

// TODO: Remove this file altogether

#include <triqs/utility/first_include.hpp>

#include <array>
#include <csignal>
#include <set>
#include <algorithm>
#include <iterator>
#include <limits>
#include <vector>

#include <triqs/utility/signal_handler.hpp>

#include "array_utility.hpp"
#include "time_expr.hpp"
#include "time_interp.hpp"
#include "solver.hpp"
#include "hs_structure.hpp"
#include "init_state.hpp"
#include "worldlines.hpp"
#include "mpi_dispatcher.hpp"
#include "wl_worker.hpp"
#include "mesh_utils.hpp"
#include "time_point_selector.hpp"

namespace signal_handler = triqs::signal_handler;

static void check_signals() {
 if(signal_handler::received()) {
  int signal = signal_handler::last();
  if(signal == SIGINT) TRIQS_KEYBOARD_INTERRUPT;
  else                 TRIQS_RUNTIME_ERROR;
 }
}

namespace realevol {

const double nan = std::numeric_limits<double>::quiet_NaN();

solver::solver(gf_struct_t const& gf_struct, chi_indices_t const& chi_indices,
               double t_max, int n_t) :
 gf_struct(gf_struct), chi_indices(chi_indices),
 t_mesh(0, t_max, n_t) {

 if(t_max <= 0) TRIQS_RUNTIME_ERROR << "t_max must be positive";

 struct index_visitor  {
  std::vector<std::string> indices;
  void operator()(int i) { indices.push_back(std::to_string(i)); }
  void operator()(std::string s) { indices.push_back(s); }
 };

 std::vector<std::string> block_names;
 std::vector<gf_2t_t> g_l_blocks, g_g_blocks;
 for (auto const& bl : gf_struct) {
  block_names.push_back(bl.first);
  int n = bl.second.size();

  index_visitor iv;
  for (auto & ind: bl.second) { std::visit(iv, ind); }
  std::vector<std::vector<std::string>> indices{{iv.indices,iv.indices}};

  g_l_blocks.push_back(gf_2t_t{{t_mesh, t_mesh}, make_shape(n, n), indices});
  g_g_blocks.push_back(gf_2t_t{{t_mesh, t_mesh}, make_shape(n, n), indices});
 }

 g_l = make_block_gf(block_names, g_l_blocks);
 g_g = make_block_gf(block_names, g_g_blocks);

 g_l() = dcomplex(nan, nan);
 g_g() = dcomplex(nan, nan);

 int chi_size = chi_indices.size();
 if(chi_size != std::set<std::pair<std::string, std::variant<int, std::string>>>
   (chi_indices.begin(), chi_indices.end()).size())
  TRIQS_RUNTIME_ERROR << "Repeated indices are met in chi_indices";

 if(chi_size > 0) {
  chi = gf_2t_t{{t_mesh, t_mesh}, {chi_size, chi_size}};
  chi() = dcomplex(nan, nan);
 }
}

template<typename HamiltonianType>
HamiltonianType try_reduce_to_constant(HamiltonianType const& op, gf_mesh<retime> const& mesh) {
 HamiltonianType res;
 for(auto const& m : op) {
  auto new_coef = try_reduce_to_constant<typename HamiltonianType::scalar_t, gf_mesh<retime>>(m.coef, mesh);
  HamiltonianType new_monomial(new_coef);
  for(auto const& c : m.monomial)
   new_monomial *= HamiltonianType::make_canonical(c.stat, c.dagger, c.indices);
  res += new_monomial;
 }
 return res;
}

void init_observable(gf_2t_view f, time_point_selector<2> const& t_selector) {
 gf_mesh<retime>::mesh_point_t t, tp;
 for(auto ttp : f.mesh()) {
  std::tie(t, tp) = ttp.components_tuple();
  f[{t, tp}]() = t_selector({t, tp}) ? .0 : dcomplex(nan, nan);
 }
}

void restore_antihermiticity(gf_2t_view f, time_point_selector<2> const& t_selector) {
 gf_mesh<retime>::mesh_point_t t, tp;

 // FIXME: Workaround for TRIQS issue #798
 matrix<dcomplex> mat(f.mesh().size(), f.mesh().size());
 for(auto ttp : f.mesh()) {
  std::tie(t, tp) = ttp.components_tuple();
  if(t_selector({t, tp})) { // Have we computed this pair?
   mat = f[{t, tp}];
   f[{tp, t}] = -dagger(mat);
  }
 }
}

template<typename HamiltonianType>
void solver::compute_2t_obs(HamiltonianType const& h_, compute_2t_obs_parameters_t const& params) {

 // Save parameters
 compute_2t_obs_params = params;

 if(!(params.compute_g_l || params.compute_g_g || params.compute_chi)) {
  if(params.verbosity >= 1 && comm.rank() == 0)
   std::cout << "Nothing to compute. Exiting." << std::endl;
  return;
 }

 // Do we have a valid initial state to start computation?
 if(initial_state == nullptr)
  TRIQS_RUNTIME_ERROR << "Initial state has not been set!";

 // Complete fundamental operator set of the problem
 auto const& fops = initial_state->get_fops();

 auto check_fops_subset = [&fops](fundamental_operator_set const& subfops,
                                  statistic_enum stat, std::string const& error) {
  for(auto it = subfops.begin(stat); it != subfops.end(stat); ++it) {
   if(!fops.has_indices(it->index,stat)) TRIQS_RUNTIME_ERROR << error
   << it->index << " incompatible with the initial problem";
  }
 };

 signal_handler::start();

 // Fundamental operator set generated by GFs
 auto gf_fops = (params.compute_g_l || params.compute_g_g) ?
                fundamental_operator_set(gf_struct, {}) : fundamental_operator_set();
 check_fops_subset(gf_fops, Fermion, "Green's function contains index ");

 // Fundamental operator set generated by susceptibilities
 fundamental_operator_set chi_fops;
 if(params.compute_chi) {
  for(auto const& ind : chi_indices) chi_fops.insert_fermion(ind.first, ind.second);
 }
 check_fops_subset(chi_fops, Fermion, "Susceptibility contains index ");

 // Fundamental operator set generated by the Hamiltonian
 auto h_fops = h_.make_fundamental_operator_set();
 check_fops_subset(h_fops, Fermion, "Hamiltonian contains fermionic index ");
 check_fops_subset(h_fops, Boson, "Hamiltonian contains bosonic index ");

 // Scan for the coefficients constant at all points of t_mesh,
 // and replace them with the corresponding constants
 auto h = try_reduce_to_constant(h_, t_mesh);

 if(params.verbosity >= 2)
  std::cout << "The Hamiltonian of the problem:" << std::endl << h << std::endl;

 // Fundamental operator set generated by all observables
 auto obs_fops = merge(gf_fops, chi_fops);

 // Analyse structure of the Hilbert space
 hilbert_space_structure<HamiltonianType> hs_struct(h,
                                                    initial_state->get_fops(),
                                                    initial_state->get_full_hs(),
                                                    obs_fops,
                                                    is_zero_on_mesh<gf_mesh<retime>>(t_mesh));
 check_signals();

 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Found " << hs_struct.sub_hilbert_spaces.size()
            << " invariant subspaces." << std::endl;
 if(params.verbosity >= 2 && comm.rank() == 0) {
  for(auto const& sp : hs_struct.sub_hilbert_spaces) {
  std::cout << " Subspace " << sp.get_index() << " [" << sp.size() << "]: ";
  for(auto f : sp.get_all_fock_states()) std::cout << f << " ";
   std::cout << std::endl;
  }
 }

 // Check on which invariant subspaces the Hamiltonian is static
 auto static_pred = [](auto const& st) {
   bool c = true;
   foreach(st,[&c](uint32_t, auto const& te){ if(!is_constant(te)) c = false; });
   return c;
 };

 auto is_static_sp = hs_struct.classify_subspaces(h, static_pred);

 if(params.verbosity >= 2 && comm.rank() == 0) {
  std::cout << "The Hamiltonian is static on the following invariant subspaces";
  std::cout << std::endl << " ";
  for(long spn = 0; spn < is_static_sp.size(); ++spn) {
   if(is_static_sp[spn]) std::cout << spn << " ";
  }
  std::cout << std::endl;
 }

 // Compute subspace branching for the initial state
 auto const& subspaces = initial_state->get_sub_hilbert_spaces();
 auto subspace_branchings = hs_struct.compute_branchings(subspaces);
 check_signals();

 if(params.verbosity >= 2 && comm.rank() == 0) {
  std::cout << "Subspace branching for the initial state:" << std::endl;
  for(long spn = 0; spn < subspaces.size(); ++spn) {
   std::cout << " " << subspaces[spn].get_index() << " -> ";
   for(long spi : subspace_branchings[spn]) std::cout << spi << " ";
   std::cout << std::endl;
  }
 }

 // Generate all contributing world lines
 worldlines_maker<HamiltonianType> wlm(*initial_state,
                                       hs_struct,
                                       subspace_branchings,
                                       params.hbar);

 auto g_g_wl = params.compute_g_g ?
               wlm.make_gf_worldlines(gf_struct, true) : std::vector<worldline_desc_t<2>>();
 check_signals();
 auto g_l_wl = params.compute_g_l ?
               wlm.make_gf_worldlines(gf_struct, false) : std::vector<worldline_desc_t<2>>();
 check_signals();
 auto chi_wl = params.compute_chi ?
               wlm.make_chi_worldlines(chi_indices) : std::vector<worldline_desc_t<2>>();
 check_signals();

 long nwl = 0;
 if(params.verbosity >= 2 && comm.rank() == 0) {
  auto print_worldlines = [&nwl](auto const& worldlines, std::string const& msg) {
   std::cout << "World lines contributing to " << msg << " ("
             << worldlines.size() << " in total)" << std::endl;
   for(long i : range(worldlines.size())) {
    std::cout << "[" << nwl << "] " << worldlines[i] << std::endl;
    ++nwl;
   }
  };
  if(params.compute_g_g) print_worldlines(g_g_wl, "the greater GF component:");
  if(params.compute_g_l) print_worldlines(g_l_wl, "the lesser GF component:");
  if(params.compute_chi) print_worldlines(chi_wl, "the susceptibility:");
 }

 std::vector<worldline_desc_t<2>> all_worldlines;
 all_worldlines.reserve(g_g_wl.size() + g_l_wl.size() + chi_wl.size());
 for(auto && wl : {g_g_wl, g_l_wl, chi_wl})
  std::move(wl.begin(), wl.end(), std::back_inserter(all_worldlines));
 if(all_worldlines.size() == 0) return;

 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Starting GF and susceptibility calculation" << std::endl;

 mpi_dispatcher<long> disp(comm, all_worldlines.size());

 time_point_selector<2> t_selector({params.t_range, params.tp_range},
   std::array<double, 1>{params.delta_t_max},
   true
 );

 wl_worker<HamiltonianType> worker(*initial_state, h, params.hbar, hs_struct, is_static_sp,
                                   t_selector,
                                   params.lanczos_min_matrix_size,
                                   params.lanczos_gs_energy_tol,
                                   params.lanczos_max_krylov_dim
                                  );

 auto choose_obs =
 [this](worldline_desc_t<2> const& wl, bool verbose, long nwl) -> gf_2t_t& {
  if(verbose)
   std::cout << "[Node " << comm.rank() << "] Evaluating world line " << nwl;
  switch(wl.observable) {
   case worldline_desc_t<2>::GreaterGf:
    if(verbose) std::cout << " (greater GF component)" << std::endl;
    return g_g[wl.block_index];
   case worldline_desc_t<2>::LesserGf:
    if(verbose) std::cout << " (lesser GF component)" << std::endl;
    return g_l[wl.block_index];
   case worldline_desc_t<2>::Susceptibility:
    if(verbose) std::cout << " (susceptibility component)" << std::endl;
    return chi;
  }
 };

 check_signals();

 for(auto & f : g_g) {
  if(params.compute_g_g)
   init_observable(f, t_selector);
  else
   f() = dcomplex(nan, nan);
 }
 for(auto & f : g_l) {
  if(params.compute_g_l)
   init_observable(f, t_selector);
  else
   f() = dcomplex(nan, nan);
 }
 if(params.compute_chi)
  init_observable(chi, t_selector);
 else
  chi() = dcomplex(nan, nan);

 while(true) {
  if((nwl = disp().value_or(-1)) == -1) break;
  auto const& wl = all_worldlines[nwl];
  auto & obs = choose_obs(wl, params.verbosity >= 2, nwl);
  switch(params.hamiltonian_interpol) {
   case Rectangle:
    worker(wl, obs, std::integral_constant<h_interpolation,Rectangle>());
    break;
   case Trapezoid:
    worker(wl, obs, std::integral_constant<h_interpolation,Trapezoid>());
    break;
   case Simpson:
    worker(wl, obs, std::integral_constant<h_interpolation,Simpson>());
    break;
  }
  check_signals();
 }
 comm.barrier();

 check_signals();

 // Collect results from all MPI ranks
 if(params.compute_g_g) g_g = mpi_reduce(g_g, comm, 0, true);
 if(params.compute_g_l) g_l = mpi_reduce(g_l, comm, 0, true);
 if(params.compute_chi) chi = mpi_reduce(chi, comm, 0, true);

 if(params.compute_g_g)
  for(auto & f : g_g) restore_antihermiticity(f, t_selector);
 if(params.compute_g_l)
  for(auto & f : g_l) restore_antihermiticity(f, t_selector);
 if(params.compute_chi)
  restore_antihermiticity(chi, t_selector);

 signal_handler::stop();

 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Done" << std::endl;
}

template
void solver::compute_2t_obs(time_expr_operator_t const& h, compute_2t_obs_parameters_t const& params);
template
void solver::compute_2t_obs(time_interp_operator_t const& h, compute_2t_obs_parameters_t const& params);

} // namespace realevol
