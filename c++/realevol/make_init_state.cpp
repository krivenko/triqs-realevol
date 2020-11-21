/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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
#include <triqs/utility/first_include.hpp>

#ifndef NDEBUG
#define TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
#endif

#include <cmath>
#include <limits>
#include <vector>
#include <set>
#include <array>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <type_traits>

#include <signal.h>

#include <triqs/arrays/vector.hpp>
#include <triqs/arrays/linalg/eigenelements.hpp>
#include <triqs/utility/signal_handler.hpp>
#include <mpi/mpi.hpp>

#include <realevol/hilbert_space/space_partition.hpp>
#include <realevol/hilbert_space/state_view.hpp>

#include <ezarpack/storages/triqs.hpp>
#include <ezarpack/arpack_solver.hpp>
#include <ezarpack/version.hpp>

#include "init_state.hpp"
#include "mpi_dispatcher.hpp"
#include "global_index.hpp"
#include "sort2.hpp"

namespace signal_handler = triqs::signal_handler;

#define CHECK_SIGNALS                           \
if(signal_handler::received()) {                \
 int signal = signal_handler::last();           \
 if(signal == SIGINT) TRIQS_KEYBOARD_INTERRUPT; \
 else                 TRIQS_RUNTIME_ERROR;      \
}

namespace realevol {

template<typename OperatorType>
static_operator_t make_static_op(OperatorType const& op, std::string const& error_message) {
 static_operator_t static_op;
 for(auto const& m : op) {
  if(!is_constant(m.coef)) TRIQS_RUNTIME_ERROR << error_message;
  auto new_coef = m.coef(0);
  static_operator_t new_monomial(new_coef);
  for(auto const& c : m.monomial)
   new_monomial *= static_operator_t::make_canonical(c.stat, c.dagger, c.indices);
  static_op += new_monomial;
 }
 return static_op;
}

template<typename OperatorType>
init_state make_pure_init_state(OperatorType const& generator,
                                fundamental_operator_set const& fops,
                                std::map<operators::indices_t, int> const& bits_per_boson) {

 // Static version of generator
 auto gen = make_static_op(generator, "Generating operator must be time-independent!");

 signal_handler::start();

 init_state ist(fops, bits_per_boson);

 // Add one subspace with index 0
 ist.sub_hilbert_spaces.emplace_back(0);
 auto & sp = ist.sub_hilbert_spaces[0];

 // Vaccum state in the full Hilbert space
 state_on_space_t vacuum(ist.full_hs);
 vacuum(fock_state_t(0)) = 1;

 // Imperative version of generator
 static_op_on_space_t imp_gen(gen, ist.fops, ist.full_hs);

 auto psi = imp_gen(vacuum);

 // Fill st.sub_hilbert_spaces[0]
 foreach(psi, [&sp](fock_state_t f, dcomplex) { sp.add_fock_state(f); });
 CHECK_SIGNALS;

 // Add a weighted state
 state_on_subspace_t st_on_sp(sp);
 double norm2 = 0;
 foreach(psi, [&st_on_sp,&sp,&norm2](fock_state_t f, dcomplex a) {
  st_on_sp(sp.get_state_index(f)) = a;
  norm2 += abs(a*a);
 });
 if(norm2==0) TRIQS_RUNTIME_ERROR << "Generated state is trivial!";
 st_on_sp /= std::sqrt(norm2);

 ist.weighted_states.emplace_back(std::move(st_on_sp), 1.0);

 signal_handler::stop();

 return ist;
}

template
init_state make_pure_init_state(time_expr_operator_t const&,
                                fundamental_operator_set const&,
                                std::map<operators::indices_t, int> const&);
template
init_state make_pure_init_state(time_interp_operator_t const&,
                                fundamental_operator_set const&,
                                std::map<operators::indices_t, int> const&);

// -----------------------------------------------------------------

using namespace triqs::arrays;
using sp_levels_t = std::pair<int,vector<double>>;
using eigensystem_t = std::pair<vector<double>,matrix<dcomplex>>;
using real_state_on_subspace_t = state<sub_hilbert_space,double,false>;
using real_static_op_on_subspace_t = imperative_operator<sub_hilbert_space,double,false>;

/////////////////////////////////////////////////////////
/// Solver procedures for make_equilibrium_init_state ///
/////////////////////////////////////////////////////////

auto make_arpack_solver_params(int nev, bool eigenvectors, std::true_type) {
 using namespace ezarpack;
 using params_t = arpack_solver<Complex, triqs_storage>::params_t;
 return params_t(nev, params_t::SmallestReal, eigenvectors ? params_t::Ritz : params_t::None);
}

auto make_arpack_solver_params(int nev, bool eigenvectors, std::false_type) {
 using namespace ezarpack;
 using params_t = arpack_solver<Symmetric, triqs_storage>::params_t;
 return params_t(nev, params_t::Smallest, eigenvectors);
}

// -----------------------------------------------------------------

template<typename T>
void find_lowest_levels_on_subspace(sub_hilbert_space const& sp,
                                    imperative_operator<sub_hilbert_space,T,false> const& h,
                                    std::vector<sp_levels_t> & lowest_levels,
                                    double & gs_energy,
                                    double energy_window,
                                    int verbosity,
                                    int arpack_min_matrix_size,
                                    int arpack_tol, int arpack_ncv,
                                    int rank) {

 using triqs::utility::real;
 using is_complex_t = std::is_same<T,dcomplex>;

 int sp_n = sp.get_index();
 int N = sp.size();

 if(N == 1) { // 1-dimensional subspace -> no need to diagonalize

  if(verbosity >= 2)
   std::cout << "[Node " << rank << "] Subspace "
             << sp_n << " is trivial" << std::endl;

  state<sub_hilbert_space,T,false> st(sp);
  st(0) = 1.0;
  double energy = real(h(st)(0));
  gs_energy = std::min(gs_energy, energy);
  if(energy <= gs_energy + energy_window)
   lowest_levels.emplace_back(sp_levels_t{sp_n, {energy}});
  return;

 } else if(sp.size() < arpack_min_matrix_size) { // Small problem -> use LAPACK

  if(verbosity >= 2)
   std::cout << "[Node " << rank << "] Using LAPACK on subspace "
             << sp_n << std::endl;

  matrix<T> M(N,N);
  state<sub_hilbert_space,T,false> from(sp), to(sp);
  for(long i : range(N)) {
   from(i) = 1;
   h.apply(from, to);
   M(range(),i) = to.amplitudes();
   from(i) = 0;
  }

  auto eig = linalg::eigenvalues_in_place(&M);
  gs_energy = std::min(gs_energy, eig(0));
  int i = 0;
  for(; i < N && eig(i) <= gs_energy + energy_window; ++i);
  if(i > 0) lowest_levels.emplace_back(sp_levels_t{sp_n, eig(range(i))});
  return;
 }

 if(verbosity >= 2)
  std::cout << "[Node " << rank << "] Using ARPACK on subspace "
            << sp_n << std::endl;

 // Set up ARPACK solver
 using namespace ezarpack;
 arpack_solver<is_complex_t::value ? Complex : Symmetric,
               triqs_storage> arps(sp.size());

 using state_view_t = state_view<sub_hilbert_space,T>;
 std::array<state_view_t,3> state_views = {
  state_view_t{arps.workspace_vector(0), sp},
  state_view_t{arps.workspace_vector(1), sp},
  state_view_t{arps.workspace_vector(2), sp}
 };

 auto params = make_arpack_solver_params(0, false, is_complex_t());
 params.tolerance = arpack_tol;
 params.ncv = arpack_ncv;

 // Partially diagonalize op within sp iteratively increasing
 // the number of eigenvalues to be found

 // Temporary GS energy, reset after each call to the ARPACK solver
 double gs_energy_tmp = gs_energy;
 while(true) {
  ++params.n_eigenvalues;
  if(verbosity >= 2)
  std::cout << "[Node " << rank << "] Calling ARPACK to find the lowest "
            << params.n_eigenvalues << " eigenvalues" << std::endl;

  auto apply_h = [&h,&state_views,&arps](auto, auto) {
   auto in_n = arps.in_vector_n();
   auto out_n = arps.out_vector_n();
   h.apply(state_views[in_n], state_views[out_n]);
  };

  // Run ARPACK solver
  arps(apply_h, params);

  auto const& eig = arps.eigenvalues();

  gs_energy_tmp = std::min(gs_energy, real(eig(params.n_eigenvalues-1)));

  // ARPACK cannot calculate the last eigenvalue. Stopping ...
  if(params.n_eigenvalues == sp.size()-(is_complex_t::value ? 2 : 1)) {
   std::cout << "[Node " << rank << "] WARNING: ARPACK has found "
             << params.n_eigenvalues
             << " relevant eigenvalues on subspace "
             << sp.get_index() << ". "
             << (is_complex_t::value ? "Two eigenvalues" : "One eigenvalue")
             << " can still be missing" << std::endl;
   break;
  }

  // The most recently found eigenvalue was too high
  if(real(eig(0)) > gs_energy_tmp + energy_window) break;
 }

 gs_energy = std::min(gs_energy, gs_energy_tmp); // Use the GS energy from the last iteration

 if(params.n_eigenvalues == 1) return;

 auto const& eig = arps.eigenvalues();
 vector<double> ev(eig.size()-1);
 for(int i : range(ev.size())) ev(i) = real(eig(ev.size()-i));
 lowest_levels.emplace_back(sp_levels_t{sp_n, std::move(ev)});
}

// -----------------------------------------------------------------

template<typename T>
void permute_matrix_cols(matrix<T> & M, std::vector<std::size_t> const& p) {
 int N = second_dim(M);
 assert(N == p.size());
 std::vector<bool> done(N);
 for (std::size_t i = 0; i < N; ++i) {
  if (done[i]) continue;
  done[i] = true;
  std::size_t prev_j = i;
  std::size_t j = p[i];
  while (i != j) {
   using std::swap;
   for(int r : range(first_dim(M)))
    swap(M(r,prev_j), M(r,j));
   done[j] = true;
   prev_j = j;
   j = p[j];
  }
 }
}

// -----------------------------------------------------------------

template<typename T>
void compute_eigenvectors(sub_hilbert_space const& sp,
                          imperative_operator<sub_hilbert_space,T,false> const& h,
                          int n_vectors_to_compute,
                          std::vector<eigensystem_t> & eigensystems,
                          int verbosity,
                          int arpack_min_matrix_size,
                          int arpack_tol, int arpack_ncv,
                          int rank) {

 using triqs::utility::real;
 using is_complex_t = std::is_same<T,dcomplex>;

 int sp_n = sp.get_index();
 int N = sp.size();

 if(N == 1) { // 1-dimensional subspace -> no need to diagonalize

  if(verbosity >= 2)
   std::cout << "[Node " << rank << "] Subspace "
             << sp_n << " is trivial" << std::endl;

  state<sub_hilbert_space,T,false> st(sp);
  st(0) = 1.0;
  double energy = real(h(st)(0));
  eigensystems.emplace_back(vector<double>{energy}, matrix<dcomplex>{{1.0}});
  return;

 } else if(sp.size() < arpack_min_matrix_size) { // Small problem -> use LAPACK

  if(verbosity >= 2)
   std::cout << "[Node " << rank << "] Using LAPACK on subspace "
             << sp_n << " (" << n_vectors_to_compute
             << " eigenpairs to compute)"  << std::endl;

  matrix<T> M(N,N);
  state<sub_hilbert_space,T,false> from(sp), to(sp);
  for(long i : range(N)) {
   from(i) = 1;
   h.apply(from, to);
   M(range(),i) = to.amplitudes();
   from(i) = 0;
  }

  auto eig = linalg::eigenelements_in_place(&M);
  auto r = range(n_vectors_to_compute);
  eigensystems.emplace_back(eig.first(r), eig.second(r,range()).transpose());
  return;

 }

 if(verbosity >= 2)
  std::cout << "[Node " << rank << "] Using ARPACK on subspace "
            << sp_n << " (" << n_vectors_to_compute
            << " eigenpairs to compute)" << std::endl;

 // Set up ARPACK solver
 using namespace ezarpack;
 arpack_solver<is_complex_t::value ? Complex : Symmetric, triqs_storage> arps(sp.size());

 using state_view_t = state_view<sub_hilbert_space,T>;
 std::array<state_view_t,3> state_views = {
  state_view_t{arps.workspace_vector(0), sp},
  state_view_t{arps.workspace_vector(1), sp},
  state_view_t{arps.workspace_vector(2), sp}
 };

 auto params = make_arpack_solver_params(n_vectors_to_compute, true, is_complex_t());
 params.tolerance = arpack_tol;
 params.ncv = arpack_ncv;

 auto apply_h = [&h,&state_views,&arps](auto, auto) {
  auto in_n = arps.in_vector_n();
  auto out_n = arps.out_vector_n();
  h.apply(state_views[in_n], state_views[out_n]);
 };

 // Run ARPACK solver
 arps(apply_h, params);

 eigensystems.emplace_back(real(arps.eigenvalues()),
                           arps.eigenvectors()(range(),range(n_vectors_to_compute)));

 // Sort eigenvalues and eigenvectors as complex ARPACK is not trustworthy in this respect...
 if(is_complex_t::value) {
  auto & es = eigensystems.back();
  auto p = sort_permutation(es.first);
  apply_permutation_in_place(es.first, p);
  permute_matrix_cols(es.second, p);
 }
}

// -----------------------------------------------------------------

struct JobWithGsEnergy {
 long index;
 double gs_energy;
};

} //namespace realevol

namespace mpi {

template<> struct mpi_type<realevol::JobWithGsEnergy> {
 static MPI_Datatype get() noexcept {
  static bool type_committed = false;
  static MPI_Datatype dt;
  if(!type_committed) {
    int blocklengths[] = {1,1};
    MPI_Aint displacements[] = {0,sizeof(long)};
    MPI_Datatype types[] = {MPI_LONG,MPI_DOUBLE};
    MPI_Type_create_struct(2, blocklengths, displacements, types, &dt);
    MPI_Type_commit(&dt);
    type_committed = true;
  }
  return dt;
 }
};

} // namespace mpi

// -----------------------------------------------------------------

namespace realevol {

template<typename OperatorType>
init_state make_equilibrium_init_state(OperatorType const& h,
                                       fundamental_operator_set const& fops,
                                       double temperature,
                                       eq_solver_parameters_t const& params,
                                       std::map<operators::indices_t, int> const& bits_per_boson,
                                       mpi::communicator const& comm) {

 // Static version the equilibrium Hamiltonian
 auto h_ = make_static_op(h, "Initial Hamiltonian must be time-independent!");

 // Check that h is Hermitian
 if(!(h_ - dagger(h_)).is_zero())
  TRIQS_RUNTIME_ERROR << "Supplied Hamitonian is not Hermitian!";

 if(params.arpack_min_matrix_size < 4)
  TRIQS_RUNTIME_ERROR << "arpack_min_matrix_size must be >= 4!";

 double beta, energy_window;
 if(temperature == 0)
  // XXX If this is too low, n_relevant_ev is "0" for zero temperature with ARPACK and MPI.
  // The minimum tolerance that worked on my machine was "1e-13", for "1e-14" there was this problem. mdanilov
  energy_window = (params.arpack_tolerance != 0 ? params.arpack_tolerance : 1e-12);
 else {
  beta = 1 / temperature;
  energy_window = -temperature * std::log(params.min_rel_weight);
 }

 signal_handler::start();

 // Initial state object
 init_state ist(fops, bits_per_boson);

 CHECK_SIGNALS;

 // Partition the Hilbert space
 space_partition<state_on_space_t, static_op_on_space_t>
 partition(state_on_space_t(ist.get_full_hs()), static_op_on_space_t(h_, fops, ist.get_full_hs()));

 CHECK_SIGNALS;

 // Fill subspaces
 std::vector<sub_hilbert_space> subspaces;
 subspaces.reserve(partition.n_subspaces());
 for (long n = 0; n < partition.n_subspaces(); ++n) subspaces.emplace_back(n);
 foreach(partition, [&](fock_state_t s, int spn) { subspaces[spn].add_fock_state(s); });

 CHECK_SIGNALS;

 bool h_is_real = imag(h_).is_zero();
 if(h_is_real) h_ = real(h_);

 if(params.verbosity >= 1 && comm.rank() == 0) {
  std::cout << "Found " << subspaces.size() << " invariant subspaces." << std::endl;
  if(params.verbosity >= 2) {
   std::cout << "Dimensions of subspaces: ";
   for(auto const& sp : subspaces) std::cout << sp.size() << " ";
   std::cout << std::endl;
  }
  std::cout << "Using " << (h_is_real ? "real" : "complex") <<
               " arithmetics in diagonalization." << std::endl;
 }

 // Diagonalization phase 1
 // Distribute subspaces between MPI ranks, and find the lowest energy levels
 // in each subspace without computing the corresponding eigenvectors (to save memory).

 // For solvers: ground state energy on this MPI rank
 // For master: current ground state energy reported by the solvers
 double gs_energy = std::numeric_limits<double>::infinity();

 // Will be called only on master
 auto job_generator = [&gs_energy](long i, double new_gs_energy) {
  gs_energy = std::min(gs_energy, new_gs_energy);
  return JobWithGsEnergy{i, gs_energy};
 };
 mpi_dispatcher<JobWithGsEnergy, double> disp(comm, job_generator, subspaces.size());

 // Rank-local list of processed subspaces and their lowest levels
 std::vector<sp_levels_t> sp_lowest_levels;

 comm.barrier();
 while(true) {
  auto job_or_stop = disp(gs_energy);
  if(!job_or_stop) break;
  auto job = job_or_stop.value();
  auto const& sp = subspaces[job.index];
  gs_energy = job.gs_energy;

  if(params.verbosity >= 2)
   std::cout << "[Node " << comm.rank() << "] Searching the lowest eigenvalues on subspace "
             << sp.get_index() << " (dimension " << sp.size() << "), "
             << "starting E_gs = " << gs_energy << std::endl;

  auto ncv_it = params.arpack_ncv.find(sp.get_index());
  int ncv = ncv_it == params.arpack_ncv.end() ? -1 : ncv_it->second;

  if(h_is_real) {
   real_static_op_on_subspace_t op(h_, fops, ist.get_full_hs());
   find_lowest_levels_on_subspace(sp, op, sp_lowest_levels,
                                  gs_energy, energy_window,
                                  params.verbosity,
                                  params.arpack_min_matrix_size,
                                  params.arpack_tolerance, ncv,
                                  comm.rank());
  } else {
   static_op_on_subspace_t op(h_, fops, ist.get_full_hs());
   find_lowest_levels_on_subspace(sp, op, sp_lowest_levels,
                                  gs_energy, energy_window,
                                  params.verbosity,
                                  params.arpack_min_matrix_size,
                                  params.arpack_tolerance, ncv,
                                  comm.rank());
  }

  CHECK_SIGNALS;

  if(params.verbosity >= 2 && sp_lowest_levels.back().first == sp.get_index())
   std::cout << "[Node " << comm.rank() << "] Provisionally relevant levels on subspace "
             << sp.get_index() << ": "
             << sp_lowest_levels.back().second << std::endl;
 }

 // Find the global energy minimum
 gs_energy = mpi::all_reduce(gs_energy, comm, 0, MPI_MIN);

 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Ground state energy: " << gs_energy << std::endl;

 // Diagonalization phase 2
 // Now we know the exact ground state energy, and can properly select
 // the relevant subspaces. On those we compute the eigenvectors.
 double max_energy = gs_energy + energy_window;

 // Rank-local quantities
 double Z = 0;                               // Partition function
 std::vector<global_index> rel_sp_i;         // List of relevant subspace indices
 std::vector<eigensystem_t> eigensystems;    // List of eigensystems

 int rank_local_n = std::accumulate(sp_lowest_levels.begin(), sp_lowest_levels.end(),
                                    0, [](long s, auto x){ return s + x.second.size(); });
 rel_sp_i.reserve(rank_local_n);
 eigensystems.reserve(rank_local_n);

 CHECK_SIGNALS;

 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Computing eigenvectors ..." << std::endl;

 for(auto const& l : sp_lowest_levels) {
  int n_relevant_ev = std::count_if(l.second.begin(), l.second.end(),
                                    [max_energy](double e){ return e <= max_energy; });
  if(n_relevant_ev == 0) continue;

  auto const& sp = subspaces[l.first];
  if(params.verbosity >= 2)
   std::cout << "[Node " << comm.rank() << "] Computing eigenvectors on subspace "
             << sp.get_index() << " (dimension " << sp.size() << ")" << std::endl;

  rel_sp_i.emplace_back(l.first, comm.rank(), rel_sp_i.size());

  auto ncv_it = params.arpack_ncv.find(l.first);
  int ncv = ncv_it == params.arpack_ncv.end() ? -1 : ncv_it->second;

  if(h_is_real) {
   real_static_op_on_subspace_t op(h_, fops, ist.get_full_hs());
   compute_eigenvectors(sp, op,
                        n_relevant_ev,
                        eigensystems,
                        params.verbosity,
                        params.arpack_min_matrix_size,
                        params.arpack_tolerance, ncv,
                        comm.rank());
  } else {
   static_op_on_subspace_t op(h_, fops, ist.get_full_hs());
   compute_eigenvectors(sp, op,
                        n_relevant_ev,
                        eigensystems,
                        params.verbosity,
                        params.arpack_min_matrix_size,
                        params.arpack_tolerance, ncv,
                        comm.rank());
  }

  CHECK_SIGNALS;

  // Rank-local checks
  assert(rel_sp_i.size() == eigensystems.size());
  assert(first_dim(eigensystems.back().first) == second_dim(eigensystems.back().second));

  // Contribution to the partition function
  if(temperature == 0)
   Z += n_relevant_ev;
  else
   for(auto e : l.second) Z += std::exp(-beta*(e - gs_energy));
 }

 sp_lowest_levels.clear();

 // Complete partition function
 Z = mpi::all_reduce(Z, comm, 0, MPI_SUM);

 // Gather relevant eigenpairs
 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Gathering eigensystems ..." << std::endl;

 // Collect information about relevant subspaces from all ranks
 auto tmp = mpi::all_gather(rel_sp_i, comm, 0);
 std::set<global_index> all_relevant_sp_i(tmp.begin(), tmp.end());

 auto & shs = ist.sub_hilbert_spaces;

 // Ensure that no reallocations occur in ist.sub_hilbert_spaces
 shs.reserve(all_relevant_sp_i.size());

 CHECK_SIGNALS;

 // Finally, fill the init_state object
 eigensystem_t eig;
 for(auto const& spi : all_relevant_sp_i) {
  shs.emplace_back(std::move(subspaces[spi.global]));
  shs.back().set_index(shs.size()-1); // Renumber invariant subspaces in ist

  bool my_rank = comm.rank() == spi.rank;
  auto & energies = my_rank ? eigensystems[spi.local].first : eig.first;
  auto & evec = my_rank ? eigensystems[spi.local].second : eig.second;
  if(comm.size() > 1) {
   mpi_broadcast(energies, comm, spi.rank);
   mpi_broadcast(evec, comm, spi.rank);
  }
  assert(first_dim(energies) > 0 && first_dim(energies) == second_dim(evec));

  for(long i : range(energies.size())) {
   double weight = std::exp(-beta*(energies[i] - gs_energy)) / Z;
   ist.weighted_states.emplace_back(state_on_subspace_t(shs.back()), weight);
   ist.weighted_states.back().state.amplitudes() = evec(range(),i);
  }

  CHECK_SIGNALS;
 }

 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Done. " << ist.weighted_states.size()
            << " weighted pure states in the initial state." << std::endl;

 signal_handler::stop();

 return ist;
}

template
init_state make_equilibrium_init_state(time_expr_operator_t const&,
                                       fundamental_operator_set const&,
                                       double,
                                       eq_solver_parameters_t const&,
                                       std::map<operators::indices_t, int> const&,
                                       mpi::communicator const&);

template
init_state make_equilibrium_init_state(time_interp_operator_t const&,
                                       fundamental_operator_set const&,
                                       double,
                                       eq_solver_parameters_t const&,
                                       std::map<operators::indices_t, int> const&,
                                       mpi::communicator const&);
}
