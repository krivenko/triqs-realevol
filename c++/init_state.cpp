/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 I. Krivenko
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
#include <triqs/hilbert_space/space_partition.hpp>
#include <triqs/hilbert_space/state_view.hpp>
#include <triqs/arrays/vector.hpp>
#include <triqs/arrays/linalg/eigenelements.hpp>
#include <triqs/mpi/vector.hpp>

#include "arpack/arpack_worker.hpp"
#include "init_state.hpp"
#include "mpi_dispatcher.hpp"

namespace realevol {

std::ostream & operator<<(std::ostream & os, init_state const& st) {
 os << "Fundamental operator set:" << std::endl;
 os << " Fermions: ";
 for(auto it = st.fops.begin(Fermion); it != st.fops.end(Fermion); ++it)
  os << "(" << it->index << ") ";
 os << std::endl << " Bosons: ";
 for(auto it = st.fops.begin(Boson); it != st.fops.end(Boson); ++it)
  os << "(" << it->index << ") ";

 os << std::endl << "Dimension of fermionic space: " << st.full_hs.size(Fermion) << std::endl;
 os << "Dimension of bosonic space: " << st.full_hs.size(Boson) << std::endl;

 os << st.sub_hilbert_spaces.size() << " relevant invariant subspaces:" << std::endl;
 for(auto const& sp : st.sub_hilbert_spaces) {
  os << " Subspace " << sp.get_index() << " [" << sp.size() << "]: ";
  for(auto f : sp.get_all_fock_states()) os << f << " ";
  os << std::endl;
 }
 os << std::endl;

 for(int i = 0; i < st.weighted_states.size(); ++i) {
  auto const& s = st.weighted_states[i];
  os << "State " << i << " within subspace " << s.state.get_hilbert().get_index()
     << ", weight = " << s.weight << ":" << s.state << std::endl;
 }

 return os;
}

// -----------------------------------------------------------------

void h5_write(h5::group fg, std::string const &name, init_state const& st) {
 auto gr = fg.create_group(name);
 gr.write_triqs_hdf5_data_scheme(st);

 h5_write_attribute(gr, "fops", st.fops);
 h5_write(gr, "full_hs", st.full_hs);
 h5_write(gr, "sub_hilbert_spaces", st.sub_hilbert_spaces);

 auto states_gr = gr.create_group("weighted_states");
 for(int i = 0; i < st.weighted_states.size(); ++i) {
  auto wst_gr = states_gr.create_group(std::to_string(i));
  auto const& wst = st.weighted_states[i];

  h5_write(wst_gr, "weight", wst.weight);
  h5_write(wst_gr, "sp_index", wst.state.get_hilbert().get_index());
  h5_write(wst_gr, "amplitudes", wst.state.amplitudes());
 }
}

// -----------------------------------------------------------------

void h5_read(h5::group fg, std::string const &name, init_state & st) {
 auto gr = fg.open_group(name);

 h5_read_attribute(gr, "fops", st.fops);
 h5_read(gr, "full_hs", st.full_hs);
 h5_read(gr, "sub_hilbert_spaces", st.sub_hilbert_spaces);

 auto states_gr = gr.open_group("weighted_states");
 int n_wst = states_gr.get_all_subgroup_names().size();
 for(int i = 0; i < n_wst; ++i) {
  auto wst_gr = states_gr.open_group(std::to_string(i));

  double weight;
  h5_read(wst_gr, "weight", weight);
  unsigned long sp_index;
  h5_read(wst_gr, "sp_index", sp_index);
  st.weighted_states.emplace_back(state_on_subspace_t(st.sub_hilbert_spaces[sp_index]), weight);
  h5_read(wst_gr, "amplitudes", st.weighted_states.back().state.amplitudes());
 }
}

// -----------------------------------------------------------------

static_operator_t make_static_op(operator_t const& op, std::string const& error_message) {
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

init_state make_pure_init_state(operator_t const& generator,
                                fundamental_operator_set const& fops,
                                std::map<operators::indices_t, int> const& bits_per_boson) {

 // Static version of generator
 auto gen = make_static_op(generator, "Generating operator must be time-independent!");

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
 foreach(psi, [&sp](fock_state_t f, dcomplex){ sp.add_fock_state(f); });

 // Add a weighted state
 state_on_subspace_t st_on_sp(sp);
 foreach(psi, [&st_on_sp,&sp](fock_state_t f, dcomplex a){ st_on_sp(sp.get_state_index(f)) = a; });
 ist.weighted_states.emplace_back(std::move(st_on_sp), 1.0);

 return ist;
}

// -----------------------------------------------------------------

using namespace triqs::arrays;
using namespace triqs::arrays::arpack;
using sp_levels_t = std::pair<int,vector<double>>;
using eigensystem_t = std::pair<vector<double>,matrix<dcomplex>>;
using real_state_on_subspace_t = state<sub_hilbert_space,double,false>;
using real_static_op_on_subspace_t = imperative_operator<sub_hilbert_space,double,false>;

/////////////////////////////////////////////////////////
/// Solver procedures for make_equilibrium_init_state ///
/////////////////////////////////////////////////////////

auto make_arpack_worker_params(int nev, bool eigenvectors, std::true_type) {
 using params_t = arpack_worker<Complex>::params_t;
 return params_t(nev, params_t::SmallestReal, eigenvectors ? params_t::Ritz : params_t::None);
}

auto make_arpack_worker_params(int nev, bool eigenvectors, std::false_type) {
 using params_t = arpack_worker<Symmetric>::params_t;
 return params_t(nev, params_t::Smallest, eigenvectors);
}

// -----------------------------------------------------------------
// Little structure to describe a mapping from a communicator-global index
// to a (rank,local index) pair
struct global_index {
 int global; // Global index
 int rank;   // Rank of the owning process
 int local;  // Local index within the owning rank

 global_index() = default;
 global_index(int global, int rank, int local) : global(global), rank(rank), local(local) {}

 friend bool operator<(global_index const& i1, global_index const& i2) {
  return i1.global < i2.global;
 }
};

} // namespace realevol

namespace triqs { namespace mpi {
template<> inline MPI_Datatype mpi_datatype<realevol::global_index>() {
 static bool type_committed = false;
 static MPI_Datatype dt;
 if(!type_committed) {
  int blocklengths[] = {1,1,1};
  MPI_Aint displacements[] = {0,sizeof(int),2*sizeof(int)};
  MPI_Datatype types[] = {MPI_INT,MPI_INT,MPI_INT};
  MPI_Type_create_struct(3, blocklengths, displacements, types, &dt);
  MPI_Type_commit(&dt);
  type_committed = true;
 }
 return dt;
}
}}

// -----------------------------------------------------------------

namespace realevol {

template<typename T>
double find_lowest_levels_on_subspace(sub_hilbert_space const& sp,
                                      imperative_operator<sub_hilbert_space,T,false> const& h,
                                      std::vector<sp_levels_t> & lowest_levels,
                                      double exc_energy_threshold,
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
  lowest_levels.emplace_back(sp_levels_t{sp_n, {energy}});
  return energy;

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
  double e0 = eig(0);
  for(int i : range(1,N)) {
   if(eig(i) - e0 > exc_energy_threshold) {
    lowest_levels.emplace_back(sp_levels_t{sp_n, eig(range(i))});
    return e0;
   }
  }
  lowest_levels.emplace_back(sp_levels_t{sp_n, eig});
  return e0;

 }

 if(verbosity >= 2)
  std::cout << "[Node " << rank << "] Using ARPACK on subspace "
            << sp_n << std::endl;

 // Set up ARPACK worker
 arpack_worker<is_complex_t::value ? Complex : Symmetric> arpw(sp.size());

 using state_view_t = state_view<sub_hilbert_space,T>;
 std::array<state_view_t,3> state_views = {
  state_view_t{arpw.workspace_vector(0), sp},
  state_view_t{arpw.workspace_vector(1), sp},
  state_view_t{arpw.workspace_vector(2), sp}
 };

 auto params = make_arpack_worker_params(1, false, is_complex_t());
 params.tolerance = arpack_tol;
 params.ncv = arpack_ncv;

 int first_ev_index;

 // Partially diagonalize op within sp iteratively increasing
 // the number of eigenvalues to be found
 while(true) {
  ++params.n_eigenvalues;
  if(verbosity >= 2)
  std::cout << "[Node " << rank << "] Calling ARPACK to find the lowest "
            << params.n_eigenvalues << " eigenvalues" << std::endl;

  auto apply_h = [&h,&state_views](vector_view<T>, int from_n,
                                   vector_view<T>, int to_n) {
   h.apply(state_views[from_n], state_views[to_n]);
  };

  // Run ARPACK solver
  arpw(apply_h, params);

  auto const& eig = arpw.eigenvalues();

   // ARPACK cannot calculate the last eigenvalue. Stopping ...
  if(params.n_eigenvalues == sp.size()-(is_complex_t::value ? 2 : 1)) {
   std::cout << "[Node " << rank << "] WARNING: ARPACK has found "
             << params.n_eigenvalues
             << " relevant eigenvalues on subspace "
             << sp.get_index() << ". "
             << (is_complex_t::value ? "Two eigenvalues" : "One eigenvalue")
             << " can still be missing" << std::endl;
   first_ev_index = 0;
   break;
  }

  // The most recently found eigenvalue was too high
  if(real(eig(0) - eig(params.n_eigenvalues-1)) > exc_energy_threshold) {
   first_ev_index = 1;
   break;
  }
 }

 auto const& eig = arpw.eigenvalues();

 vector<double> ev(eig.size() - first_ev_index);
 for(int i = 0; i < eig.size() - first_ev_index; ++i)
  ev(i) = real(eig(eig.size()-1-i));

 lowest_levels.emplace_back(sp_levels_t{sp_n, std::move(ev)});
 auto const& e = lowest_levels.back().second;
 return lowest_levels.back().second(0);
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

 // Set up ARPACK worker
 arpack_worker<is_complex_t::value ? Complex : Symmetric> arpw(sp.size());

 using state_view_t = state_view<sub_hilbert_space,T>;
 std::array<state_view_t,3> state_views = {
  state_view_t{arpw.workspace_vector(0), sp},
  state_view_t{arpw.workspace_vector(1), sp},
  state_view_t{arpw.workspace_vector(2), sp}
 };

 auto params = make_arpack_worker_params(n_vectors_to_compute, true, is_complex_t());
 params.tolerance = arpack_tol;
 params.ncv = arpack_ncv;

 auto apply_h = [&h,&state_views](vector_view<T>, int from_n,
                                  vector_view<T>, int to_n) {
  h.apply(state_views[from_n], state_views[to_n]);
 };

 // Run ARPACK solver
 arpw(apply_h, params);

 eigensystems.emplace_back(real(arpw.eigenvalues()),
                           arpw.eigenvectors()(range(),range(n_vectors_to_compute)));
}

// -----------------------------------------------------------------
init_state make_equilibrium_init_state(operator_t const& h,
                                       fundamental_operator_set const& fops,
                                       double temperature,
                                       eq_solver_parameters const& params,
                                       std::map<operators::indices_t, int> const& bits_per_boson,
                                       triqs::mpi::communicator const& comm) {

 // Static version the equilibrium Hamiltonian
 auto h_ = make_static_op(h, "Initial Hamiltonian must be time-independent!");

 // Check that h is Hermitian
 if(!(h_ - dagger(h_)).is_zero())
  TRIQS_RUNTIME_ERROR << "Supplied Hamitonian is not Hermitian!";

 if(params.arpack_min_matrix_size < 4)
  TRIQS_RUNTIME_ERROR << "arpack_min_matrix_size must be >= 4!";

 double beta, exc_energy_threshold;
 if(temperature == 0)
  exc_energy_threshold = 1e-14;
 else {
  beta = 1 / temperature;
  exc_energy_threshold = -temperature * std::log(params.min_rel_weight);
 }

 // Initial state object
 init_state ist(fops, bits_per_boson);

 // Partition the Hilbert space
 space_partition<state_on_space_t, static_op_on_space_t>
 partition(state_on_space_t(ist.get_full_hs()), static_op_on_space_t(h_, fops, ist.get_full_hs()));

 // Fill subspaces
 std::vector<sub_hilbert_space> subspaces;
 subspaces.reserve(partition.n_subspaces());
 for (long n = 0; n < partition.n_subspaces(); ++n) subspaces.emplace_back(n);
 foreach(partition, [&](fock_state_t s, int spn) { subspaces[spn].add_fock_state(s); });

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

 mpi_dispatcher disp(comm, subspaces.size());
 int jobid;

 // Ground state anergy on this MPI rank
 double gs_energy = std::numeric_limits<double>::infinity();

 // Rank-local list of processed subspaces and their lowest levels
 std::vector<sp_levels_t> sp_lowest_levels;

 while((jobid = disp()) != mpi_dispatcher::no_jobs_left) {
  auto const& sp = subspaces[jobid];
  if(params.verbosity >= 2)
   std::cout << "[Node " << comm.rank() << "] Searching the lowest eigenvalues on subspace "
             << sp.get_index() << std::endl;

  double new_gs_energy;

  auto ncv_it = params.arpack_ncv.find(sp.get_index());
  int ncv = ncv_it == params.arpack_ncv.end() ? -1 : ncv_it->second;

  if(h_is_real) {
   real_static_op_on_subspace_t op(h_, fops, ist.get_full_hs());
   new_gs_energy = find_lowest_levels_on_subspace(sp, op, sp_lowest_levels,
                                                  exc_energy_threshold,
                                                  params.verbosity,
                                                  params.arpack_min_matrix_size,
                                                  params.arpack_tolerance, ncv,
                                                  comm.rank());
  } else {
   static_op_on_subspace_t op(h_, fops, ist.get_full_hs());
   new_gs_energy = find_lowest_levels_on_subspace(sp, op, sp_lowest_levels,
                                                  exc_energy_threshold,
                                                  params.verbosity,
                                                  params.arpack_min_matrix_size,
                                                  params.arpack_tolerance, ncv,
                                                  comm.rank());
  }

  if(params.verbosity >= 2)
   std::cout << "[Node " << comm.rank() << "] The lowest levels on subspace "
             << sp.get_index() << ": "
             << sp_lowest_levels.back().second << std::endl;

  // New rank-local energy minimum
  gs_energy = std::min(gs_energy, new_gs_energy);
 }

 // Find the global energy minimum
 gs_energy = mpi_all_reduce(gs_energy, comm, 0, MPI_MIN);

 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Ground state energy: " << gs_energy << std::endl;

 // Diagonalization phase 2
 // Now we know the exact ground state energy, and can properly select
 // the relevant subspaces. On those we compute the eigenvectors.
 double max_energy = gs_energy + exc_energy_threshold;

 // Rank-local quantities
 double Z = 0;                               // Partition function
 std::vector<global_index> rel_sp_i;         // List of relevant subspace indices
 std::vector<eigensystem_t> eigensystems;    // List of eigensystems

 int rank_local_n = std::accumulate(sp_lowest_levels.begin(), sp_lowest_levels.end(),
                                    0, [](long s, auto x){ return s + x.second.size(); });
 rel_sp_i.reserve(rank_local_n);
 eigensystems.reserve(rank_local_n);

 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Computing eigenvectors ..." << std::endl;

 for(auto const& l : sp_lowest_levels) {
  int n_relevant_ev = std::count_if(l.second.begin(), l.second.end(),
                                    [max_energy](double e){ return e <= max_energy; });
  if(n_relevant_ev == 0) continue;

  rel_sp_i.emplace_back(l.first, comm.rank(), rel_sp_i.size());

  auto ncv_it = params.arpack_ncv.find(l.first);
  int ncv = ncv_it == params.arpack_ncv.end() ? -1 : ncv_it->second;

  if(h_is_real) {
   real_static_op_on_subspace_t op(h_, fops, ist.get_full_hs());
   compute_eigenvectors(subspaces[l.first], op,
                        n_relevant_ev,
                        eigensystems,
                        params.verbosity,
                        params.arpack_min_matrix_size,
                        params.arpack_tolerance, ncv,
                        comm.rank());
  } else {
   static_op_on_subspace_t op(h_, fops, ist.get_full_hs());
   compute_eigenvectors(subspaces[l.first], op,
                        n_relevant_ev,
                        eigensystems,
                        params.verbosity,
                        params.arpack_min_matrix_size,
                        params.arpack_tolerance, ncv,
                        comm.rank());
  }

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
 Z = mpi_all_reduce(Z, comm, 0, MPI_SUM);

 // Gather relevant eigenpairs
 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Gathering eigensystems ..." << std::endl;

 // Collect information about relevant subspaces from all ranks
 auto tmp = mpi_gather(rel_sp_i, comm, 0, true, std::true_type());
 std::set<global_index> all_relevant_sp_i(tmp.begin(), tmp.end());

 auto & shs = ist.sub_hilbert_spaces;

 // Ensure that no reallocations occur in ist.sub_hilbert_spaces
 shs.reserve(all_relevant_sp_i.size());

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
 }

 if(params.verbosity >= 1 && comm.rank() == 0)
  std::cout << "Done. " << ist.weighted_states.size()
            << " weighted pure states in the initial state." << std::endl;

 return ist;
}

}
