/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2022, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <triqs/utility/first_include.hpp>
#include <triqs/test_tools/arrays.hpp>

#include <cstdlib>
#include <vector>
#include <iostream>
#include <complex>

#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <triqs/arrays/asserts.hpp>

using namespace triqs::arrays;

vector<double> make_zero_state(vector<double> const& st) {
 vector<double> zero_st(st.size());
 zero_st() = 0;
 return zero_st;
}

template<typename VT>
auto dot_product(VT const& a, VT const& b) -> typename VT::value_type { return dotc(a,b);}

#include <realevol/lanczos_worker.hpp>

using realevol::lanczos_worker;

TEST(lanczos_worker, real) {
 matrix<double> h(5,5);  // Hamiltonian matrix
 h() = 0;
 for(int n : {0,1,2,3,4}) h(n,n) = n;

 auto H = [&h](vector<double> const& v){ return h*v; };

 // Initial vectors psi_0
 std::vector<vector<double>> psi0;

 // Eigenstate of H
 psi0.push_back({1.,.0,.0,.0,.0});
 // Mixture of 2 eigenstates
 psi0.push_back({1.0/sqrt(3.0),sqrt(2.0/3.0),.0,.0,.0});
 // Mixture of 3 eigenstates
 psi0.push_back({1.0/sqrt(6.0),1.0/sqrt(3.0),1.0/sqrt(2.0),.0,.0});
 // Mixture of 4 eigenstates
 psi0.push_back({1.0/sqrt(10.0),1.0/sqrt(5.0),sqrt(3.0/10.0),sqrt(2.0/5.0),0.0});
 // Mixture of all 5 eigenstates
 psi0.push_back({1.0/sqrt(15.0),sqrt(2.0/15.0),1.0/sqrt(5.0),2.0/sqrt(15.0),1.0/sqrt(3.0)});

 lanczos_worker<decltype(H), vector<double>> kw(H,1e-10);

 for(int n = 0; n < 5; ++n) {
  // Check dimensions of Krylov's subspaces
  kw(psi0[n]);
  EXPECT_EQ(n+1, kw.values().size());
 }
}

TEST(lanczos_worker, complex) {
 std::complex<double> I(0,1);

 matrix<std::complex<double>> h(5,5);  // Hamiltonian matrix
 h() = 0;
 for(int n : {0,1,2,3,4}) h(n,n) = n;

 auto H = [&h](vector<std::complex<double>> const& v){ return h*v; };

 // Initial vectors psi_0
 std::vector<vector<std::complex<double>>> psi0;

 // Eigenstate of H
 psi0.push_back({I,{.0},{.0},{.0},{.0}});
 // Mixture of 2 eigenstates
 psi0.push_back({I/sqrt(3.0),{sqrt(2.0/3.0)},{.0},{.0},{.0}});
 // Mixture of 3 eigenstates
 psi0.push_back({I/sqrt(6.0),{1.0/sqrt(3.0)},{1.0/sqrt(2.0)},{.0},{.0}});
 // Mixture of 4 eigenstates
 psi0.push_back({I/sqrt(10.0),{1.0/sqrt(5.0)},{sqrt(3.0/10.0)},{sqrt(2.0/5.0)},{.0}});
 // Mixture of all 5 eigenstates
 psi0.push_back({I/sqrt(15.0),{sqrt(2.0/15.0)},{1.0/sqrt(5.0)},{2.0/sqrt(15.0)},{1.0/sqrt(3.0)}});

 lanczos_worker<decltype(H), vector<std::complex<double>>> lw(H,1e-10);

 for(int n = 0; n < 5; ++n) {
  // Check dimensions of Krylov's subspaces
  lw(psi0[n]);
  EXPECT_EQ(n+1, lw.values().size());
 }
}

MAKE_MAIN;
