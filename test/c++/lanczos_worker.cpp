#include <triqs/utility/first_include.hpp>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <complex>

#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <triqs/arrays/asserts.hpp>

using namespace triqs::arrays;

vector<double> make_zero_state(vector<double> const& st)
{
    vector<double> zero_st(st.size());
    zero_st() = 0;
    return zero_st;
}

template<typename VT>
auto dot_product(VT const& a, VT const& b) -> typename VT::value_type { return dotc(a,b);}

#include "lanczos_worker.hpp"

using realevol::lanczos_worker;

int main() {

    // test lanczos_worker

    { // real scalar
    matrix<double> h(5,5);  // Hamiltonian matrix
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

    for(int n = 0; n < 5; ++n){
        // Check dimensions of Krylov's subspaces
        kw(psi0[n]);
        if(kw.values().size() != n+1) return EXIT_FAILURE;
    }
    }

    { // complex scalar
    std::complex<double> I(0,1);

    matrix<std::complex<double>> h(5,5);  // Hamiltonian matrix
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

    for(int n = 0; n < 5; ++n){
        // Check dimensions of Krylov's subspaces
        lw(psi0[n]);
        if(lw.values().size() != n+1) return EXIT_FAILURE;
    }
    }

    return EXIT_SUCCESS;
}
