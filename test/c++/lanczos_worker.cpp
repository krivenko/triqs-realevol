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

double dot_product( vector<double> const& a, vector<double> const& b) { return dotc(a,b);}

#include "lanczos_worker.hpp"

using realevol::lanczos_worker;

int main() {

    // test krylov_worker
    matrix<double> h(5,5);  // Hamiltonian matrix
    for(int n : {0,1,2,3,4}) h(n,n) = n;

    auto H = [h](vector<double> const& v){ return h*v; };

    // Initial vectors psi_0
    std::vector<vector<double>> psi0;

    // Eigenstate of H
    psi0.emplace_back(vector<double>{1.,.0,.0,.0,.0});
    // Mixture of 2 eigenstates
    psi0.emplace_back(vector<double>{1.0/sqrt(3.0),sqrt(2.0/3.0),.0,.0,.0});
    // Mixture of 3 eigenstates
    psi0.emplace_back(vector<double>{1.0/sqrt(6.0),1.0/sqrt(3.0),1.0/sqrt(2.0),.0,.0});
    // Mixture of 4 eigenstates
    psi0.emplace_back(vector<double>{1.0/sqrt(10.0),1.0/sqrt(5.0),sqrt(3.0/10.0),sqrt(2.0/5.0),0.0});
    // Mixture of all 5 eigenstates
    psi0.emplace_back(vector<double>{1.0/sqrt(15.0),sqrt(2.0/15.0),1.0/sqrt(5.0),2.0/sqrt(15.0),1.0/sqrt(3.0)});

    lanczos_worker<decltype(H), vector<double>> kw(H,1e-10);

    for(int n = 0; n < 5; ++n){ 
        // Check dimensions of Krylov's subspaces
        kw(psi0[n]);
        if(kw.values().size() != n+1) return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
