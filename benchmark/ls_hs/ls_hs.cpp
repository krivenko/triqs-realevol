#include <triqs/utility/first_include.hpp>

#include <array>
#include <string>
#include <map>
#include <utility>
#include <iostream>
#include <cstdlib>

#include <triqs/h5.hpp>

#include <solver.hpp>
#include <time_expr_r.hpp>

#include "moments.hpp"

using namespace realevol;
namespace h5 = triqs::h5;
using triqs::utility::c;
using triqs::utility::c_dag;
using triqs::utility::n;

struct write_obs : public boost::static_visitor<> {

    h5::group & gr;
    std::string name;
    write_obs(h5::group & gr, std::string name) : gr(gr), name(name) {}

    template<typename MeshContainer>
    void operator()(MeshContainer const& sol) const { h5_write(gr,name,sol); }
};

void write_results(var_results_t const& res, h5::group gr) {
    for(auto const& r : res)
        boost::apply_visitor(write_obs(gr,r.first), r.second);
}


int main(int argc, char **argv) {

    std::cout << "This is a simulation of the dynamic low spin/high spin transition" << std::endl;

    // Input parameters
    constexpr int L = 2; // d-orbital
    double mu = 26;

    // Slater parameters (from Haverkort's thesis, page 9)
    double F0 = 5.0;
    double F2 = 10.0;
    double F4 = 6.25;

    // Crystal field splitting
    double _10Dq = 1.0;

    // Bath parameters
    double epsilon = 0;
    double V = 1.0;

    // Time mesh
    triqs::gfs::segment_mesh mesh(0,50.0,1001);

    // Solver method
    ode_solve_method method = method_lanczos;

    constexpr int N_comp = 2*L+1;

    // Indices
    std::array<std::string,2> spin_names = {"up","dn"};
    std::array<std::string,N_comp> cubic_names = {"xy","yz","z^2","xz","x^2-y^2"};
    std::array<std::string,N_comp> bath_sites = {"1","2","3","4","5"};

    // Hamiltonian
    triqs::utility::many_body_operator<time_expr_r> H; 

    // Chemical potential
    for(auto const& s_n : spin_names){
        for(auto const& c_n : cubic_names){
            H += -mu*n<time_expr_r>(s_n,c_n);
        }
    }

    // Full interaction matrix
    // Basis of spherical harmonics Y_{-2}, Y_{-1}, Y_{0}, Y_{1}, Y_{2}
    triqs::arrays::array<double,4> U_matrix_sph(N_comp,N_comp,N_comp,N_comp);

    for(int mp1 = -L; mp1 <= L; ++mp1)
    for(int mp2 = -L; mp2 <= L; ++mp2)
    for(int m1 = -L; m1 <= L; ++m1)
    for(int m2 = -L; m2 <= L; ++m2){
        // m1 <-> m2: 'Field theory' notation
        U_matrix_sph(mp1+L,mp2+L,m1+L,m2+L) = 
            F0 * angular_matrix_element(L,0,mp1,mp2,m2,m1) +
            F2 * angular_matrix_element(L,2,mp1,mp2,m2,m1) +
            F4 * angular_matrix_element(L,4,mp1,mp2,m2,m1);
    }

    // Full interaction matrix
    // Basis of cubic harmonics d_{xy}, d_{yz}, d_{3z^2-r^2}, d_{xz}, d_{x^2-y^2}
    triqs::arrays::array<std::complex<double>,4> U_matrix_cube(N_comp,N_comp,N_comp,N_comp);

    // Basis transformation matrix
    triqs::arrays::matrix<std::complex<double>> W(N_comp,N_comp);
    W() = 0;

    const double inv_sqrt2 = 1.0/std::sqrt(2);
    const std::complex<double> I(0,1.0);

    W(0,0) = I*inv_sqrt2; W(0,4) = -I*inv_sqrt2;
    W(1,1) = I*inv_sqrt2; W(1,3) = I*inv_sqrt2;
    W(2,2) = 1.0;
    W(3,1) = inv_sqrt2; W(3,3) = -inv_sqrt2;
    W(4,0) = inv_sqrt2; W(4,4) = inv_sqrt2;

    auto inv_W = inverse(W);

    // Transform the U-matrix
    for(int ap1 = 0; ap1 < N_comp; ++ap1)
    for(int ap2 = 0; ap2 < N_comp; ++ap2)
    for(int a1 = 0; a1 < N_comp; ++a1)
    for(int a2 = 0; a2 < N_comp; ++a2){
        U_matrix_cube(ap1,ap2,a1,a2) = 0;
        for(int mp1 = 0; mp1 < N_comp; ++mp1)
        for(int mp2 = 0; mp2 < N_comp; ++mp2)
        for(int m1 = 0; m1 < N_comp; ++m1)
        for(int m2 = 0; m2 < N_comp; ++m2) {
            U_matrix_cube(ap1,ap2,a1,a2) +=
                W(ap1,mp1)*W(ap2,mp2)*
                U_matrix_sph(mp1,mp2,m1,m2)*
                inv_W(m1,a1)*inv_W(m2,a2);
        }
        if(std::abs(U_matrix_cube(ap1,ap2,a1,a2)) > 1e-10){

            if(imag(U_matrix_cube(ap1,ap2,a1,a2)) > 1e-10)
                TRIQS_RUNTIME_ERROR << "Cubic harmonics are real, so should be the matrix elements of U.";

            for(auto const & s1 : spin_names)
            for(auto const & s2 : spin_names){

                H += 0.5*real(U_matrix_cube(ap1,ap2,a1,a2))*
                    c_dag<time_expr_r>(s1,cubic_names[ap1])*
                    c_dag<time_expr_r>(s2,cubic_names[ap2])*
                    c<time_expr_r>(s2,cubic_names[a1])*
                    c<time_expr_r>(s1,cubic_names[a2]);
            }
        }
    }

    // Bath sites
    for(auto const& s_n : spin_names){
        for(auto const& bs_n : bath_sites){
            H += -epsilon*n<time_expr_r>(s_n,bs_n);
        }
    }

    // Coupling with the bath
    for(auto const& s_n : spin_names){
        for(auto const& c_n : cubic_names){
            for(auto const& bs_n : bath_sites){
                H += V * c_dag<time_expr_r>(s_n,c_n) * c<time_expr_r>(s_n,bs_n);
                H += V * c_dag<time_expr_r>(s_n,bs_n) * c<time_expr_r>(s_n,c_n);
            }
        }
    }

    // All indices
    std::set<typename solver<false>::indices_t> all_indices;
    for(auto const& s_n : spin_names){
        for(auto const& c_n : cubic_names) all_indices.insert({s_n,c_n});
        for(auto const& bs_n : bath_sites) all_indices.insert({s_n,bs_n});
    }

    // Solver object
    solver<false> S(all_indices);

    // Parameters
    auto params = solve_parameters_t<false>(H,mesh);

    params.verbosity = 2;
    params.stored_psi_values = 20;
    params.method = method;

    // Observables
    std::map<std::string,operator_t<false>> observables;
    params.observables["unity"] = operator_t<false>() + 1.0; // ugly
    params.observables["N_eg"] = n<time_expr_r>("up-z^2") + n<time_expr_r>("dn-z^2")
                               + n<time_expr_r>("up-x^2-y^2") + n<time_expr_r>("dn-x^2-y^2");
    params.observables["N_t2g"] = n<time_expr_r>("up-xy") + n<time_expr_r>("dn-xy")
                                + n<time_expr_r>("up-yz") + n<time_expr_r>("dn-yz")
                                + n<time_expr_r>("up-xz") + n<time_expr_r>("dn-xz");

    h5::file output_file("ls_hs.h5",H5F_ACC_TRUNC);
    h5::group root_gr(output_file);

    // Initial state
    S.psi0({
        {"up","xy"},{"up","yz"},{"up","xz"},    // LS
        {"dn","xy"},{"dn","yz"},{"dn","xz"},    // All t2g states are occupied
        {"up","1"},{"up","2"},{"up","3"},       // 3 of 5 bath levels are occupied
        {"dn","1"},{"dn","2"},{"dn","3"}
    }) = 1.0;

    S.solve(params);

    auto res = S.get_results();
    write_results(res,root_gr);

    return EXIT_SUCCESS;
}
