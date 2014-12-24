#pragma once

#include <vector>
#include <type_traits>
#include <triqs/arrays.hpp>
#include <triqs/arrays/blas_lapack/stev.hpp>
#include <triqs/utility/is_complex.hpp>
#include <triqs/utility/draft/numeric_ops.hpp>

using namespace triqs::arrays;
using triqs::arrays::blas::tridiag_worker;
using triqs::is_complex;

namespace realevol {

template <typename OperatorType, typename StateType> class lanczos_worker {

    template<typename T>
    struct extract_real_t {
        template<typename U> struct identity { using type = U; };
        template<typename U> struct value_type_of { using type = typename U::value_type; };

        using type = typename std::conditional<
            is_complex<T>::value,value_type_of<T>,identity<T>
        >::type::type;
    };

    using scalar_t = typename StateType::value_type;
    using real_scalar_t = typename extract_real_t<scalar_t>::type;

    OperatorType const& H;

    // Krylov basis states
    // H \approx V * T * V^+
    // Elements of 'basisstates' are columns of V
    std::vector<StateType> basisstates;

    // The tridiagonal matrix T has the following form
    // | alpha[0]    beta[0]     0       ...     |
    // | beta[0]     alpha[1]    beta[1]     ... |
    // | 0            beta[1]    alpha[2]    ... |
    // |             ...                         |
    std::vector<real_scalar_t> alpha; // diagonal matrix elements
    std::vector<real_scalar_t> beta;  // superdiagonal matrix elements

    // Temporaries
    StateType res_vector;

    static constexpr unsigned int reserved_krylov_dim = 20;

    // Convergence threshold for the GS energy
    real_scalar_t gs_energy_convergence;

    // Tridiagonal matrix diagonalizer
    tridiag_worker<false> tdw;

    // For complex numbers: extract the real part and make sure that the imaginary part is negligible
    // For real numbers: returns the argument
    template<bool IsComplex = is_complex<scalar_t>::value>
    real_scalar_t checked_real(scalar_t const& x, typename std::enable_if<IsComplex,void*>::type = 0)
    {
        assert(triqs::utility::is_zero(std::imag(x)));
        return std::real(x);
    }
    template<bool IsComplex = is_complex<scalar_t>::value>
    real_scalar_t checked_real(scalar_t x, typename std::enable_if<!IsComplex,void*>::type = 0)
    {
        return x;
    }

    // Returns the only matrix element of the 1x1 Krylov-projected matrix
    real_scalar_t first_iteration(StateType const& initial_state)
    {
        basisstates.push_back(initial_state);
        res_vector = H(basisstates.back());
        alpha.push_back(checked_real(dot_product(initial_state, res_vector)));
        res_vector -= alpha.back() * initial_state;
        return alpha.back();
    }

    // Calculates the next state in Krylov's basis.
    // Returns false if the previous state was an eigenstate of H
    bool advance()
    {
        real_scalar_t new_beta = std::sqrt(checked_real(dot_product(res_vector, res_vector)));
        // We don't really want to divide by zero
        if(triqs::utility::is_zero(new_beta,gs_energy_convergence)) return false;
        beta.push_back(new_beta);
        basisstates.push_back(res_vector / new_beta);
        res_vector = H(basisstates.back());
        alpha.push_back(checked_real(dot_product(basisstates.back(), res_vector)));
        res_vector -= alpha.back() * basisstates.back();
        res_vector -= beta.back() * basisstates[basisstates.size() - 2];
        return true;
    }

public:

    using state_type = StateType ;

    lanczos_worker(OperatorType const& H, real_scalar_t gs_energy_convergence)
        : H(H), gs_energy_convergence(gs_energy_convergence), tdw(reserved_krylov_dim)
    {
        alpha.reserve(reserved_krylov_dim);
        beta.reserve(reserved_krylov_dim - 1);
        basisstates.reserve(reserved_krylov_dim);
    }
    lanczos_worker(lanczos_worker const&) = default;
    lanczos_worker& operator=(lanczos_worker const&) = delete;

    // initial_state MUST be of norm 1
    void operator()(StateType const& initial_state) {
        reset();

        // First iteration
        real_scalar_t gs_energy = first_iteration(initial_state);
        tdw(alpha, beta); // FIXME: no need to call this... but otherwise tdw.values() can not be called

        while (advance()) {
            tdw(alpha, beta);
            if(triqs::utility::is_zero(tdw.values()[0] - gs_energy,gs_energy_convergence)) break;
            gs_energy = tdw.values()[0];
        }
    }

    // Access eigenvalues and eigenvectors of the Krylov-projected operator
    vector_view<double> values() const { return tdw.values(); }
    matrix_view<double> vectors() const { return tdw.vectors(); }

    void reset() {
        alpha.clear();
        beta.clear();
        basisstates.clear();
    }

    template <typename KrylovCoeffs> StateType krylov_2_fock(KrylovCoeffs const& phi) {
        state_type st = make_zero_state(res_vector);
        for (std::size_t i = 0; i < phi.size(); ++i) st += phi(i) * basisstates[i];
        return st;
    }
};

}

