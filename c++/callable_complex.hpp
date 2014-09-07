#pragma once

#include <complex>
#include <utility>
#include <ostream>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/base_object.hpp>

namespace realevol {

template<class T>
class callable_complex : public std::complex<T> {
public:

    using base_t = std::complex<T>;
    using base_t::base_t; // Inherit constructors

    // Call operator: forward the call to the real and imaginary parts
    template<typename... Args>
    auto operator()(Args&&... args) const -> std::complex<typename std::result_of<T(Args...)>::type>
    {
        return {base_t::real()(std::forward<Args>(args)...),
                base_t::imag()(std::forward<Args>(args)...)};
    }

private:
    // Boost serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // serialize base class
        ar & boost::serialization::base_object<base_t>(*this);
    }
};

template<class T>
bool is_constant(callable_complex<T> const& z){
    return is_constant(std::real(z)) && is_constant(std::imag(z));
}

}

namespace triqs { namespace utility {

template<typename T>
bool is_zero(realevol::callable_complex<T> const& cte) {
    return triqs::utility::is_zero(cte.real()) &&
           triqs::utility::is_zero(cte.imag());
}

template<typename T>
realevol::callable_complex<T> _conj(realevol::callable_complex<T> const& cte) {
    return {cte.real(),-cte.imag()};
}

}}