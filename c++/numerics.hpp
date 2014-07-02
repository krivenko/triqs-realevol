#pragma once

#include <complex>
#include <limits>
#include <type_traits>
#include <boost/type_traits/is_complex.hpp>

namespace realevol {

// Generic tests for being zero 
template<class T>
typename std::enable_if<std::is_integral<T>::value,bool>::type is_zero(T x)
{
    return x==0;
}
template<class T>
typename std::enable_if<std::is_floating_point<T>::value,bool>::type is_zero(T x)
{
    return std::fabs(x) < 100*std::numeric_limits<T>::epsilon();
}

template<class T>
typename std::enable_if<boost::is_complex<T>::value,bool>::type is_zero(T x)
{
    return is_zero(std::real(x)) && is_zero(std::imag(x));
}

// Generic complex conjugate
template<class T>
typename std::enable_if<std::is_arithmetic<T>::value,T>::type _conj(T x)
{
    return x;
}

template<class T>
typename std::enable_if<boost::is_complex<T>::value,T>::type _conj(T x)
{
    return std::conj(x);
}
    
}