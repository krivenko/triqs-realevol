#pragma once

#include <complex>
#include <limits>
#include <type_traits>
#include <utility>
#include <tuple>
#include <ostream>
#include <triqs/utility/is_complex.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/operators.hpp>

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
typename std::enable_if<triqs::is_complex<T>::value,bool>::type is_zero(T x)
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
typename std::enable_if<triqs::is_complex<T>::value,T>::type _conj(T x)
{
    return std::conj(x);
}

// Check for being constant
template<class T> bool is_constant(T x) { return true; }

// Callable complex type
template<class T>
class callable_complex : public boost::operators<callable_complex<T>> {
    
    T re, im;
    
    template<typename X>
    struct _T_compatible : public std::is_convertible<typename std::decay<X>::type,T> {};
    
    template<typename X>
    struct _is_complex : public std::integral_constant<bool,
                std::is_same<typename std::decay<X>::type,callable_complex>::value ||
                (triqs::is_complex<X>::value && _T_compatible<typename X::value_type>::value)
                > {};
    
public:
    
    using value_type = T;

    // Constructors
    callable_complex() : re(.0), im(.0) {}
    template<typename R, typename I,
             typename std::enable_if<_T_compatible<R>::value && _T_compatible<I>::value,bool>::type = false>
    callable_complex(R r, I i) : re(r), im(i)
    {}
    template<typename C,
             typename std::enable_if<triqs::is_complex<C>::value && _T_compatible<typename C::value_type>::value,
             bool>::type = false>
    callable_complex(C z) : re(z.real()), im(z.imag())
    {}    
    template<typename R,
             typename std::enable_if<_T_compatible<R>::value,bool>::type = false>
    callable_complex(R r) : re(r), im(.0)
    {}

    // Copy-constructor, move-constructor & assignment
    callable_complex(callable_complex const&) = default;
    callable_complex(callable_complex &&) = default;
    callable_complex & operator=(callable_complex const&) = default;
    
    template<typename R,
             typename std::enable_if<_T_compatible<R>::value,bool>::type = false>
    callable_complex & operator=(R r)
    {
        re = r;
        im = .0;
        return *this;
    }
    
    // Forward the call to the real and imaginary parts
    template<typename... Args>
    auto operator()(Args&&... args) const -> std::complex<typename std::result_of<T(Args...)>::type>
    {
        return {re(std::forward<Args>(args)...),im(std::forward<Args>(args)...)};
    }
        
    // Unary minus
    callable_complex operator-() const
    {
        return callable_complex(-re,-im);
    }
    
    // Addition
    template<typename C,typename std::enable_if<_is_complex<C>::value,bool>::type = false>
    callable_complex operator+=(C z)
    {
        re += z.real();
        im += z.imag();
        return *this;        
    }
    
    template<typename U,typename std::enable_if<_T_compatible<U>::value,bool>::type = false>
    callable_complex operator+=(U x)
    {
        re += x.re;
        return *this;
    }
   
    // Subtraction
    template<typename C,typename std::enable_if<_is_complex<C>::value,bool>::type = false>
    callable_complex operator-=(C z)
    {
        re -= z.real();
        im -= z.imag();
        return *this;        
    }
    
    template<typename U,typename std::enable_if<_T_compatible<U>::value,bool>::type = false>
    callable_complex operator-=(U x)
    {
        re -= x.re;
        return *this;
    }
    
    // Multiplication
    template<typename C,typename std::enable_if<_is_complex<C>::value,bool>::type = false>
    callable_complex operator*=(C z)
    {
        std::tie(re,im) = std::make_tuple(re*z.real() - im*z.imag(),
                                          re*z.imag() + im*z.real());
        return *this;        
    }
    
    template<typename U,typename std::enable_if<_T_compatible<U>::value,bool>::type = false>
    callable_complex operator*=(U x)
    {
        re *= x;
        im *= x;
        return *this;
    }
    
    // Division
    template<typename C,typename std::enable_if<_is_complex<C>::value,bool>::type = false>
    callable_complex operator/=(C z)
    {
        T abs2 = z.real()*z.real() + z.imag()*z.imag();
        std::tie(re,im) = std::make_tuple((re*z.real() + im*z.imag())/abs2,
                                          (im*z.real() - re*z.imag())/abs2);
        return *this;        
    }
    
    template<typename U,typename std::enable_if<_T_compatible<U>::value,bool>::type = false>
    callable_complex operator/=(U x)
    {
        re /= x;
        im /= x;
        return *this;
    }    

    value_type real() const { return re; }
    value_type imag() const { return im; }
    
    bool operator==(callable_complex const& z) const
    {
        return re==z.re && im == z.im;
    }
    
    friend bool is_constant(callable_complex const& z)
    {
        return is_constant(z.re) && is_constant(z.im);
    }

    friend bool is_zero(callable_complex const& z)
    {
        return is_zero(z.re) && is_zero(z.im);
    }

    friend std::ostream& operator<<(std::ostream& os, callable_complex const& z)
    {
        os << "(" << z.re << "," << z.im << ")";
        return os;
    }
    
    friend callable_complex _conj(callable_complex z)
    {
        return callable_complex(z.re,-z.im);
    }

private:

    // Boost serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // serialize base class information
        ar & re;
        ar & im;
    }
};

}