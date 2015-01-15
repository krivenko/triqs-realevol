#pragma once

#include <boost/math/special_functions.hpp>

using boost::math::pow;
using boost::math::factorial;;

double three_j_symbol(int j1, int m1,
                      int j2, int m2,
                      int j3, int m3)
{
    // Check if the arguments are physical
    if(m1+m2+m3 != 0 ||
        m1 < -j1 || m1 > j1 ||
        m2 < -j2 || m2 > j2 ||
        m3 < -j3 || m3 > j3 ||
        j3 > j1 + j2 || j3 < std::abs(j1 - j2)
    ) return .0;
        
    double result = ((j1-j2-m3) % 2) ? -1 : 1;
    result *= std::sqrt(
        factorial<double>(j1+j2-j3)*
        factorial<double>(j1-j2+j3)*
        factorial<double>(-j1+j2+j3)/
        factorial<double>(j1+j2+j3+1)
    );
    result *= std::sqrt(
        factorial<double>(j1-m1)*factorial<double>(j1+m1)*
        factorial<double>(j2-m2)*factorial<double>(j2+m2)*
        factorial<double>(j3-m3)*factorial<double>(j3+m3)
    );
    
    int t_min = std::max(std::max(j2-j3-m1,j1-j3+m2),0);
    int t_max = std::min(std::min(j1-m1,j2+m2),j1+j2-j3);
    double t_sum = 0;
    for(int t = t_min; t <= t_max; ++t)
        t_sum += ((t % 2) ? -1 : 1) /
            (
                factorial<double>(t)*
                factorial<double>(j3-j2+m1+t)*
                factorial<double>(j3-j1-m2+t)*
                factorial<double>(j1+j2-j3-t)*
                factorial<double>(j1-m1-t)*
                factorial<double>(j2+m2-t)
            );
       
    result *= t_sum; 
    return result;
}

// Angular matrix elements A_k(m'_1, m'_2, m_2, m_1)
// (l is an angular momentum of a single electron)  
double angular_matrix_element(int l, int k, int mp1, int mp2, int m2, int m1)
{
    double result = 0;
    
    for(int q=-k; q <= k; ++q)
        result += 
            three_j_symbol(l,-mp1,k,-q,l,m2)*
            three_j_symbol(l,-mp2,k,q,l,m1)*
            ((mp1+q+mp2) % 2 ? -1 : 1);
    
    result *= pow<2>(2*l+1)*pow<2>(three_j_symbol(l,0,k,0,l,0));
    
    return result;
}
