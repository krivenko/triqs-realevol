#include <cstdlib>
#include <cmath>
#include <sstream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/concept_check.hpp>

#include "uniform_mesh.hpp"

using namespace realevol;

template<typename Mesh>
bool check_uniform_mesh(Mesh const& m){

    double ref_v[] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

    auto it = m.begin();
    for(; it != m.end(); ++it){
        unsigned int index = std::distance(m.begin(), it);
        if(std::fabs(it->value - ref_v[index]) >= 1e-10) return false;
    }

    return true;
}

int main(void)
{
    uniform_mesh<unsigned int,double> m1(0,1,11);

    if(m1[5] != 0.5) return EXIT_FAILURE;

    if(!check_uniform_mesh(m1)) return EXIT_FAILURE;

    std::stringstream ss;
    boost::archive::text_oarchive oar(ss);
    oar << m1;
    boost::archive::text_iarchive iar(ss);
    uniform_mesh<unsigned int,double> m2(2,30,21);
    iar >> m2;

    if(!check_uniform_mesh(m2)) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}