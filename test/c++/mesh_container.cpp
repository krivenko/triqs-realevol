#include <cstdlib>
#include <cmath>
#include <sstream>

#include "uniform_mesh.hpp"
#include "mesh_container.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

using namespace realevol;

template<class T> T sqr(T x) { return x*x; } 

bool check_mesh_container(mesh_container<double,uniform_mesh<>> & f)
{    
    mesh_container<double,uniform_mesh<>>::arg_value_iterator i = f.arg_value_begin();
    for(; i != f.arg_value_end(); i++){
        if(std::fabs(sqr(i->get<0>()) - i->get<1>()) > 1e-10)
            return false;
    }
    
    return true;
}

int main(void)
{
    uniform_mesh<> m(0,0.5,6);
    mesh_container<double,uniform_mesh<>> f1(m);
    
    f1[0] = sqr(f1.get_mesh()[0]);
    f1[1] = sqr(f1.get_mesh()[1]);
    f1[2] = sqr(f1.get_mesh()[2]);
    f1[3] = sqr(f1.get_mesh()[3]);
    f1[4] = sqr(f1.get_mesh()[4]);
    f1[5] = sqr(f1.get_mesh()[5]);

    if(!check_mesh_container(f1)) return EXIT_FAILURE;
    
    std::stringstream ss;
    boost::archive::text_oarchive oar(ss);
    oar << f1;
    boost::archive::text_iarchive iar(ss);
    mesh_container<double,uniform_mesh<>> f2(m);
    iar >> f2;
    
    if(!check_mesh_container(f2)) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}