#include <cstdlib>

#include <triqs/gfs/meshes/segment.hpp>
#include "mesh_container_cyclic.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

using namespace realevol;
using triqs::gfs::segment_mesh;

template<class T> T sqr(T x) { return x*x; } 

using mesh_t = segment_mesh;
using cont_t = mesh_container_cyclic<int,mesh_t>;

int main()
{
    mesh_t m(0,9.5,20);
    cont_t c(m,6,13);

    std::cout << "size = " << c.size() << std::endl;
    std::cout << "storage_size = " << c.storage_size() << std::endl;

    std::cout << "Container filled with a default value:" << std::endl;
    for(auto p : c) std::cout << double(p.mesh_point) << '\t' << p.value << std::endl;

    c[0] = 0;
    c[1] = 1;
    c[2] = 2;
    c[3] = 3;
    c[4] = 4;
    c[5] = 5;
    c[6] = 6;
    c[7] = 7;

    std::cout << "Updated container I:" << std::endl;
    for(auto const& p : c) std::cout << double(p.mesh_point) << '\t' << p.value << std::endl;

    std::cout << "Updated container II:" << std::endl;
    int n = 0;
    for(auto p : c) p.value = n++;
    for(n=0; n<c.size(); n++) std::cout << c[n] << std::endl;

    return EXIT_SUCCESS;
}