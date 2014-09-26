#include <cstdlib>
#include <cmath>
#include <sstream>

#include <triqs/gfs/meshes/segment.hpp>
#include "mesh_container.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

using namespace realevol;
using triqs::gfs::segment_mesh;

double sqr(double x) { return x*x; } 

bool check_mesh_container(mesh_container<double,segment_mesh> & f)
{
    for(auto e : f)
        if(std::fabs(sqr(double(e.mesh_point)) - e.value) > 1e-10)
            return false;
    return true;
}

int main()
{
    segment_mesh m(0,0.5,6);
    mesh_container<double,segment_mesh> f1(m);

    f1[0] = sqr(f1.get_mesh()[0]);
    f1[1] = sqr(f1.get_mesh()[1]);
    f1[2] = sqr(f1.get_mesh()[2]);
    f1[3] = sqr(f1.get_mesh()[3]);
    f1[4] = sqr(f1.get_mesh()[4]);
    f1[5] = sqr(f1.get_mesh()[5]);

    if(!check_mesh_container(f1)) return EXIT_FAILURE;

    // Test serialization
    std::stringstream ss;
    boost::archive::text_oarchive oar(ss);
    oar << f1;
    boost::archive::text_iarchive iar(ss);
    mesh_container<double,segment_mesh> f2(m);
    iar >> f2;

    if(!check_mesh_container(f2)) return EXIT_FAILURE;

    // Test HDF5
    triqs::h5::file hdf_file1("mesh_container.h5",H5F_ACC_TRUNC);
    h5_write(hdf_file1, "container", f1);
    hdf_file1.close();

    mesh_container<double,segment_mesh> f3(m);
    triqs::h5::file hdf_file2("mesh_container.h5",H5F_ACC_RDONLY);
    h5_read(hdf_file2, "container", f3);
    hdf_file2.close();

    if(!check_mesh_container(f3)) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}