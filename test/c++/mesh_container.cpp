#include <triqs/test_tools/arrays.hpp>

#include <h5/h5.hpp>

#include <triqs/gfs.hpp>
#include <triqs/utility/numeric_ops.hpp>

#include <realevol/mesh_container.hpp>

using namespace realevol;
using mesh_t = triqs::gfs::gf_mesh<triqs::gfs::retime>;

double sqr(double x) { return x*x; }

TEST(mesh_container, simple) {
 mesh_t m(0,0.5,6);
 mesh_container<double, mesh_t> f1(m);

 for(int i : {0,1,2,3,4,5}) f1[i] = sqr(f1.get_mesh()[i]);

 using triqs::utility::is_zero;
 for(auto e : f1)
  EXPECT_TRUE(is_zero(sqr(double(e.mesh_point)) - e.value, 1e-10));
}

TEST(mesh_container, HDF5) {
 mesh_t m(0,0.5,6);
 mesh_container<double, mesh_t> f1(m);

 for(int i : {0,1,2,3,4,5}) f1[i] = sqr(f1.get_mesh()[i]);

 h5::file hdf_file1("mesh_container.h5", 'w');
 h5_write(hdf_file1, "container", f1);
 hdf_file1.close();

 mesh_container<double, mesh_t> f2(m);
 h5::file hdf_file2("mesh_container.h5", 'r');
 h5_read(hdf_file2, "container", f2);
 hdf_file2.close();

 EXPECT_EQ(f1, f2);
}

MAKE_MAIN;
