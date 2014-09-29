#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <complex>

#include <triqs/arrays/vector.hpp>
#include <triqs/arrays/blas_lapack/dot.hpp>
#include <triqs/arrays/asserts.hpp>
#include <triqs/gfs/meshes/segment.hpp>

#include "mesh_container.hpp"
#include "mesh_container_cyclic.hpp"

using namespace realevol;
using namespace triqs::arrays;
using triqs::gfs::segment_mesh;

const std::complex<double> I(0,1.0);

using Vector = vector<std::complex<double>>;
std::complex<double> dot_product(Vector const& a, Vector const& b) { return dotc(a,b);}
Vector make_zero_state(Vector const& V) { return Vector{0,0}; }

#include "lanczos.hpp"

double A1 = 1.0;
double A2 = 2.0;
double B = 0.7;

// Right-hand side of the equation
auto RHS = [](Vector const& sol, double t){
    Vector result(2);
    result[0] = A1*sol[0] + B*sol[1];
    result[1] = B*sol[0] + A2*sol[1];
    return result;
};

// Printing function
auto printer = [](segment_mesh::mesh_point_t mp, Vector const& x){
    std::cout << std::fixed << std::setprecision(5)
    << double(mp) << "     " << x[0].real() << "     " << x[0].imag()
    << "     " << x[1].real() << "     " << x[1].imag()
    << std::endl; };

int main()
{
    // Time mesh on a segment [0;10.0] (201 points)
    segment_mesh mesh(0,10.0,201);

    // Initial conditions
    Vector initial_values(2);
    initial_values[0] = 1.0;
    initial_values[1] = 0.0;

    // Sought solution (normal mesh container)
    mesh_container<Vector,decltype(mesh)> solution1(mesh,2);

    // Solve the equation
    lanczos<decltype(solution1),decltype(RHS)> solver1(RHS);
    solution1[0] = initial_values;
    solver1(begin(solution1),end(solution1),I);

    // Print the solution
    std::cout << "========== Test I ==========" << std::endl;
    for(auto e : solution1) printer(e.mesh_point,e.value);

    // Sought solution (cyclic container)
    mesh_container_cyclic<Vector,decltype(mesh)> solution2(mesh,10,2);

    // Solve the equation & print the solution
    std::cout << "========== Test II =========" << std::endl;
    lanczos<decltype(solution2),decltype(RHS)> solver2(RHS);
    solution2[0] = initial_values;

    auto st_size = solution2.storage_size();
    auto it1 = std::begin(solution2);
    auto it2 = it1 + st_size;
    while(true){
        solver2(it1,it2,I);
        std::for_each(it1,it2-1,[](decltype(solution2)::pair_t e){ printer(e.mesh_point,e.value); });

        if(it2 == std::end(solution2)) break;

        it1 += st_size-1;
        it2 += st_size-1;
        if(it2 > std::end(solution2)) it2 = std::end(solution2);
    }
    printer((it2-1)->mesh_point,(it2-1)->value);

    return EXIT_SUCCESS;
}