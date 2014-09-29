#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <functional>

#include <triqs/arrays/vector.hpp>

#include <triqs/gfs/meshes/segment.hpp>
#include "mesh_container.hpp"
#include "mesh_container_cyclic.hpp"
#include "runge_kutta.hpp"

using namespace realevol;
using namespace triqs::arrays;
using triqs::gfs::segment_mesh;

// Solve a motion equation of a classical harmonic oscillator with damping

double w0 = 1.0; // frequency
double g = 0.2; // damping

using Vector = vector<double>;
enum EqType {Coordinate = 0, Velocity = 1};

// Right-hand side of the equation
auto RHS = [](Vector const& sol, double t){
    Vector result(2);
    result[Coordinate] = sol[Velocity];
    result[Velocity] = -w0*w0*sol[Coordinate] - 2*g*w0*sol[Velocity];
    return result;
};

// Printing function
auto printer = [](segment_mesh::mesh_point_t mp, Vector const& x){
    std::cout << std::fixed << std::setprecision(5) << double(mp) << "    " << x[Coordinate] << std::endl; }; 

int main()
{
    // Time mesh on a segment [0;20.0] (201 points)
    segment_mesh mesh(0,20.0,201);

    // Initial conditions
    Vector initial_values(2);
    initial_values[Coordinate] = 1.0;
    initial_values[Velocity] = 0.0;

    // Sought solution (normal mesh container)
    mesh_container<Vector,decltype(mesh)> solution1(mesh,2);

    // Solve the equation
    runge_kutta<decltype(solution1),decltype(RHS)> solver1(RHS);
    solution1[0] = initial_values;
    solver1(begin(solution1),end(solution1));

    // Print the solution
    std::cout << "========== Test I ==========" << std::endl;
    for(auto e : solution1) printer(e.mesh_point,e.value);

    // Sought solution (cyclic container)
    mesh_container_cyclic<Vector,decltype(mesh)> solution2(mesh,10,2);

    // Solve the equation & print the solution
    std::cout << "========== Test II =========" << std::endl;
    runge_kutta<decltype(solution2),decltype(RHS)> solver2(RHS);
    solution2[0] = initial_values;

    auto st_size = solution2.storage_size();
    auto it1 = std::begin(solution2);
    auto it2 = it1 + st_size;
    while(true){
        solver2(it1,it2);
        std::for_each(it1,it2-1,[](decltype(solution2)::pair_t e){ printer(e.mesh_point,e.value); });

        if(it2 == std::end(solution2)) break;

        it1 += st_size-1;
        it2 += st_size-1;
        if(it2 > std::end(solution2)) it2 = std::end(solution2);
    }
    printer((it2-1)->mesh_point,(it2-1)->value);

    return EXIT_SUCCESS;
}