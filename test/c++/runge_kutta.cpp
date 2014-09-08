#include <cstdlib>
#include <iostream>
#include <iomanip>

#include <triqs/arrays/vector.hpp>

#include "uniform_mesh.hpp"
#include "mesh_container.hpp"
#include "fifo_mesh_container.hpp"
#include "runge_kutta.hpp"

using namespace realevol;
using namespace triqs::arrays;

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
auto printer = [](double t, Vector const& x){ std::cout << std::fixed << std::setprecision(5)
                                              << t << "    " << x[Coordinate] << std::endl; }; 

int main()
{
    // Time mesh on a segment [0;20.0] (201 points)
    uniform_mesh<> mesh(0,20.0,201);

    // Initial conditions
    Vector initial_values(2);
    initial_values[Coordinate] = 1.0;
    initial_values[Velocity] = 0.0;

    // Sought solution (normal mesh container)
    mesh_container<Vector,decltype(mesh)> solution1(mesh,2);

    // Solve the equation
    runge_kutta<decltype(solution1),decltype(RHS)> solver1(RHS);
    solver1(solution1,initial_values);

    // Print the solution
    std::cout << "========== Test I ==========" << std::endl;
    for(auto it=solution1.arg_value_begin(); it!=solution1.arg_value_end(); ++it)
        printer(it->get<0>(),it->get<1>());

    // Sought solution (fifo mesh container)
    fifo_mesh_container<Vector,decltype(mesh)> solution2(mesh,{printer,10,5},2);

    // Solve the equation & print the solution
    runge_kutta<decltype(solution2),decltype(RHS)> solver2(RHS);
    std::cout << "========== Test II =========" << std::endl;
    solver2(solution2,initial_values);

    return EXIT_SUCCESS;
}