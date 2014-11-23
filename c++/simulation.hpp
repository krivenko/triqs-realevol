#pragma once

#include <string>
#include <map>

#include <mesh_container.hpp>

namespace realevol {

template<typename OperatorType, typename Mesh>
class simulation {

public:

    using results_t = std::map<std::string,mesh_container<double,Mesh>>;;

    simulation() {}

    results_t const& get_results() { return results; }

private:
    results_t results;
};

}

// Construct workers to solve Schroedinger equation
//     for(auto sp : subspaces) {
//         schroedinger_worker<> worker();
//     }