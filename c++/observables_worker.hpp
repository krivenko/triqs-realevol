#pragma once

#include <string>
#include <map>
#include <algorithm>

#include <triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp>
#include <triqs/draft/hilbert_space_tools/hilbert_space.hpp>
#include <triqs/draft/hilbert_space_tools/imperative_operator.hpp>

#include "solve_parameters.hpp"

namespace realevol {

using namespace triqs::utility;

using block_mapping_t = std::set<std::pair<uint32_t,uint32_t>>;

template<typename ObservablesType, typename SpacePartitionType>
block_mapping_t find_all_connections(ObservablesType const& observables,
                                     fundamental_operator_set const& fops,
                                     SpacePartitionType & SP) {
    block_mapping_t obs_connections;
    for(auto const& obs : observables) {
        auto conn = SP.find_mappings({obs.second,fops});

        std::set_union(std::begin(obs_connections), std::end(obs_connections),
                       std::begin(conn), std::end(conn),
                       std::inserter(obs_connections, std::begin(obs_connections)));
    }
    return obs_connections;
}

void filter_connections(block_mapping_t & connections, std::vector<std::pair<int,int>> const& subspace_disp_table) {
    for(auto it = std::begin(connections); it != std::end(connections); ++it) {
        if(subspace_disp_table[it->first].first == -1 ||
           subspace_disp_table[it->second].first == -1) {
            connections.erase(it);
        }
    }
}

template<bool ComplexOperators>
class observables_worker {

public:

};

}
