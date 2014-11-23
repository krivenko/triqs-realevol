#pragma once

#include <complex>
#include <string>
#include <vector>
#include <map>
#include <boost/mpi/communicator.hpp>

#include <triqs/draft/hilbert_space_tools/state.hpp>
#include <triqs/utility/draft/numeric_ops.hpp>

#include <mesh_container.hpp>
#include <lpt_scheduler.hpp>

namespace realevol {

using triqs::utility::project;
using triqs::utility::hilbert_space;
using triqs::utility::sub_hilbert_space;

template<typename OperatorType, typename Mesh>
class simulation {

    boost::mpi::communicator comm;
    std::vector<int> subspace_to_rank;
    const int missing_rank = -1;

public:

    using results_t = std::map<std::string,mesh_container<double,Mesh>>;;

    simulation(boost::mpi::communicator const& comm,
               std::vector<sub_hilbert_space> subspaces,
               state<hilbert_space,std::complex<double>,false> init_state,
               solve_parameters_t<OperatorType> const& p
              ) :
        comm(comm),
        subspace_to_rank(subspaces.size(),missing_rank)
    {

        std::vector<job_t<int>> scheduler_jobs;
        scheduler_jobs.reserve(subspaces.size());

        for(int spn = 0; spn < subspaces.size(); ++spn) {
            // Project the initial state on this subspace
            auto proj_st = project<state<sub_hilbert_space, std::complex<double>, false>>(init_state,subspaces[spn]);

            // Any non-zero amplitudes in this subspace?
            bool nontrivial = false;
            foreach(proj_st,[&nontrivial](int i, std::complex<double> a){ if(!triqs::utility::is_zero(a)) nontrivial = true; });

            // Complexity is roughly O(subspaces[spn].size())
            if(nontrivial) scheduler_jobs.push_back({spn,static_cast<uint>(subspaces[spn].size())});
        }

        auto schedules = lpt_scheduler(scheduler_jobs,comm.size());

        if(p.verbosity >= 2) {
            std::cout << "Distribution of " << scheduler_jobs.size() << " nontrivial subspaces over "
                      << comm.size() << " nodes:" << std::endl;
            for(int nsch = 0; nsch < schedules.size(); ++nsch){
                if(schedules[nsch].empty()) continue;
                std::cout << "Node " << nsch << ": ";
                for(auto const& job : schedules[nsch]) std::cout << job.index << " ";
                std::cout << std::endl;
            }
        }

        // Fill subspace_to_rank
        for(int nsch = 0; nsch < schedules.size(); ++nsch) {
            for(auto const& job : schedules[nsch]) subspace_to_rank[job.index] = nsch;
        }

        // TODO
        // MPI rank separation starts here
    }

    results_t const& get_results() { return results; }

private:
    results_t results;
};

}