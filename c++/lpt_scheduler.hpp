#pragma once

#include <vector>
#include <deque>
#include <utility>
#include <algorithm>

namespace realevol {

template<typename IndexType>
struct job_t {

    IndexType index;
    unsigned int weight;

    friend bool operator==(job_t const& j1, job_t const& j2){ return (j1.index==j2.index); }
    friend bool operator>(job_t const& j1, job_t const& j2){ return j1.weight>j2.weight; }
};

template<typename IndexType>
using schedule_t = std::deque<job_t<IndexType>>;

// Longest Processing Time (LPT) algorithm to solve a multiprocessor scheduling problem
template<typename IndexType>
std::vector<schedule_t<IndexType>> lpt_scheduler(std::vector<job_t<IndexType>> const& jobs, int n_processors) {

    // vector of jobs sorted by their length
    std::vector<job_t<IndexType>> sorted_jobs(jobs);
    std::stable_sort(sorted_jobs.begin(),sorted_jobs.end(),std::greater<job_t<IndexType>>());

    std::vector<schedule_t<IndexType>> schedules(n_processors);
    std::vector<unsigned int> total_lengths(n_processors,0);

    for(auto const& j : sorted_jobs){
        auto least_busy_p = std::min_element(total_lengths.begin(), total_lengths.end());
        *least_busy_p += j.weight;

        schedules[least_busy_p - total_lengths.begin()].push_back(j);
    }

    return schedules;
}

}