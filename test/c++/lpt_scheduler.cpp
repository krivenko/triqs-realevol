#include <triqs/utility/first_include.hpp>

#include <vector>

#include "lpt_scheduler.hpp"

using namespace realevol;

int main() {

    std::vector<job_t<int>> jobs {{0,15},{1,22},{2,14},{3,9},{4,17},{5,31},{6,27},{7,8},{8,15},{9,11},{10,11}};
    auto schedules = lpt_scheduler(jobs,3);

    std::vector<schedule_t<int>> ref_schedules {{{5,31},{8,15},{10,11}},
                                                {{6,27},{0,15},{9,11},{3,9}},
                                                {{1,22},{4,17},{2,14},{7,8}}};

    return schedules != ref_schedules;
}