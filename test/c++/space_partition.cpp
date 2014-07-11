#include <triqs/utility/first_include.hpp>

#include <set>
#include <vector>
#include <algorithm>
#include <sstream>

#include <triqs/operators/many_body_operator.hpp>

#include "space_partition.hpp"
#include "triqs/draft/hilbert_space_tools/fundamental_operator_set.hpp"
#include "triqs/draft/hilbert_space_tools/hilbert_space.hpp"
#include "triqs/draft/hilbert_space_tools/imperative_operator.hpp"
#include "triqs/draft/hilbert_space_tools/state.hpp"

using namespace realevol;
using namespace triqs::utility;

int main() {

 // 3 bands Kanamori

 // Parameters
 double mu = 0.7;
 double U = 3.0;
 double J = 0.3;

 // basis of operators to use
 fundamental_operator_set fops;
 for (int o = 0; o < 3; ++o) {
  fops.insert("up", o);
  fops.insert("dn", o);
 }

 // Hamiltonian
 many_body_operator<double> H;
 for (int o = 0; o < 3; ++o) {
  H += -mu * (n("up", o) + n("dn", o));
 }
 for (int o = 0; o < 3; ++o) {
  H += U * n("up", o) * n("dn", o);
 }
 for (int o1 = 0; o1 < 3; ++o1)
  for (int o2 = 0; o2 < 3; ++o2) {
   if (o1 == o2) continue;
   H += (U - 2 * J) * n("up", o1) * n("dn", o2);
  }
 for (int o1 = 0; o1 < 3; ++o1)
  for (int o2 = 0; o2 < 3; ++o2) {
   if (o2 >= o1) continue;
   H += (U - 3 * J) * n("up", o1) * n("up", o2);
   H += (U - 3 * J) * n("dn", o1) * n("dn", o2);
  }

 for (int o1 = 0; o1 < 3; ++o1)
  for (int o2 = 0; o2 < 3; ++o2) {
   if (o1 == o2) continue;
   H += -J * c_dag("up", o1) * c_dag("dn", o1) * c("up", o2) * c("dn", o2);
   H += -J * c_dag("up", o1) * c_dag("dn", o2) * c("up", o2) * c("dn", o1);
  }

 // Hilbert space
 hilbert_space hs(fops);

 // Imperative operator for H
 imperative_operator<hilbert_space, double, false> Hop(H, fops);

 // Sample state
 state<hilbert_space, double, true> st(hs);

 // Space partition
 using space_partition_t = space_partition<state<hilbert_space, double, true>, imperative_operator<hilbert_space, double, false>>;
 space_partition_t SP(st, Hop);

 /////////////////////////////
 // Part I: Check subspaces //
 /////////////////////////////

 // Calculated classification of states
 // sets are used to neglect order of subspaces and of states within a subspace
 std::vector<std::set<int>> v_cl(SP.n_subspaces());
 foreach(SP, [&v_cl](int st, int spn) { v_cl[spn].insert(st); });
 std::set<std::set<int>> cl{v_cl.cbegin(), v_cl.cend()};

 int u0 = 1 << fops[{"up", 0}];
 int u1 = 1 << fops[{"up", 1}];
 int u2 = 1 << fops[{"up", 2}];
 int d0 = 1 << fops[{"dn", 0}];
 int d1 = 1 << fops[{"dn", 1}];
 int d2 = 1 << fops[{"dn", 2}];

 // Expected classification of states
 std::set<std::set<int>> ref_cl{
     // N=0
     {0},
     // N=1
     {d0},
     {d1},
     {d2},
     {u0},
     {u1},
     {u2},
     // N=2, same spin
     {d0 + d1},
     {d0 + d2},
     {d1 + d2},
     {u0 + u1},
     {u0 + u2},
     {u1 + u2},
     // N=2, pair hopping
     {d0 + u0, d1 + u1, d2 + u2},
     // N=2, spin flip
     {d0 + u1, d1 + u0},
     {d0 + u2, d2 + u0},
     {d1 + u2, d2 + u1},
     // N=3
     {d0 + d1 + d2},
     {u0 + u1 + u2},
     {d0 + d1 + u0, d1 + d2 + u2},
     {d0 + d2 + u0, d1 + d2 + u1},
     {d0 + d1 + u1, d0 + d2 + u2},
     {d0 + u0 + u1, d2 + u1 + u2},
     {d1 + u0 + u1, d2 + u0 + u2},
     {d0 + u0 + u2, d1 + u1 + u2},
     {d1 + d2 + u0, d0 + d2 + u1, d0 + d1 + u2},
     {d2 + u0 + u1, d0 + u1 + u2, d1 + u0 + u2},
     // N=4, 2 holes with the same spin
     {d2 + u0 + u1 + u2},
     {d1 + u0 + u1 + u2},
     {d0 + u0 + u1 + u2},
     {d0 + d1 + d2 + u2},
     {d0 + d1 + d2 + u1},
     {d0 + d1 + d2 + u0},
     // N=4, pair hopping
     {d1 + d2 + u1 + u2, d0 + d2 + u0 + u2, d0 + d1 + u0 + u1},
     // N=4, spin flip
     {d1 + d2 + u0 + u2, d0 + d2 + u1 + u2},
     {d1 + d2 + u0 + u1, d0 + d1 + u1 + u2},
     {d0 + d2 + u0 + u1, d0 + d1 + u0 + u2},
     // N=5
     {d1 + d2 + u0 + u1 + u2},
     {d0 + d2 + u0 + u1 + u2},
     {d0 + d1 + u0 + u1 + u2},
     {d0 + d1 + d2 + u1 + u2},
     {d0 + d1 + d2 + u0 + u2},
     {d0 + d1 + d2 + u0 + u1},
     // N=6
     {d0 + d1 + d2 + u0 + u1 + u2}};

 if (cl != ref_cl) return -1;

 return 0;
}
