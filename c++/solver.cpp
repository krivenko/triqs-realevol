#include <triqs/utility/first_include.hpp>

#include "solver.hpp"

namespace realevol {

struct index_visitor  {
 std::vector<std::string> indices;
 void operator()(int i) { indices.push_back(std::to_string(i)); }
 void operator()(std::string s) { indices.push_back(s); }
};

solver::solver(std::map<std::string,indices_type> const& gf_struct, std::pair<double,double> time_window, int n_t) :
  gf_struct(gf_struct) {

 std::vector<std::string> block_names;
 std::vector<g_2t_t> g_ret_blocks;
 std::vector<g_2t_t> g_adv_blocks;

 for (auto const& bl : gf_struct) {
  block_names.push_back(bl.first);
  int n = bl.second.size();

  index_visitor iv;
  for (auto & ind: bl.second) { apply_visitor(iv, ind); }
  std::vector<std::vector<std::string>> indices{{iv.indices,iv.indices}};

  auto t_mesh = gf_mesh<retime>{time_window.first, time_window.second, n_t};

  g_ret_blocks.push_back(g_2t_t{{t_mesh, t_mesh}, {n, n}, indices});
  g_adv_blocks.push_back(g_2t_t{{t_mesh, t_mesh}, {n, n}, indices});
 }

 g_ret = make_block_gf(block_names, g_ret_blocks);
 g_adv = make_block_gf(block_names, g_adv_blocks);
}

void solver::solve(solve_parameters_t const& p) {
}

}
