/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2023, I. Krivenko, M. Danilov, P. Kubiczek
 *
 * realevol is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * realevol is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * realevol. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <benchmark/benchmark.h>

#include <triqs/gfs.hpp>

#include <realevol/time_expr.hpp>

using namespace realevol;

triqs::gfs::gf_mesh<triqs::gfs::retime> mesh(1.0, 5.0, 41);

static void time_expr_real_const(benchmark::State &state) {
  int n_terms = state.range(0);

  time_expr expr(.0);
  for(int n = 0; n < n_terms; ++n)
    expr += time_expr(2.0*n);

  for(auto _ : state) {
    for(auto t : mesh)
      benchmark::DoNotOptimize(expr(double(t)));
  }
}

static void time_expr_complex_const(benchmark::State &state) {
  int n_terms = state.range(0);

  time_expr expr(.0);
  for(int n = 0; n < n_terms; ++n)
    expr += time_expr(2.0*n - 2.0i*n);

  for(auto _ : state) {
    for(auto t : mesh)
      benchmark::DoNotOptimize(expr(double(t)));
  }
}

static void time_expr_real(benchmark::State &state) {
  int n_terms = state.range(0);

  time_expr expr(.0);
  for(int n = 0; n < n_terms; ++n)
    expr += time_expr(1.0/(n+1)) * time_expr("cos(pi*" + std::to_string(n) + "*t)");

  for(auto _ : state) {
    for(auto t : mesh)
      benchmark::DoNotOptimize(expr(double(t)));
  }
}

static void time_expr_complex(benchmark::State &state) {
  int n_terms = state.range(0);

  time_expr expr(.0);
  for(int n = 0; n < n_terms; ++n) {
    auto n_str = std::to_string(n);
    expr += time_expr(1.0/(n+1)) *
            time_expr("cos(pi*" + n_str + "*t)", "sin(pi*" + n_str + "*t)");
  }

  for(auto _ : state) {
    for(auto t : mesh)
      benchmark::DoNotOptimize(expr(double(t)));
  }
}

int max_n_terms = 1024;

// Register benchmarks
BENCHMARK(time_expr_real_const)->Range(0, max_n_terms);
BENCHMARK(time_expr_complex_const)->Range(0, max_n_terms);
BENCHMARK(time_expr_real)->Range(0, max_n_terms);
BENCHMARK(time_expr_complex)->Range(0, max_n_terms);

BENCHMARK_MAIN();
