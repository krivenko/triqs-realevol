/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2021, I. Krivenko, M. Danilov, P. Kubiczek
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

#include <triqs/test_tools/arrays.hpp>

#include <chrono>
#include <thread>

#include <realevol/mpi_dispatcher.hpp>

using namespace realevol;

int main(int argc, char **argv) {
 mpi::environment env(argc, argv);
 ::testing::InitGoogleTest(&argc, argv);
 return RUN_ALL_TESTS();
}

using namespace std::chrono;
using namespace std::chrono_literals;
using std::this_thread::sleep_for;

template<typename Job = long, typename... GeneratorArgs>
Job test_sum(mpi_dispatcher<Job,GeneratorArgs...> & disp,
             mpi::communicator const& comm,
             GeneratorArgs... args) {
 Job sum = 0;
 milliseconds sleep_duration(0);
 while(true) {
  auto job_or_stop = disp(args...);
  if(!job_or_stop) break;
  auto job = job_or_stop.value();
  sleep_for(sleep_duration);
  sleep_duration += 1ms;
  sum += job;
 }
 sum = mpi::all_reduce(sum, comm, 0, MPI_SUM);

 return sum;
}

TEST(mpi_dispatcher, length) {
 mpi::communicator comm;
 mpi_dispatcher<long> disp(comm, 50);
 EXPECT_EQ(1225, test_sum(disp, comm));
}

TEST(mpi_dispatcher, lengths_undersub) {
 mpi::communicator comm;
 mpi_dispatcher<long> disp(comm, 2);
 EXPECT_EQ(1, test_sum(disp, comm));
}

TEST(mpi_dispatcher, vector) {
 mpi::communicator comm;
 std::vector<long> v(50);
 for(int n = 0; n < 50; ++n) v[n] = n*n;
 mpi_dispatcher<long> disp(comm, v);
 EXPECT_EQ(40425, test_sum(disp, comm));
}

TEST(mpi_dispatcher, vector_undersub) {
 mpi::communicator comm;
 std::vector<long> v(2);
 for(int n = 0; n < 2; ++n) v[n] = n*n;
 mpi_dispatcher<long> disp(comm, v);
 EXPECT_EQ(1, test_sum(disp, comm));
}

TEST(mpi_dispatcher, generator) {
 mpi::communicator comm;
 auto generator = [](long n) { return n*n; };
 mpi_dispatcher<long> disp(comm, generator, 50);
 EXPECT_EQ(40425, test_sum(disp, comm));
}

TEST(mpi_dispatcher, generator_undersub) {
 mpi::communicator comm;
 auto generator = [](long n) { return n*n; };
 mpi_dispatcher<long> disp(comm, generator, 2);
 EXPECT_EQ(1, test_sum(disp, comm));
}

TEST(mpi_dispatcher, generator_args) {
 mpi::communicator comm;
 auto generator = [](long n, int p, int r) { return p*n + 10*r; };
 int N = 50, p = 3, r = 8;
 mpi_dispatcher<long, int, int> disp(comm, generator, N);
 EXPECT_EQ(p*N*(N-1)/2 + 10*r*N, test_sum(disp, comm, p, r));
}

TEST(mpi_dispatcher, generator_args_undersub) {
 mpi::communicator comm;
 auto generator = [](long n, int p, int r) { return p*n + 10*r; };
 int N = 2, p = 3, r = 8;
 mpi_dispatcher<long, int, int> disp(comm, generator, N);
 EXPECT_EQ(p*N*(N-1)/2 + 10*r*N, test_sum(disp, comm, p, r));
}

TEST(mpi_dispatcher, complex) {
 mpi::communicator comm;
 mpi_dispatcher<dcomplex, double> disp(comm, 50);
 EXPECT_CLOSE(dcomplex(1225,150), test_sum(disp, comm, 3.0));
}

TEST(mpi_dispatcher, complex_undersub) {
 mpi::communicator comm;
 mpi_dispatcher<dcomplex, double> disp(comm, 2);
 EXPECT_CLOSE(dcomplex(1,6), test_sum(disp, comm, 3.0));
}
