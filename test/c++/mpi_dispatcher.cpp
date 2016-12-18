#include <triqs/test_tools/arrays.hpp>

#include <chrono>
#include <thread>

#include "mpi_dispatcher.hpp"

using namespace realevol;

int main(int argc, char **argv) {
 triqs::mpi::environment env(argc, argv);
 ::testing::InitGoogleTest(&argc, argv);
 return RUN_ALL_TESTS();
}

using namespace std::chrono;
using namespace std::chrono_literals;
using std::this_thread::sleep_for;

template<typename Job = long, typename... GeneratorArgs>
Job test_sum(mpi_dispatcher<Job,GeneratorArgs...> & disp,
             triqs::mpi::communicator const& comm,
             GeneratorArgs... args) {
 Job sum = 0;
 milliseconds sleep_duration(0);
 try {
  while(true) {
   auto job = disp(args...);
   sleep_for(sleep_duration);
   sleep_duration += 1ms;
   sum += job;
  }
 } catch(typename mpi_dispatcher<Job,GeneratorArgs...>::no_jobs_left &) {}
 sum = mpi_all_reduce(sum, comm, 0, MPI_SUM);

 return sum;
}

TEST(mpi_dispatcher, length) {
 triqs::mpi::communicator comm;
 mpi_dispatcher<long> disp(comm, 50);
 EXPECT_EQ(1225, test_sum(disp, comm));
}

TEST(mpi_dispatcher, lengths_undersub) {
 triqs::mpi::communicator comm;
 mpi_dispatcher<long> disp(comm, 2);
 EXPECT_EQ(1, test_sum(disp, comm));
}

TEST(mpi_dispatcher, vector) {
 triqs::mpi::communicator comm;
 std::vector<long> v(50);
 for(int n = 0; n < 50; ++n) v[n] = n*n;
 mpi_dispatcher<long> disp(comm, v);
 EXPECT_EQ(40425, test_sum(disp, comm));
}

TEST(mpi_dispatcher, vector_undersub) {
 triqs::mpi::communicator comm;
 std::vector<long> v(2);
 for(int n = 0; n < 2; ++n) v[n] = n*n;
 mpi_dispatcher<long> disp(comm, v);
 EXPECT_EQ(1, test_sum(disp, comm));
}

TEST(mpi_dispatcher, generator) {
 triqs::mpi::communicator comm;
 auto generator = [](long n) { return n*n; };
 mpi_dispatcher<long> disp(comm, generator, 50);
 EXPECT_EQ(40425, test_sum(disp, comm));
}

TEST(mpi_dispatcher, generator_undersub) {
 triqs::mpi::communicator comm;
 auto generator = [](long n) { return n*n; };
 mpi_dispatcher<long> disp(comm, generator, 2);
 EXPECT_EQ(1, test_sum(disp, comm));
}

TEST(mpi_dispatcher, generator_args) {
 triqs::mpi::communicator comm;
 auto generator = [](long n, int p, int r) { return p*n + 10*r; };
 int N = 50, p = 3, r = 8;
 mpi_dispatcher<long, int, int> disp(comm, generator, N);
 EXPECT_EQ(p*N*(N-1)/2 + 10*r*N, test_sum(disp, comm, p, r));
}

TEST(mpi_dispatcher, generator_args_undersub) {
 triqs::mpi::communicator comm;
 auto generator = [](long n, int p, int r) { return p*n + 10*r; };
 int N = 2, p = 3, r = 8;
 mpi_dispatcher<long, int, int> disp(comm, generator, N);
 EXPECT_EQ(p*N*(N-1)/2 + 10*r*N, test_sum(disp, comm, p, r));
}

TEST(mpi_dispatcher, complex) {
 triqs::mpi::communicator comm;
 mpi_dispatcher<dcomplex, double> disp(comm, 50);
 EXPECT_CLOSE(dcomplex(1225,150), test_sum(disp, comm, 3.0));
}

TEST(mpi_dispatcher, complex_undersub) {
 triqs::mpi::communicator comm;
 mpi_dispatcher<dcomplex, double> disp(comm, 2);
 EXPECT_CLOSE(dcomplex(1,6), test_sum(disp, comm, 3.0));
}
