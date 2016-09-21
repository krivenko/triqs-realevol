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

int test_sum(mpi_dispatcher & disp, triqs::mpi::communicator const& comm) {
 int jobid;
 int sum = 0;
 std::chrono::milliseconds sleep_duration(0);
 while((jobid = disp()) != mpi_dispatcher::no_jobs_left) {
  std::this_thread::sleep_for(sleep_duration);
  using namespace std::chrono_literals;
  sleep_duration += 1ms;
  sum += jobid;
 }
 sum = mpi_all_reduce(sum, comm, 0, MPI_SUM);

 return sum;
}

TEST(mpi_dispatcher, length) {
 triqs::mpi::communicator comm;
 mpi_dispatcher disp(comm, 50);
 EXPECT_EQ(1225, test_sum(disp, comm));
}

TEST(mpi_dispatcher, lengths_undersub) {
 triqs::mpi::communicator comm;
 mpi_dispatcher disp(comm, 2);
 EXPECT_EQ(1, test_sum(disp, comm));
}

TEST(mpi_dispatcher, vector) {
 triqs::mpi::communicator comm;
 std::vector<long> v(50);
 for(int n = 0; n < 50; ++n) v[n] = n*n;
 mpi_dispatcher disp(comm, v);
 EXPECT_EQ(40425, test_sum(disp, comm));
}

TEST(mpi_dispatcher, vector_undersub) {
 triqs::mpi::communicator comm;
 std::vector<long> v(2);
 for(int n = 0; n < 2; ++n) v[n] = n*n;
 mpi_dispatcher disp(comm, v);
 EXPECT_EQ(1, test_sum(disp, comm));
}
