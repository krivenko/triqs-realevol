/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2016 I. Krivenko
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include <exception>

#include <triqs/mpi/base.hpp>
#include <triqs/utility/exceptions.hpp>

namespace realevol {

template<typename Job> class mpi_dispatcher {

 std::function<Job(long)> job_generator; // To be called only on master!
 long n_jobs;
 long next_job_index;
 triqs::mpi::communicator comm;

 static constexpr int request_job_tag = 1024;
 static constexpr int send_job_tag = 1025;
 static constexpr int no_jobs_left_tag = 1026;

 long serial() {
  if(next_job_index < n_jobs)
   return job_generator(next_job_index++);
  else
   throw(no_jobs_left());
 }

 long master() {
  using triqs::mpi::mpi_datatype;
  int n_workers = comm.size()-1;
  while(next_job_index < n_jobs + n_workers) {
   MPI_Request recv_r, send_r;
   MPI_Irecv(nullptr, 0, MPI_INT, MPI_ANY_SOURCE, request_job_tag, comm.get(), &recv_r);
   MPI_Status stat;
   MPI_Wait(&recv_r,&stat);
   if(next_job_index < n_jobs) { // Sill have jobs to send
    auto job = job_generator(next_job_index);
    MPI_Isend(&job, 1, mpi_datatype<Job>(),
              stat.MPI_SOURCE, send_job_tag, comm.get(), &send_r);
    ++next_job_index;
   } else { // Sending no_jobs_left
    MPI_Isend(nullptr, 0, MPI_INT,
              stat.MPI_SOURCE, no_jobs_left_tag, comm.get(), &send_r);
    ++next_job_index;
   }
  }
  comm.barrier();
  throw(no_jobs_left());
 }

 long worker() {
  using triqs::mpi::mpi_datatype;
  MPI_Request send_r;
  MPI_Request recv_r[2];
  MPI_Isend(nullptr, 0, MPI_INT, 0, request_job_tag, comm.get(), &send_r);
  Job job;
  MPI_Irecv(&job, 1, mpi_datatype<Job>(), 0, send_job_tag, comm.get(), &recv_r[0]);
  MPI_Irecv(nullptr, 0, MPI_INT, 0, no_jobs_left_tag, comm.get(), &recv_r[1]);
  int i;
  MPI_Waitany(2,recv_r,&i,MPI_STATUS_IGNORE);
  if(i==0) {
   MPI_Cancel(&recv_r[1]);
   return job;
  } else {
   MPI_Cancel(&recv_r[0]);
   comm.barrier();
  }
  throw(no_jobs_left());
 }

public:

 using job_t = Job;

 /// Construct using a general job generating function
 template<typename JobGenerator>
 mpi_dispatcher(triqs::mpi::communicator const& comm, JobGenerator job_generator, long n_jobs) :
  job_generator(job_generator), n_jobs(n_jobs), next_job_index(0), comm(comm) {
  TRIQS_ASSERT(n_jobs > 0);
 }

 /// Construct and generate jobs as Job::Job(i), i = 0,..., n_jobs-1
 mpi_dispatcher(triqs::mpi::communicator const& comm, long n_jobs) :
  job_generator([](long i){ return job_t(i); }),
  n_jobs(n_jobs), next_job_index(0), comm(comm) {
  TRIQS_ASSERT(n_jobs > 0);
 }

 /// Construct on a list of jobs
 mpi_dispatcher(triqs::mpi::communicator const& comm, std::vector<Job> const& jobs) :
  job_generator([jobs](long i){ return jobs[i]; }),
  n_jobs(jobs.size()), next_job_index(0), comm(comm) {
  TRIQS_ASSERT(jobs.size() > 0);
 }

 Job operator()() {
  if(comm.size() == 1) return serial();
  else {
   if(comm.rank() == 0) return master();
   else                 return worker();
  }
 }

 class no_jobs_left : std::exception {};
};

}
