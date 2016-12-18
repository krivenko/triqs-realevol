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
#include <triqs/utility/tuple_tools.hpp>
#include <triqs/utility/exceptions.hpp>

namespace realevol {

template<typename Job, typename... GeneratorArgs>
class mpi_dispatcher {

 // Job generating function; To be called only on master!
 std::function<Job(long, GeneratorArgs...)> job_generator;
 // Number of jobs to be generated
 long n_jobs;
 // Job index to be sent next
 long next_job_index = 0;
 // MPI datatype for Job
 MPI_Datatype job_datatype = triqs::mpi::mpi_datatype<Job>();

 triqs::mpi::communicator comm;

 enum tags : int {request_job_tag = 1024,
                  send_job_tag,
                  generator_arg_tag};

 Job serial(GeneratorArgs... args) {
  if(next_job_index < n_jobs)
   return job_generator(next_job_index++, args...);
  else
   throw(no_jobs_left());
 }

 Job master(GeneratorArgs... args) {
  int n_workers = comm.size()-1;
  MPI_Request recv_r;
  MPI_Recv_init(nullptr, 0, MPI_INT, MPI_ANY_SOURCE, request_job_tag, comm.get(), &recv_r);
  while(next_job_index < n_jobs + n_workers) {
   MPI_Status stat;
   MPI_Start(&recv_r);
   MPI_Wait(&recv_r, &stat);

   auto args = std::tuple<GeneratorArgs...>();
   triqs::tuple::for_each(args, [this,&stat](auto & arg) {
    using triqs::mpi::mpi_datatype;
    auto arg_datatype = mpi_datatype<std14::decay_t<decltype(arg)>>();
    MPI_Recv(&arg, 1, arg_datatype, stat.MPI_SOURCE, generator_arg_tag, comm.get(),
             MPI_STATUS_IGNORE);
   });

   MPI_Request send_r;
   if(next_job_index < n_jobs) { // Still have jobs to send
    auto job = triqs::tuple::apply(job_generator,
                                   std::tuple_cat(std::make_tuple(next_job_index), args));
    MPI_Isend(&job, 1, job_datatype, stat.MPI_SOURCE, send_job_tag, comm.get(), &send_r);
   } else {                      // Sending empty message
    MPI_Isend(nullptr, 0, job_datatype, stat.MPI_SOURCE, send_job_tag, comm.get(), &send_r);
   }
   ++next_job_index;
  }
  MPI_Request_free(&recv_r);

  comm.barrier();
  throw(no_jobs_left());
 }

 Job worker(GeneratorArgs... args) {
  // Notify master node
  MPI_Request send_r;
  MPI_Isend(nullptr, 0, MPI_INT, 0, request_job_tag, comm.get(), &send_r);
  MPI_Wait(&send_r,MPI_STATUS_IGNORE);

  // Send arguments of the job generator
  using triqs::mpi::mpi_datatype;
  triqs::tuple::for_each(std::make_tuple(args...), [this](auto arg) {
   MPI_Send(&arg, 1, mpi_datatype<decltype(arg)>(), 0, generator_arg_tag, comm.get());
  });

  // Receive new job
  Job job;
  MPI_Request recv_r;
  MPI_Status recv_s;
  MPI_Irecv(&job, 1, job_datatype, 0, send_job_tag, comm.get(), &recv_r);
  MPI_Wait(&recv_r,&recv_s);
  int c;
  MPI_Get_count(&recv_s, job_datatype,&c);
  if(c == 1) return job;

  comm.barrier();
  throw(no_jobs_left());
 }

public:

 using job_t = Job;

 /// Construct using a general job generator
 template<typename JobGenerator>
 mpi_dispatcher(triqs::mpi::communicator const& comm, JobGenerator job_generator, long n_jobs) :
  job_generator(job_generator), n_jobs(n_jobs), comm(comm) {
  TRIQS_ASSERT(n_jobs > 0);
 }

 /// Construct and generate jobs as Job::Job(i, args...), i = 0,..., n_jobs-1
 mpi_dispatcher(triqs::mpi::communicator const& comm, long n_jobs) :
  job_generator([](long i, GeneratorArgs... args){ return job_t(i, args...); }),
  n_jobs(n_jobs), comm(comm) {
  TRIQS_ASSERT(n_jobs > 0);
 }

 /// Construct on a list of jobs
 mpi_dispatcher(triqs::mpi::communicator const& comm, std::vector<Job> const& jobs) :
  job_generator([jobs](long i, GeneratorArgs...){ return jobs[i]; }),
  n_jobs(jobs.size()), comm(comm) {
  TRIQS_ASSERT(jobs.size() > 0);
 }

 Job operator()(GeneratorArgs... args) throw (no_jobs_left) {
  if(comm.size() == 1) return serial(args...);
  else {
   if(comm.rank() == 0) return master(args...);
   else                 return worker(args...);
  }
 }

 // Restart job generation
 void reset() { next_job_index = 0; }

 // Throw it from operator()
 // C++17's std::optional would be a better solution.
 class no_jobs_left : std::exception {};
};

}
