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

#include <triqs/mpi/base.hpp>

namespace realevol {

class mpi_dispatcher {

 std::vector<long> jobs;
 long next_job_index;
 triqs::mpi::communicator comm;

 static constexpr int request_tag = 1024;
 static constexpr int answer_tag = 1025;

 long serial() {
  if(next_job_index >= jobs.size())
   return no_jobs_left;
  else
   return jobs[next_job_index++];
 }

 long master() {
  while(next_job_index < jobs.size()) {
   MPI_Request r1, r2;
   MPI_Irecv(nullptr, 0, MPI_INT, MPI_ANY_SOURCE, request_tag, comm.get(), &r1);
   MPI_Status stat;
   MPI_Wait(&r1,&stat);
   MPI_Isend(&jobs[next_job_index], 1, MPI_LONG, stat.MPI_SOURCE, answer_tag, comm.get(), &r2);
   ++next_job_index;
  }
  comm.barrier();
  return no_jobs_left;
 }

 long worker() {
  MPI_Request r1, r2;
  MPI_Isend(nullptr, 0, MPI_INT, 0, request_tag, comm.get(), &r1);
  long jobid;
  MPI_Irecv(&jobid, 1, MPI_LONG, 0, answer_tag, comm.get(), &r2);
  MPI_Wait(&r2,MPI_STATUS_IGNORE);
  if(jobid == no_jobs_left) comm.barrier();
  return jobid;
 }

public:

 static constexpr long no_jobs_left = -1;

 mpi_dispatcher(triqs::mpi::communicator const& comm, long n_jobs) :
  next_job_index(0), comm(comm) {
  assert(n_jobs > 0);
  if(comm.rank() == 0) {
   int n_workers = comm.size()-1;
   jobs.resize(n_jobs+n_workers);
   std::iota(jobs.begin(), jobs.end()-n_workers, 0);
   std::fill(jobs.end()-n_workers, jobs.end(), long(no_jobs_left));
  }
 }

 mpi_dispatcher(triqs::mpi::communicator const& comm, std::vector<long> const& jobs_) :
  next_job_index(0), comm(comm) {
  assert(n_jobs > 0);
  if(comm.rank() == 0) {
   int n_workers = comm.size()-1;
   jobs.resize(jobs_.size()+n_workers);
   std::copy(jobs_.begin(), jobs_.end(), jobs.begin());
   std::fill(jobs.end()-n_workers, jobs.end(), long(no_jobs_left));
  }
 }

 long operator()() {
  if(comm.size() == 1) return serial();
  else {
   if(comm.rank() == 0) return master();
   else                 return worker();
  }
 }

};

}
