/*******************************************************************************
 *
 * realevol - Real time evolution solver based on TRIQS
 *
 * Copyright (C) 2014-2020, I. Krivenko, M. Danilov, P. Kubiczek
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
#pragma once

#include <mpi/mpi.hpp>

namespace realevol {

// -----------------------------------------------------------------
// Little structure to describe a mapping from a communicator-global index
// to a (rank,local index) pair
struct global_index {
 long global; // Global index
 int rank;    // Rank of the owning process
 long local;  // Local index within the owning rank

 global_index() = default;
 global_index(long global, int rank, long local) : global(global), rank(rank), local(local) {}

 friend bool operator<(global_index const& i1, global_index const& i2) {
  return i1.global < i2.global;
 }
};

} // namespace realevol

namespace mpi {

template<> struct mpi_type<realevol::global_index> {
 static MPI_Datatype get() noexcept {
  static bool type_committed = false;
  static MPI_Datatype dt;
  if(!type_committed) {
   int blocklengths[] = {1,1,1};
   MPI_Aint displacements[] = {0,sizeof(long),sizeof(long)+sizeof(int)};
   MPI_Datatype types[] = {MPI_LONG,MPI_INT,MPI_LONG};
   MPI_Type_create_struct(3, blocklengths, displacements, types, &dt);
   MPI_Type_commit(&dt);
   type_committed = true;
  }
  return dt;
 }
};

} // namespace mpi
