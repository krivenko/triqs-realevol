#pragma once

#include "solver.hpp"
#include <triqs/gfs/meshes/segment.hpp>

namespace realevol {

// Explicit instantiations
// Add more mesh types later if needed
template class solver<triqs::gfs::segment_mesh,true>;

}