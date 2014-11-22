#pragma once

#include <boost/variant.hpp>
#include <boost/mpl/transform.hpp>

#include <triqs/gfs/meshes/segment.hpp>

#include "./mesh_container.hpp"
#include "./mesh_container_cyclic.hpp"

namespace realevol {

using any_mesh_t = boost::variant<triqs::gfs::segment_mesh>;

template<typename T>
struct mesh_to_container_impl {
    template<typename Mesh> struct apply {
        using type = mesh_container<T,Mesh>;
    };
};

template<typename T>
using any_mesh_container_t = typename boost::make_variant_over<
                                typename boost::mpl::transform<
                                    typename any_mesh_t::types,
                                    mesh_to_container_impl<T>
                                >::type
                             >::type;

template<typename T>
struct mesh_to_container_cyclic_impl {
    template<typename Mesh> struct apply {
        using type = mesh_container<T,Mesh>;
    };
};

template<typename T>
using any_mesh_container_cyclic_t = typename boost::make_variant_over<
                                        typename boost::mpl::transform<
                                        typename any_mesh_t::types,
                                        mesh_to_container_cyclic_impl<T>
                                    >::type
                                >::type;

}