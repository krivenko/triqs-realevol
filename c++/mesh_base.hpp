#pragma once

#include <type_traits>
#include <boost/serialization/access.hpp>

#include "mesh_iterator.hpp"

namespace realevol {

// Abstract mesh
template<
    class NodeNumber,   // Node index type (unsigned)
    class MeshPoint,    // Node value type (floating point)
    class Mesh          // CRTP
> 
struct mesh_base {
    
    using node_number_t = NodeNumber;
    using mesh_point_t = MeshPoint;
    
    static_assert(std::is_unsigned<node_number_t>::value,"NodeNumber is not an unsigned integral type.");
    static_assert(std::is_floating_point<mesh_point_t>::value,"MeshPoint is not a floating-point type.");
    
    mesh_base(node_number_t nodes) : nodes(nodes) {}
    
    // Number of nodes of the mesh
    inline node_number_t get_nodes() const { return nodes; }
    
    inline bool is_on_mesh(node_number_t node) const
    {
        return node < nodes;
    }

    // Mesh iterator
    using const_iterator = mesh_iterator<Mesh>;
    
    const_iterator begin(void) const noexcept
    {
        return const_iterator(static_cast<Mesh const*>(this));
    }
    
    const_iterator cbegin(void) const noexcept
    {
        return const_iterator(static_cast<Mesh const*>(this));
    }
    
    const_iterator end(void) const noexcept
    {  
        return const_iterator(static_cast<Mesh const*>(this),nodes);
    }

    const_iterator cend(void) const noexcept
    {  
        return const_iterator(static_cast<Mesh const*>(this),nodes);
    }
    
protected:
    node_number_t nodes;
    
    // Methods for boost::serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & nodes;
    }
    
};

}