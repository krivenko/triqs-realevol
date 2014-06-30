#pragma once

#include"mesh_base.hpp"
#include"mesh_iterator.hpp"

#include<boost/serialization/access.hpp>
#include<boost/serialization/base_object.hpp>

namespace realevol {

// Uniform mesh
template<
    class NodeNumber = std::size_t,    // Node index type (unsigned)
    class MeshPoint = double           // Node value type (floating point)
>
struct uniform_mesh :
    public mesh_base<NodeNumber,MeshPoint,uniform_mesh<NodeNumber,MeshPoint>> {
    
    using node_number_t = NodeNumber;
    using mesh_point_t = MeshPoint;
    using base_t =  mesh_base<node_number_t,mesh_point_t,uniform_mesh<node_number_t,mesh_point_t>>;
        
    uniform_mesh(mesh_point_t start, mesh_point_t end, node_number_t nodes) :
        base_t(nodes),
        start(start),
        step((end-start)/(nodes-1)) {}
                
    inline mesh_point_t operator[](node_number_t node) const
    {   
        return start + step*node;
    }
    
    inline mesh_point_t get_step() const { return step; }
    
private:
    mesh_point_t start;
    mesh_point_t step;
    
    // Methods for boost::serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<base_t>(*this);
        ar & start;
        ar & step;
    }
};

}