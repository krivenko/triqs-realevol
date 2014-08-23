#pragma once

#include <boost/iterator/iterator_facade.hpp>

template<class Mesh> class mesh_iterator;

// Auxiliary MPL-compatible class (metafunction)
template<class Mesh>
struct mesh_iterator_impl {

    using type = boost::iterator_facade<
        mesh_iterator<Mesh>,
        typename Mesh::mesh_point_t const,
        boost::random_access_traversal_tag,
        typename Mesh::mesh_point_t const>
        ;
};

// Iterator over a mesh
template<class Mesh>
class mesh_iterator : public mesh_iterator_impl<Mesh>::type
{
    const Mesh* mesh;
    typename Mesh::node_number_t node;

public:
    using difference_type = typename mesh_iterator_impl<Mesh>::type::difference_type;

    mesh_iterator() = delete;
    explicit mesh_iterator(const Mesh* mesh) : mesh(mesh), node(0) {}
    explicit mesh_iterator(const Mesh* mesh, typename Mesh::node_number_t node) : mesh(mesh), node(node) {}

private:
    friend class boost::iterator_core_access;

    inline const typename Mesh::mesh_point_t dereference() const
    {
        return (*mesh)[node];
    }

    inline void increment() { ++node; }
    inline void decrement() { --node; }

    inline void advance(difference_type n)
    {
        node += n;
    }

    inline difference_type distance_to(mesh_iterator const& other) const
    {
        return difference_type(other.node) - difference_type(node);
    }

    inline bool equal(mesh_iterator const& other) const
    {
        return node == other.node;
    }
};