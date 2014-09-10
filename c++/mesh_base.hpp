#pragma once

#include <type_traits>
#include <boost/serialization/access.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace realevol {

// Abstract mesh
template<
    class NodeIndex,    // Node index type (unsigned)
    class Value,        // Node value type (floating point)
    class Mesh          // CRTP
>
struct mesh_base {

    using node_index_t = NodeIndex;
    using value_t = Value;

    static_assert(std::is_unsigned<node_index_t>::value,"NodeIndex is not an unsigned integral type.");
    static_assert(std::is_floating_point<value_t>::value,"Value is not a floating-point type.");

    mesh_base(node_index_t nodes) : nodes(nodes) {}

    // Number of nodes of the mesh
    inline node_index_t size() const { return nodes; }

    inline bool is_on_mesh(node_index_t node) const {
        return node < nodes;
    }

    struct point_t {
        typename Mesh::node_index_t index;
        typename Mesh::value_t value;
    };

private:

    template<typename Derived>
    struct iterator_base_t {
        using type = boost::iterator_facade<Derived,point_t const,boost::random_access_traversal_tag,point_t const>;
    };

public:

    // Iterator over a mesh
    class const_iterator : public iterator_base_t<const_iterator>::type
    {
        Mesh const& mesh;
        typename Mesh::node_index_t node;

    public:
        using difference_type = typename iterator_base_t<const_iterator>::type::difference_type;
        using point_t = mesh_base::point_t;

        const_iterator() = delete;
        const_iterator(const_iterator const&) = default;
        const_iterator(const_iterator &&) = default;
        explicit const_iterator(Mesh const& mesh) : mesh(mesh), node(0) {}
        explicit const_iterator(Mesh const& mesh, typename Mesh::node_index_t node) : mesh(mesh), node(node) {}

    private:
        friend class boost::iterator_core_access;
        inline const point_t dereference() const {
            return {node,mesh[node]};
        }
        inline void increment() { ++node; }
        inline void decrement() { --node; }
        inline void advance(difference_type n){
            node += n;
        }
        inline difference_type distance_to(const_iterator const& other) const {
            return other.node - node;
        }
        inline bool equal(const_iterator const& other) const {
            return node == other.node;
        }
    };

    const_iterator begin(void) const noexcept {
        return const_iterator(static_cast<Mesh const&>(*this));
    }
    const_iterator cbegin(void) const noexcept {
        return const_iterator(static_cast<Mesh const&>(*this));
    }
    const_iterator end(void) const noexcept {
        return const_iterator(static_cast<Mesh const&>(*this),nodes);
    }
    const_iterator cend(void) const noexcept {
        return const_iterator(static_cast<Mesh const&>(*this),nodes);
    }

protected:
    node_index_t nodes;

    // Methods for boost::serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & nodes;
    }
};

}