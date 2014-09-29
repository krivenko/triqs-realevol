#pragma once

#include <vector>
#include <ostream>

#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace realevol {

template<class T, class Mesh>
class mesh_container_cyclic : public std::vector<T> {

public:

    using mesh_t = Mesh;
    using base_t = std::vector<T>;
    using value_type = typename base_t::value_type;

    // vector-compatible constructors
    template<class... ValueConstructorArgs>
    explicit mesh_container_cyclic(const mesh_t& mesh, std::size_t storage_size, ValueConstructorArgs && ...vc_args) :
        base_t(storage_size,value_type(vc_args...)), mesh(mesh) {}

    mesh_container_cyclic(const mesh_t& mesh, base_t const& value) : base_t(value), mesh(mesh) {}

    mesh_t const& get_mesh() const { return mesh; }
    std::size_t size() const { return mesh.size(); }
    std::size_t storage_size() const { return base_t::size(); }

    // Iterator over pairs mesh point-value pairs (const)
    struct const_pair_t {
        typename mesh_t::mesh_point_t mesh_point;
        value_type const& value;
    };

    value_type const& operator[](std::size_t n) const { return base_t::operator[](n % storage_size()); }

    class const_iterator : 
    public boost::iterator_facade<const_iterator,const_pair_t,std::random_access_iterator_tag,const_pair_t>
    {
    public:
        const_iterator(mesh_container_cyclic const& container, std::size_t n) :
        container(&container), n(n) {}

    private:
        std::size_t n;
        mesh_container_cyclic *const container;

        friend class boost::iterator_core_access;
        inline bool equal(const_iterator const& it) const { return n == it.n; }
        inline void increment() { ++n; }
        inline void decrement() { --n; }
        inline void advance(std::size_t i) { n+=i; }
        inline const_pair_t dereference() const { return {container->get_mesh()[n],(*container)[n]}; }
        using diff_t = typename base_t::difference_type;
        inline diff_t distance_to(const_iterator const& it) const { return diff_t(it.n)-diff_t(n); }
    };

    const_iterator cbegin() const noexcept { return {*this,0}; }
    const_iterator cend() const noexcept { return {*this,size()}; }

    // Iterator over pairs mesh point-value pairs
    struct pair_t {
        typename mesh_t::mesh_point_t mesh_point;
        value_type & value;
    };

    value_type & operator[](std::size_t n) { return base_t::operator[](n % storage_size()); }

    class iterator : 
    public boost::iterator_facade<iterator,pair_t,std::random_access_iterator_tag,pair_t>
    {
    public:
        iterator(mesh_container_cyclic & container, std::size_t n) :
        container(&container), n(n) {}

    private:
        std::size_t n;
        mesh_container_cyclic * container;

        friend class boost::iterator_core_access;
        inline bool equal(iterator const& it) const { return n == it.n; }
        inline void increment() { ++n; }
        inline void decrement() { --n; }
        inline void advance(std::size_t i) { n+=i; }
        inline pair_t dereference() const { return {container->get_mesh()[n],(*container)[n]}; }
        using diff_t = typename base_t::difference_type;
        inline diff_t distance_to(iterator const& it) const { return diff_t(it.n)-diff_t(n); }
    };

    iterator begin() noexcept { return {*this,0}; }
    iterator end() noexcept { return {*this,size()}; }

private:
    mesh_t mesh;

    // Methods for boost::serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<base_t>(*this);
        ar & mesh;
    }
};

}