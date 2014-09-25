#pragma once

#include <vector>
#include <ostream>
#include <type_traits>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <triqs/utility/dressed_iterator.hpp>

namespace realevol {

template<class T, class Mesh>
class mesh_container : public std::vector<T> {

public:

    using mesh_t = Mesh;
    using base_t = std::vector<T>;
    using value_type = typename base_t::value_type;

    // vector-compatible constructors
    template<class... ValueConstructorArgs>
    explicit mesh_container(const mesh_t& mesh, ValueConstructorArgs && ...vc_args) :
        base_t(mesh.size(),value_type(vc_args...)), mesh(mesh) {}

    mesh_container(const mesh_t& mesh, base_t const& value) : base_t(value), mesh(mesh) {}

    const mesh_t get_mesh() const {
        return mesh;
    }

private:

    using _const_iter = boost::zip_iterator<boost::tuple<
                typename mesh_t::const_iterator,
                typename base_t::const_iterator
    >>;
    using _iter = boost::zip_iterator<boost::tuple<
                typename mesh_t::const_iterator,
                typename base_t::iterator
    >>;

public:

    // Iterator over pairs mesh point-value pairs (const)
    struct const_pair_t {
        typename mesh_t::point_t mesh_point;
        value_type const& value;
        const_pair_t(_const_iter const& it) : mesh_point(boost::get<0>(*it)), value(boost::get<1>(*it)) {}
    };

    using const_iterator = triqs::utility::dressed_iterator<_const_iter,const_pair_t>;
    const_iterator cbegin() const noexcept {
        return const_iterator(boost::make_zip_iterator(
            boost::make_tuple(std::begin(mesh), base_t::cbegin()))
        );
    }
    const_iterator cend() const noexcept {
        return const_iterator(boost::make_zip_iterator(
            boost::make_tuple(std::end(mesh), base_t::cend()))
        );
    }

    // Iterator over pairs mesh point-value pairs
    struct pair_t {
        typename mesh_t::mesh_point_t mesh_point;
        value_type & value;
        pair_t(_iter const& it) : mesh_point(boost::get<0>(*it)), value(boost::get<1>(*it)) {}
    };

    using iterator = triqs::utility::dressed_iterator<_iter,pair_t>;
    iterator begin() noexcept {
        return iterator(boost::make_zip_iterator(
            boost::make_tuple(std::begin(mesh), base_t::begin()))
        );
    }
    iterator end() noexcept {
        return iterator(boost::make_zip_iterator(
            boost::make_tuple(std::end(mesh), base_t::end()))
        );
    }

    // Insert contents of the container into a stream as two columns
    friend std::ostream& operator<<(std::ostream & os, mesh_container const& MC) {
        for(auto e : MC) os << e.mesh_point.value << '\t' << e.value << std::endl;
        return os;
    }

private:
    mesh_t mesh;

    // Methods for boost::serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & mesh;
        ar & boost::serialization::base_object<base_t>(*this);
    }
};

}