#pragma once

#include <vector>
#include <ostream>
#include <type_traits>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <triqs/utility/dressed_iterator.hpp>

#include "mesh_base.hpp"

namespace realevol {

template<class T, class Mesh>
class mesh_container : public std::vector<T> {

public:

    using mesh_t = Mesh;
    using value_type = typename std::vector<T>::value_type;

    static_assert(std::is_base_of<mesh_base<typename mesh_t::node_index_t, typename mesh_t::mesh_point_t,mesh_t>,mesh_t>::value,
                  "Mesh is not derived from mesh_base");

    // vector-compatible constructors
    template<class... ValueConstructorArgs>
    explicit mesh_container(const mesh_t& mesh, ValueConstructorArgs && ...vc_args) :
        std::vector<value_type>(mesh.size(),value_type(vc_args...)), mesh(mesh) {}

    mesh_container(const mesh_t& mesh, const std::vector<value_type>& value) :
        std::vector<value_type>(value), mesh(mesh) {}

    const mesh_t get_mesh(void) const
    {
        return mesh;
    }

private:

    using _const_iter = boost::zip_iterator<boost::tuple<
                typename mesh_t::const_iterator,
                typename std::vector<T>::const_iterator
    >>;
    using _iter = boost::zip_iterator<boost::tuple<
                typename mesh_t::const_iterator,
                typename std::vector<T>::iterator
    >>;

public:

    // Iterator over pairs mesh point-value pairs (const)
    struct const_deref_result_t {
        typename mesh_t::const_iterator::deref_result_t mesh_point;
        value_type const& value;
        const_deref_result_t(_const_iter const& it) : mesh_point(boost::get<0>(*it)), value(boost::get<1>(*it)) {}
    };

    using const_iterator = triqs::utility::dressed_iterator<_const_iter,const_deref_result_t>;
    const_iterator cbegin() const noexcept {
        return const_iterator(boost::make_zip_iterator(
            boost::make_tuple(std::begin(mesh), std::vector<value_type>::cbegin()))
        );
    }
    const_iterator cend() const noexcept {
        return const_iterator(boost::make_zip_iterator(
            boost::make_tuple(std::end(mesh), std::vector<value_type>::cend()))
        );
    }

    // Iterator over pairs mesh point-value pairs
    struct deref_result_t {
        typename mesh_t::const_iterator::deref_result_t mesh_point;
        value_type & value;
        deref_result_t(_iter const& it) : mesh_point(boost::get<0>(*it)), value(boost::get<1>(*it)) {}
    };

    using iterator = triqs::utility::dressed_iterator<_iter,deref_result_t>;
    iterator begin() noexcept {
        return iterator(boost::make_zip_iterator(
            boost::make_tuple(std::begin(mesh), std::vector<value_type>::begin()))
        );
    }
    iterator end() noexcept {
        return iterator(boost::make_zip_iterator(
            boost::make_tuple(std::end(mesh), std::vector<value_type>::end()))
        );
    }

    // Insert contents of the container into a stream as two columns
    friend std::ostream& operator<<(std::ostream & os, mesh_container const& MC)
    {
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
        ar & boost::serialization::base_object<std::vector<T> >(*this);
    }
};

}