#pragma once
                   
#include <vector>
#include <ostream>
#include <type_traits>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "mesh_base.hpp"
#include "uniform_mesh.hpp"

namespace realevol {

template<class T, class Mesh = uniform_mesh<> >
struct mesh_container : public std::vector<T> {
    
public:
    
    using mesh_t = Mesh;
    using value_type = typename std::vector<T>::value_type;
    
    static_assert(std::is_base_of<mesh_base<typename mesh_t::node_number_t, typename mesh_t::mesh_point_t,mesh_t>,mesh_t>::value,
                  "Mesh is not derived from mesh_base");
    
    // vector-compatible constructors
    explicit mesh_container(const mesh_t& mesh,
                            const value_type& value = value_type()) :
        std::vector<value_type>(mesh.get_nodes(),value), mesh(mesh) {}

    template <class InputIterator>
    mesh_container(const mesh_t& mesh, InputIterator it) :
        std::vector<value_type>(it, it + mesh.get_nodes()), mesh(mesh) {}
    
    mesh_container(const mesh_t& mesh, const std::vector<value_type>& value) :
        std::vector<value_type>(value), mesh(mesh) {}

    const mesh_t get_mesh(void) const
    {
        return mesh;
    }

    // Iterator over pairs mesh point-value pairs (const)
    using arg_value_const_iterator =
        boost::zip_iterator<boost::tuple<
            typename mesh_t::const_iterator,
            typename std::vector<T>::const_iterator
        >>;
    
    arg_value_const_iterator arg_value_begin(void) const noexcept
    {
        return boost::make_zip_iterator(
            boost::make_tuple(std::begin(mesh), std::vector<value_type>::begin())
        );
    }

    arg_value_const_iterator arg_value_end(void) const noexcept
    {
        return boost::make_zip_iterator(
            boost::make_tuple(std::end(mesh), std::vector<value_type>::end())
        );
    }
    
    // Iterator over pairs mesh point-value pairs
    using arg_value_iterator =
        boost::zip_iterator<boost::tuple<
            typename mesh_t::const_iterator,
            typename std::vector<T>::iterator
        >>;
    
    arg_value_iterator arg_value_begin(void) noexcept
    {
        return boost::make_zip_iterator(
            boost::make_tuple(std::begin(mesh), std::vector<value_type>::begin())
        );
    }

    arg_value_iterator arg_value_end(void) noexcept
    {
        return boost::make_zip_iterator(
            boost::make_tuple(std::end(mesh), std::vector<value_type>::end())
        );
    }
    
    // Insert contents of the container into a stream as two columns
    friend std::ostream& operator<<(std::ostream & os, mesh_container const& MC)
    {
        for(arg_value_const_iterator it = MC.arg_value_begin();
            it != MC.arg_value_end(); ++it)
            os << std::get<0>(*it) << '\t' << std::get<1>(*it) << std::endl;

        return os;
    }
    
private:
    Mesh mesh;
    
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