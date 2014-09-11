#pragma once

#include <deque>
#include <functional>
#include <type_traits>
#include <stdexcept>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include "mesh_base.hpp"

namespace realevol {

template<class T, class Mesh>
struct fifo_mesh_container {

public:

    using mesh_t = Mesh;
    using value_type = T;
    using size_type = typename std::deque<value_type>::size_type;

    static_assert(std::is_base_of<mesh_base<typename mesh_t::node_index_t,typename mesh_t::value_t,mesh_t>,mesh_t>::value,
                  "Mesh is not derived from mesh_base");

    struct overflow_t {
        std::function<void(fifo_mesh_container&)> overflow_handler;
        size_type max_size;

        template<class OverflowHandlerType>
        overflow_t(OverflowHandlerType overflow_handler, size_type max_size = 1) :
            overflow_handler(overflow_handler), max_size(max_size)
        {}
    };

    template<class... ValueConstructorArgs>
    explicit fifo_mesh_container(mesh_t const& mesh,
                                 overflow_t const& ovfl,
                                 ValueConstructorArgs && ...vc_args
                                ) :
    mesh(mesh), ovfl(ovfl), default_value(vc_args...), front_at(0)
    {}

    const mesh_t get_mesh(void) const {
        return mesh;
    }

    struct pair_t {
        typename mesh_t::point_t mesh_point;
        value_type & value;
    };

    // Iterator
    class iterator : 
    public boost::iterator_facade<iterator,pair_t,boost::single_pass_traversal_tag,pair_t>
    {
    public:
        iterator(fifo_mesh_container & container, typename mesh_t::const_iterator mesh_it) :
        container(container), mesh_it(mesh_it)
        {}

    private:
        typename mesh_t::const_iterator mesh_it;
        fifo_mesh_container & container;

        friend class boost::iterator_core_access;
        inline bool equal(iterator const& other) const { return mesh_it == other.mesh_it; }
        inline void increment() {
            ++mesh_it;
            container.update_fifo(mesh_it);
        }
        inline pair_t dereference() const {
            return {*mesh_it,container.fifo.back()};
        }
    };

    iterator begin() noexcept {
        fifo.clear();
        if(mesh.size()>0) fifo.emplace_back(default_value);
        return iterator(*this,std::begin(mesh));
    }

    iterator end() noexcept {
        return iterator(*this,std::end(mesh));
    }

    size_type size() { return fifo.size(); }

    template<typename ProcessorType>
    void process_front(ProcessorType processor, size_type n = 1) {
        for(int i=0; i<n && fifo.size()>0; ++i){
            processor({front_at,mesh[front_at]},fifo.front());
            fifo.pop_front();
            ++front_at;
        }
    }

private:
    void update_fifo(typename mesh_t::const_iterator const& mesh_it)
    {
        if(mesh_it != std::end(mesh)){
            fifo.emplace_back(default_value);
            while(fifo.size()>ovfl.max_size) ovfl.overflow_handler(*this);
        } else {
            while(fifo.size()>0) ovfl.overflow_handler(*this);
        }
    }

    friend class iterator;

    mesh_t mesh;
    std::deque<value_type> fifo;
    size_type front_at;
    const overflow_t ovfl;
    const value_type default_value;
};

}