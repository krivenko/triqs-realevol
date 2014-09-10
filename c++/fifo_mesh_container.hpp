#pragma once

#include <deque>
#include <functional>
#include <type_traits>
#include <stdexcept>
#include <string>

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

    static_assert(std::is_base_of<mesh_base<typename mesh_t::node_index_t,typename mesh_t::mesh_point_t,mesh_t>,mesh_t>::value,
                  "Mesh is not derived from mesh_base");

    struct processing_t {
        std::function<void(typename mesh_t::mesh_point_t, value_type const&)> fifo_processor;
        size_type max_size;
        size_type proc_block_size;

        template<class FifoProcessorType>
        processing_t(FifoProcessorType fifo_processor, size_type max_size = 1, size_type proc_block_size = 1) :
            fifo_processor(fifo_processor), max_size(max_size), proc_block_size(proc_block_size) {}
    };

    template<class... ValueConstructorArgs>
    explicit fifo_mesh_container(mesh_t const& mesh,
                                 processing_t const& proc,
                                 ValueConstructorArgs && ...vc_args
                                ) :
    mesh(mesh), proc(proc), default_value(vc_args...)
    {
        if(proc.proc_block_size > proc.max_size) throw(std::range_error(
            "proc_block_size > max_size in fifo_mesh_container (" +
            std::to_string(proc.proc_block_size) + ">" +
            std::to_string(proc.max_size) + ")."
        ));
    }

    const mesh_t get_mesh(void) const
    {
        return mesh;
    }

    // Iterator
    using iterator_deref_type = boost::tuple<typename mesh_t::mesh_point_t const,value_type&>;

    class arg_value_iterator : 
    public boost::iterator_facade<
        arg_value_iterator,
        iterator_deref_type,
        boost::single_pass_traversal_tag,
        iterator_deref_type>
    {
    public:
        arg_value_iterator(fifo_mesh_container & container, typename mesh_t::const_iterator mesh_it) :
        container(container), mesh_it(mesh_it)
        {
        }

    private:
        typename mesh_t::const_iterator mesh_it;
        fifo_mesh_container & container;

        friend class boost::iterator_core_access;
        inline bool equal(arg_value_iterator const& other) const
        {
            return mesh_it == other.mesh_it;
        }
        inline void increment()
        {
            ++mesh_it;
            container.update_fifo(mesh_it);
        }
        inline iterator_deref_type dereference() const
        {
            return iterator_deref_type(mesh_it->value,container.fifo.back());
        }
    };

    arg_value_iterator arg_value_begin(void) noexcept
    {
        fifo.clear();
        if(mesh.size()>0) fifo.emplace_back(default_value);
        return arg_value_iterator(*this,std::begin(mesh));
    }

    arg_value_iterator arg_value_end(void) noexcept
    {
        return arg_value_iterator(*this,std::end(mesh));
    }

private:
    void update_fifo(typename mesh_t::const_iterator const& mesh_it)
    {
        if(mesh_it != std::end(mesh)){
            fifo.emplace_back(default_value);
            if(fifo.size()>proc.max_size){
                for(size_type n=0; n<proc.proc_block_size; ++n){
                    proc.fifo_processor((mesh_it-proc.max_size+n)->value,fifo.front());
                    fifo.pop_front();
                }
            }
        } else {
            for(size_type n=fifo.size(); n>0; --n){
                proc.fifo_processor((mesh_it-n)->value,fifo.front());
                fifo.pop_front();
            }
        }
    }

    friend class arg_value_iterator;

    mesh_t mesh;
    std::deque<value_type> fifo;
    const processing_t proc;
    const value_type default_value;
};

}