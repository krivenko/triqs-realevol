#include <cstdlib>
#include <cmath>
#include <sstream>

#include "uniform_mesh.hpp"
#include "fifo_mesh_container.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

using namespace realevol;

template<class T> T sqr(T x) { return x*x; } 

using mesh_t = uniform_mesh<>;
using cont_t = fifo_mesh_container<int,mesh_t>;

std::size_t pbs;
void overflow_handler(cont_t & f)
{
    auto proc = [](mesh_t::deref_result_t mp, cont_t::value_type v){
        std::cout << "Processing element " << v << " at " << mp.value << std::endl;
    };
    f.process_front(proc,pbs);
}

void run_test(std::size_t max_size, std::size_t proc_block_size, bool change_elements)
{
    pbs = proc_block_size;
    std::cout << "max_size = " << max_size << ", ";
    std::cout << "proc_block_size = " << proc_block_size << ", ";
    std::cout << "change_elements = " << change_elements << std::endl;

    mesh_t m(0,19.0,20);
    cont_t f(m,{overflow_handler,max_size},999/* default_value*/);

    int n = 0;
    for(auto it = std::begin(f); it != std::end(f); ++it){
        std::cout << "Adding element ";
        if(change_elements){
            std::cout << n;
            it->value = n;
        } else 
            std::cout << 999;
        std::cout << " at " << it->mesh_point.value  << std::endl;
        ++n;
    }
}

int main()
{
    run_test(1,1,false);
    run_test(1,1,true);
    run_test(10,5,false);
    run_test(10,5,true);

    return EXIT_SUCCESS;
}