#include <cstdlib>
#include <cmath>
#include <sstream>

#include "uniform_mesh.hpp"
#include "fifo_mesh_container.hpp"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

using namespace realevol;

template<class T> T sqr(T x) { return x*x; } 

void fifo_processor(double t, double const& x)
{
    std::cout << "Processing element " << x << " at " << t << std::endl;
}

void run_test(std::size_t max_size, std::size_t proc_block_size, bool change_elements)
{
    std::cout << "max_size = " << max_size << ", ";
    std::cout << "proc_block_size = " << proc_block_size << ", ";
    std::cout << "change_elements = " << change_elements << std::endl;

    uniform_mesh<> m(0,19.0,20);
    fifo_mesh_container<int,uniform_mesh<>> f(m,{fifo_processor,max_size,proc_block_size},999/* default_value*/);

    int n = 0;
    for(auto it = f.arg_value_begin(); it != f.arg_value_end(); ++it){
        std::cout << "Adding element ";
        if(change_elements){
            std::cout << n;
            it->get<1>() = n;
        } else 
            std::cout << 999;
        std::cout << " at " << it->get<0>()  << std::endl;
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