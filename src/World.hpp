#include <stdlib.h>

#include "ParallelBlock.hpp"
#include "Unit.hpp"

template<class T> class World {
    private:
        
    public:
        World(int resolution);

        friend ParallelBlock<T>;
        
        /*
         * Level 0 gets a parallel block for serial execution
         * 
         * @param[in,out] blocks : The number of blocks returned
        */
        virtual std::vector<ParallelBlock> getSubdivision(int level, size_t* blocks);
};