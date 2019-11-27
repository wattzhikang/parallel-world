#include "World.hpp"

template<class T> class CubeWorld : public World {
    private:
        size_t resolution;
        void Unit<T> data[6][resolution][resolution];
    public:
        CubeWorld(size_t resolution);
};