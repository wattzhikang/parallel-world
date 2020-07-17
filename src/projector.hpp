#ifndef PROJ_CUBE
#define PROJ_CUBE

#include "CubeWorld.hpp"

namespace projector
{
    void getSpherical(CubeWorld& map, char face, size_t x, size_t y, float *lat, float *lon);

    void getCartesian(CubeWorld& map, float lat, float lon, char *face, size_t *x, size_t *y);
} // namespace projector


#endif