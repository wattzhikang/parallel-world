#ifndef DIAMOND_SQUARE
#define DIAMOND_SQUARE

#include "CubeWorld.hpp"
#include "ParallelRNG.hpp"

void diamondSquare(CubeWorld& map, unsigned long seed);
void diamondSquareParallel(CubeWorld& map, unsigned long seed);

#endif