#include "DiamondSquare.hpp"

#include <omp.h>
#include <iostream>
#include <math.h>

#include "ParallelRNG.hpp"

// Returns the height, based on the supplied average height, and factored down by the current subdivision depth
double randHeight(ParallelRNG& rng, double avg, double range, size_t subdivision) {
    return avg + rng.getDoublePlusMinus() * (range / (1 << subdivision) /*pow(2, subdivision)*/);
}

// Performs the square step. This is fairly straightforward.
void squareStep(
    Array2D<double>& map,
    ParallelRNG& rng,
    size_t row,
    size_t column,
    size_t sideLength,
    double range,
    size_t subdivision)
{
    double sumAltitude = map(column*sideLength, row*sideLength);
    sumAltitude          += map(column*sideLength, (row+1)*sideLength);
    sumAltitude          += map((column+1)*sideLength, row*sideLength);
    sumAltitude          += map((column+1)*sideLength, (row+1)*sideLength);

    double avgAltitude = sumAltitude / 4;

    map.set(
        column*sideLength + sideLength/2, //diamondX
        row*sideLength + sideLength/2,    //diamondY
        randHeight(rng, avgAltitude, range, subdivision)
    );
}

#define X 0
#define Y 1
void diamondStep (
    Array2D<double>& map,
    ParallelRNG& rng,
    size_t row,
    size_t column,
    size_t sideLength,
    size_t sideSquares,
    double range,
    size_t subdivision)
{
    size_t halfEdge = sideLength / 2;
    size_t diamondX = column*sideLength + halfEdge;
    size_t diamondY = row*sideLength + halfEdge;

    size_t squareX, squareY;
    double sumAltitude;
    size_t numVerticies;

    /*
        This is a bit convoluted. The whole point is to calculate the locations of the upper-middle,
        left-middle, lower-middle, and right-middle verticies---and all the locations of their
        surrounding verticies, without having to hard-code all of them. The basic relationship between
        them is that the locations of the cells are all one half edge away from the center of the
        cell. The upper-middle vertex is at the same X coordinate as the center, but one half edge
        up. The left-middle vertex is at the same Y coordinate, but one half edge left. This hard-coded
        array records all these relationships as multiples, and the loop simply multiplies by a
        different pair in every iteration of the loop.
    */
    signed char multiples[4][2] = {{0, -1}, {-1, 0}, {0, 1}, {1, 0}};
    /*
                 (0,-1)  _
                    *   |\
             /      -     \
           |_       -
                    -
    (-1, 0) * - - - * - - - * (1, 0)
                    -       _
                    -       /|
            \       -      /
             _|     *
                 (0, 1)
    */
    for (size_t i = 0; i < 4; i++) {
        squareX = diamondX + halfEdge*multiples[i][X];
        squareY = diamondY + halfEdge*multiples[i][Y];

        sumAltitude  = 0;
        numVerticies = 0;

        if (row > 0 || i != 0) {
            sumAltitude += map(squareX, squareY - halfEdge);
            numVerticies++;
        }
        if (column > 0 || i != 1) {
            sumAltitude += map(squareX - halfEdge, squareY);
            numVerticies++;
        }
        if (row < sideSquares - 1 || i != 2) {
            sumAltitude += map(squareX, squareY + halfEdge);
            numVerticies++;
        }
        if (column < sideSquares - 1 || i != 3) {
            sumAltitude += map(squareX + halfEdge, squareY);
            numVerticies++;
        }

        map.set(
            squareX,
            squareY,
            randHeight(rng, sumAltitude / numVerticies, range, subdivision)
        );
    }
}

// This is the serial version of the main algorithm. See the parallel version for
// complete commentary.
#define RANGE 0.5
void diamondSquare(Array2D<double>& map, size_t length) {
    double t1, t2;

    ParallelRNG rng;
    rng.reinitialize();

    //number of times the map will be subdivided
    size_t subdivisions = (size_t) log2(length - 1);


    //initialize map
    map.set(0, 0, randHeight(rng, 0.5, RANGE, 0));
    map.set(0, length - 1, randHeight(rng, 0.5, RANGE, 0));
    map.set(length - 1, 0, randHeight(rng, 0.5 , RANGE, 0));
    map.set(length - 1, length - 1, randHeight(rng, 0.5 , RANGE, 0));

    for (size_t subdivision = 0; subdivision < subdivisions; subdivision++) {
        size_t sideSquares = pow(2,subdivision); /* 1 << n */
        size_t sideLength = (length - 1) / sideSquares; /* (length - 1) >> n */

        t1 = omp_get_wtime();

        rng.reinitialize();

        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                squareStep(map, rng, row, column, sideLength, RANGE, subdivision);
            }
        }
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                diamondStep(map, rng, row, column, sideLength, sideSquares, RANGE, subdivision);
            }
        }

        t2 = omp_get_wtime();
        std::cout
            << std::string("Squares: ") + std::to_string(pow(sideSquares, 2)) + std::string("\t")
            << std::to_string((pow(sideSquares, 2)) / (t2 - t1)) + std::string(" squares/second\t")
            << std::string("Generators: ") + std::to_string(rng.getNumGenerators())
            << std::endl
        ;
    }
}

// This is the parallel version of the algorithm.
void diamondSquareParallel(Array2D<double>& map, size_t length) {
    double t1, t2;

    ParallelRNG rng;
    rng.reinitialize();

    //number of times the map will be subdivided
    size_t subdivisions = (size_t) log2(length - 1);

    //initialize map
    map.set(0, 0, randHeight(rng, 0.5, RANGE, 0));
    map.set(0, length - 1, randHeight(rng, 0.5, RANGE, 0));
    map.set(length - 1, 0, randHeight(rng, 0.5 , RANGE, 0));
    map.set(length - 1, length - 1, randHeight(rng, 0.5 , RANGE, 0));

    #pragma omp parallel shared(map)
    {
    
    #pragma omp single
    {
    // The random number generator must be initialized in the same parallel region
    // that it will be used in, before it is used. Otherwise it may not have the
    // appropriate number of generators.
    rng.reinitialize();
    }

    for (size_t subdivision = 0; subdivision < subdivisions; subdivision++) {
        size_t sideSquares = pow(2,subdivision); // The number of cells in each row and column
        size_t sideLength = (length - 1) / sideSquares; // The number of verticies on each cell's side

        #pragma omp single
        {
        t1 = omp_get_wtime();
        }

        #pragma omp for schedule(auto)
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                squareStep(map, rng, row, column, sideLength, RANGE, subdivision);
            }
        }
        // These loops cannot be combined because the results of the diamond step depend on the results of
        // the square step.
        #pragma omp for schedule(auto)
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                diamondStep(map, rng, row, column, sideLength, sideSquares, RANGE, subdivision);
            }
        }

        #pragma omp single
        {
        t2 = omp_get_wtime();
        std::cout
            << std::string("Squares: ") + std::to_string(pow(sideSquares, 2)) + std::string("\t")
            << std::to_string((pow(sideSquares, 2)) / (t2 - t1)) + std::string(" squares/second\t")
            << std::string("Generators: ") + std::to_string(rng.getNumGenerators())
            << std::endl
        ;
        }
    }
    }
}