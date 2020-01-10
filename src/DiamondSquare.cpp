#include "DiamondSquare.hpp"

#include <omp.h>
#include <iostream>
#include <math.h>

#include "ParallelRNG.hpp"

// Returns the height, based on the supplied average height, and factored down by the current subdivision depth
double randHeight(ParallelRNG& rng, double avg, double range, size_t subdivision) {
    return avg + rng.getDoublePlusMinus() * (range / (1 << subdivision) /*pow(2, subdivision)*/);
}

#define RANGE 0.5
void initializeMap(CubeWorld& map, size_t subdivisions, ParallelRNG& rng) {
    //initialize map
    //face 4 corners
    double face4LowerLeft = randHeight(rng, 0.5, RANGE, 0);
    double face4LowerRight = randHeight(rng, 0.5, RANGE, 0);
    double face4UpperLeft = randHeight(rng, 0.5, RANGE, 0);
    double face4UpperRight = randHeight(rng, 0.5, RANGE, 0);
    //face 5 corners
    double face5LowerLeft = randHeight(rng, 0.5, RANGE, 0);
    double face5LowerRight = randHeight(rng, 0.5, RANGE, 0);
    double face5UpperLeft = randHeight(rng, 0.5, RANGE, 0);
    double face5UpperRight = randHeight(rng, 0.5, RANGE, 0);

    map.set(0, 0, 0, face5UpperLeft);
    map.set(0, map.getSize() - 1, 0, face5UpperRight);
    map.set(0, 0, map.getSize() - 1, face4LowerLeft);
    map.set(0, map.getSize() - 1, map.getSize() - 1, face4LowerRight);

    map.set(1, 0, 0, face5UpperRight);
    map.set(1, map.getSize() - 1, 0, face5LowerRight);
    map.set(1, 0, map.getSize() - 1, face4LowerRight);
    map.set(1, map.getSize() - 1, map.getSize() - 1, face4UpperRight);

    map.set(2, 0, 0, face5LowerRight);
    map.set(2, map.getSize() - 1, 0, face5LowerLeft);
    map.set(2, 0, map.getSize() - 1, face4UpperRight);
    map.set(2, map.getSize() - 1, map.getSize() - 1, face4UpperLeft);

    map.set(3, 0, 0, face5LowerLeft);
    map.set(3, map.getSize() - 1, 0, face5UpperLeft);
    map.set(3, 0, map.getSize() - 1, face4UpperLeft);
    map.set(3, map.getSize() - 1, map.getSize() - 1, face4LowerLeft);

    map.set(4, 0, 0, face4LowerLeft);
    map.set(4, map.getSize() - 1, 0, face4LowerRight);
    map.set(4, 0, map.getSize() - 1, face4UpperLeft);
    map.set(4, map.getSize() - 1, map.getSize() - 1, face4UpperRight);

    map.set(5, 0, 0, face5LowerLeft);
    map.set(5, map.getSize() - 1, 0, face5LowerRight);
    map.set(5, 0, map.getSize() - 1, face5UpperLeft);
    map.set(5, map.getSize() - 1, map.getSize() - 1, face5UpperRight);
}

// Performs the square step. This is fairly straightforward.
void squareStep(
    CubeWorld& map,
    ParallelRNG& rng,
    char face,
    size_t row,
    size_t column,
    size_t sideLength,
    size_t sideSquares,
    double range,
    size_t subdivision)
{
    double sumAltitude = map.get(face, column*sideLength, row*sideLength);
    sumAltitude       += map.get(face, column*sideLength, (row+1)*sideLength);
    sumAltitude       += map.get(face, (column+1)*sideLength, row*sideLength);
    sumAltitude       += map.get(face, (column+1)*sideLength, (row+1)*sideLength);

    double avgAltitude = sumAltitude / 4;

    map.set(
        face,
        column*sideLength + sideLength/2, //diamondX
        row*sideLength + sideLength/2,    //diamondY
        randHeight(rng, avgAltitude, range, subdivision)
    );
}

void diamondStep (
    CubeWorld& map,
    ParallelRNG& rng,
    char face,
    size_t row,
    size_t column,
    size_t sideLength,
    size_t sideSquares,
    double range,
    size_t subdivision
) {
    size_t halfEdge = sideLength / 2;
    size_t diamondX = column*sideLength + halfEdge;
    size_t diamondY = row*sideLength + halfEdge;

    size_t squareX, squareY;
    double sumAltitude, setAltitude;

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
    #define X 0
    #define Y 1
    signed char multiples[4][2] = {{0, 1}, {-1, 0}, {0, -1}, {1, 0}};
    /*
                 (0, 1)  _
                    *   |\
             /      -     \
           |_       -
                    -
    (-1, 0) * - - - * - - - * (1, 0)
                    -       _
                    -       /|
            \       -      /
             _|     *
                 (0, -1)
    */
    #define UP    0
    #define LEFT  1
    #define DOWN  2
    #define RIGHT 3
    for (size_t direction = 0; direction < 4; direction++) {
        squareX = diamondX + halfEdge*multiples[direction][X];
        squareY = diamondY + halfEdge*multiples[direction][Y];

        sumAltitude  = 0;

        if (direction == UP && row == sideSquares - 1) {
            sumAltitude += map.get(face, squareX, squareY + halfEdge + 1);
        } else {
            sumAltitude += map.get(face, squareX, squareY + halfEdge);
        }
        if (direction == LEFT && column == 0) {
            sumAltitude += map.get(face, squareX - halfEdge - 1, squareY);
        } else {
            sumAltitude += map.get(face, squareX - halfEdge, squareY);
        }
        if (direction == DOWN && row == 0) {
            sumAltitude += map.get(face, squareX, squareY - halfEdge - 1);
        } else {
            sumAltitude += map.get(face, squareX, squareY - halfEdge);
        }
        if (direction == RIGHT && column == sideSquares - 1) {
            sumAltitude += map.get(face, squareX + halfEdge + 1, squareY);
        } else {
            sumAltitude += map.get(face, squareX + halfEdge, squareY);
        }

        map.set(
            face,
            squareX,
            squareY,
            setAltitude = randHeight(rng, sumAltitude / 4, range, subdivision)
        );

        if (direction == UP && row == sideSquares - 1) {
            map.set(face, squareX, squareY + 1, setAltitude);
        }
        if (direction == LEFT && column == 0) {
            map.set(face, squareX - 1, squareY, setAltitude);
        }
        if (direction == DOWN && row == 0) {
            map.set(face, squareX, squareY - 1, setAltitude);
        }
        if (direction == RIGHT && column == sideSquares - 1) {
            map.set(face, squareX + 1, squareY, setAltitude);
        }
    }
}

// This is the serial version of the main algorithm. See the parallel version for
// complete commentary.
void diamondSquare(CubeWorld& map) {
    double t1, t2;

    ParallelRNG rng;
    rng.reinitialize();

    //number of times the map will be subdivided
    size_t subdivisions = (size_t) log2(map.getSize() - 1);

    initializeMap(map, subdivisions, rng);

    for (size_t subdivision = 0; subdivision < subdivisions; subdivision++) {
        size_t sideSquares = pow(2,subdivision); /* 1 << n */
        size_t sideLength = (map.getSize() - 1) / sideSquares; /* (length - 1) >> n */

        rng.reinitialize();

        t1 = omp_get_wtime();

        for (size_t face = 0; face < 6; face++) {
            for (size_t row = 0; row < sideSquares; row++) {
                for (size_t column = 0; column < sideSquares; column++) {
                    squareStep(map, rng, face, row, column, sideLength, sideSquares, RANGE, subdivision);
                }
            }
            for (size_t row = 0; row < sideSquares; row++) {
                for (size_t column = 0; column < sideSquares; column++) {
                    diamondStep(map, rng, face, row, column, sideLength, sideSquares, RANGE, subdivision);
                }
            }
        }

        t2 = omp_get_wtime();
        std::cout
            << std::string("Squares: ") + std::to_string(pow(sideSquares, 2) * 6) + std::string("\t")
            << std::to_string((pow(sideSquares, 2) * 6) / (t2 - t1)) + std::string(" squares/second\t")
            << std::string("Generators: ") + std::to_string(rng.getNumGenerators())
            << std::endl
        ;

    }
}

// This is the parallel version of the algorithm.
// void diamondSquareParallel(Array2D<double>& map, size_t length) {
//     double t1, t2;

//     ParallelRNG rng;
//     rng.reinitialize();

//     //number of times the map will be subdivided
//     size_t subdivisions = (size_t) log2(length - 1);

//     //initialize map
//     map.set(0, 0, randHeight(rng, 0.5, RANGE, 0));
//     map.set(0, length - 1, randHeight(rng, 0.5, RANGE, 0));
//     map.set(length - 1, 0, randHeight(rng, 0.5 , RANGE, 0));
//     map.set(length - 1, length - 1, randHeight(rng, 0.5 , RANGE, 0));

//     #pragma omp parallel shared(map)
//     {
    
//     #pragma omp single
//     {
//     // The random number generator must be initialized in the same parallel region
//     // that it will be used in, before it is used. Otherwise it may not have the
//     // appropriate number of generators.
//     rng.reinitialize();
//     }

//     for (size_t subdivision = 0; subdivision < subdivisions; subdivision++) {
//         size_t sideSquares = pow(2,subdivision); // The number of cells in each row and column
//         size_t sideLength = (length - 1) / sideSquares; // The number of verticies on each cell's side

//         #pragma omp single
//         {
//         t1 = omp_get_wtime();
//         }

//         #pragma omp for schedule(auto)
//         for (size_t row = 0; row < sideSquares; row++) {
//             for (size_t column = 0; column < sideSquares; column++) {
//                 squareStep(map, rng, row, column, sideLength, RANGE, subdivision);
//             }
//         }
//         // These loops cannot be combined because the results of the diamond step depend on the results of
//         // the square step.
//         #pragma omp for schedule(auto)
//         for (size_t row = 0; row < sideSquares; row++) {
//             for (size_t column = 0; column < sideSquares; column++) {
//                 diamondStep(map, rng, row, column, sideLength, sideSquares, RANGE, subdivision);
//             }
//         }

//         #pragma omp single
//         {
//         t2 = omp_get_wtime();
//         std::cout
//             << std::string("Squares: ") + std::to_string(pow(sideSquares, 2)) + std::string("\t")
//             << std::to_string((pow(sideSquares, 2)) / (t2 - t1)) + std::string(" squares/second\t")
//             << std::string("Generators: ") + std::to_string(rng.getNumGenerators())
//             << std::endl
//         ;
//         }
//     }
//     }
// }