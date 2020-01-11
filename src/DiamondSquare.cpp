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
    size_t row,
    size_t column,
    size_t sideLength,
    size_t sideSquares,
    double range,
    size_t subdivision)
{
    for (char face = 0; face < 6; face++)
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
}

void diamondStep (
    CubeWorld& map,
    ParallelRNG& rng,
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
    double sumAltitude[6], setAltitude[6];

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

        sumAltitude[0] = 0;
        sumAltitude[1] = 0;
        sumAltitude[2] = 0;
        sumAltitude[3] = 0;
        sumAltitude[4] = 0;
        sumAltitude[5] = 0;

        if (direction == UP && row == sideSquares - 1) {
            sumAltitude[0] += map.get(0, squareX, squareY + halfEdge + 1);
            sumAltitude[1] += map.get(1, squareX, squareY + halfEdge + 1);
            sumAltitude[2] += map.get(2, squareX, squareY + halfEdge + 1);
            sumAltitude[3] += map.get(3, squareX, squareY + halfEdge + 1);
            sumAltitude[4] += map.get(4, squareX, squareY + halfEdge + 1);
            sumAltitude[5] += map.get(5, squareX, squareY + halfEdge + 1);
        } else {
            sumAltitude[0] += map.get(0, squareX, squareY + halfEdge);
            sumAltitude[1] += map.get(1, squareX, squareY + halfEdge);
            sumAltitude[2] += map.get(2, squareX, squareY + halfEdge);
            sumAltitude[3] += map.get(3, squareX, squareY + halfEdge);
            sumAltitude[4] += map.get(4, squareX, squareY + halfEdge);
            sumAltitude[5] += map.get(5, squareX, squareY + halfEdge);
        }
        if (direction == LEFT && column == 0) {
            sumAltitude[0] += map.get(0, squareX - halfEdge - 1, squareY);
            sumAltitude[1] += map.get(1, squareX - halfEdge - 1, squareY);
            sumAltitude[2] += map.get(2, squareX - halfEdge - 1, squareY);
            sumAltitude[3] += map.get(3, squareX - halfEdge - 1, squareY);
            sumAltitude[4] += map.get(4, squareX - halfEdge - 1, squareY);
            sumAltitude[5] += map.get(5, squareX - halfEdge - 1, squareY);
        } else {
            sumAltitude[0] += map.get(0, squareX - halfEdge, squareY);
            sumAltitude[1] += map.get(1, squareX - halfEdge, squareY);
            sumAltitude[2] += map.get(2, squareX - halfEdge, squareY);
            sumAltitude[3] += map.get(3, squareX - halfEdge, squareY);
            sumAltitude[4] += map.get(4, squareX - halfEdge, squareY);
            sumAltitude[5] += map.get(5, squareX - halfEdge, squareY);
        }
        if (direction == DOWN && row == 0) {
            sumAltitude[0] += map.get(0, squareX, squareY - halfEdge - 1);
            sumAltitude[1] += map.get(1, squareX, squareY - halfEdge - 1);
            sumAltitude[2] += map.get(2, squareX, squareY - halfEdge - 1);
            sumAltitude[3] += map.get(3, squareX, squareY - halfEdge - 1);
            sumAltitude[4] += map.get(4, squareX, squareY - halfEdge - 1);
            sumAltitude[5] += map.get(5, squareX, squareY - halfEdge - 1);
        } else {
            sumAltitude[0] += map.get(0, squareX, squareY - halfEdge);
            sumAltitude[1] += map.get(1, squareX, squareY - halfEdge);
            sumAltitude[2] += map.get(2, squareX, squareY - halfEdge);
            sumAltitude[3] += map.get(3, squareX, squareY - halfEdge);
            sumAltitude[4] += map.get(4, squareX, squareY - halfEdge);
            sumAltitude[5] += map.get(5, squareX, squareY - halfEdge);
        }
        if (direction == RIGHT && column == sideSquares - 1) {
            sumAltitude[0] += map.get(0, squareX + halfEdge + 1, squareY);
            sumAltitude[1] += map.get(1, squareX + halfEdge + 1, squareY);
            sumAltitude[2] += map.get(2, squareX + halfEdge + 1, squareY);
            sumAltitude[3] += map.get(3, squareX + halfEdge + 1, squareY);
            sumAltitude[4] += map.get(4, squareX + halfEdge + 1, squareY);
            sumAltitude[5] += map.get(5, squareX + halfEdge + 1, squareY);
        } else {
            sumAltitude[0] += map.get(0, squareX + halfEdge, squareY);
            sumAltitude[1] += map.get(1, squareX + halfEdge, squareY);
            sumAltitude[2] += map.get(2, squareX + halfEdge, squareY);
            sumAltitude[3] += map.get(3, squareX + halfEdge, squareY);
            sumAltitude[4] += map.get(4, squareX + halfEdge, squareY);
            sumAltitude[5] += map.get(5, squareX + halfEdge, squareY);
        }

        map.set(
            0,
            squareX,
            squareY,
            setAltitude[0] = randHeight(rng, sumAltitude[0] / 4, range, subdivision)
        );
        map.set(
            1,
            squareX,
            squareY,
            setAltitude[1] = randHeight(rng, sumAltitude[1] / 4, range, subdivision)
        );
        map.set(
            2,
            squareX,
            squareY,
            setAltitude[2] = randHeight(rng, sumAltitude[2] / 4, range, subdivision)
        );
        map.set(
            3,
            squareX,
            squareY,
            setAltitude[3] = randHeight(rng, sumAltitude[3] / 4, range, subdivision)
        );
        map.set(
            4,
            squareX,
            squareY,
            setAltitude[4] = randHeight(rng, sumAltitude[4] / 4, range, subdivision)
        );
        map.set(
            5,
            squareX,
            squareY,
            setAltitude[5] = randHeight(rng, sumAltitude[5] / 4, range, subdivision)
        );

        if (direction == UP && row == sideSquares - 1) {
            map.set(0, squareX, squareY + 1, setAltitude[0]);
            map.set(1, squareX, squareY + 1, setAltitude[1]);
            map.set(2, squareX, squareY + 1, setAltitude[2]);
            map.set(3, squareX, squareY + 1, setAltitude[3]);
            map.set(4, squareX, squareY + 1, setAltitude[4]);
            map.set(5, squareX, squareY + 1, setAltitude[5]);
        }
        if (direction == LEFT && column == 0) {
            map.set(0, squareX - 1, squareY, setAltitude[0]);
            map.set(1, squareX - 1, squareY, setAltitude[1]);
            map.set(2, squareX - 1, squareY, setAltitude[2]);
            map.set(3, squareX - 1, squareY, setAltitude[3]);
            map.set(4, squareX - 1, squareY, setAltitude[4]);
            map.set(5, squareX - 1, squareY, setAltitude[5]);
        }
        if (direction == DOWN && row == 0) {
            map.set(0, squareX, squareY - 1, setAltitude[0]);
            map.set(1, squareX, squareY - 1, setAltitude[1]);
            map.set(2, squareX, squareY - 1, setAltitude[2]);
            map.set(3, squareX, squareY - 1, setAltitude[3]);
            map.set(4, squareX, squareY - 1, setAltitude[4]);
            map.set(5, squareX, squareY - 1, setAltitude[5]);
        }
        if (direction == RIGHT && column == sideSquares - 1) {
            map.set(0, squareX + 1, squareY, setAltitude[0]);
            map.set(1, squareX + 1, squareY, setAltitude[1]);
            map.set(2, squareX + 1, squareY, setAltitude[2]);
            map.set(3, squareX + 1, squareY, setAltitude[3]);
            map.set(4, squareX + 1, squareY, setAltitude[4]);
            map.set(5, squareX + 1, squareY, setAltitude[5]);
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

        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                squareStep(map, rng, row, column, sideLength, sideSquares, RANGE, subdivision);
            }
        }
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                diamondStep(map, rng, row, column, sideLength, sideSquares, RANGE, subdivision);
            }
        }

        t2 = omp_get_wtime();
        std::cout
            << std::string("Squares: ") + std::to_string(((size_t)pow(sideSquares, 2)) * 6) + std::string("\t")
            << std::to_string((pow(sideSquares, 2) * 6) / (t2 - t1)) + std::string(" squares/second\t")
            << std::string("Generators: ") + std::to_string(rng.getNumGenerators())
            << std::endl
        ;

    }
}

void diamondSquareParallel(CubeWorld& map) {
    double t1, t2;

    ParallelRNG rng;
    rng.reinitialize();

    //number of times the map will be subdivided
    size_t subdivisions = (size_t) log2(map.getSize() - 1);

    initializeMap(map, subdivisions, rng);

    #pragma omp parallel
    {
    for (size_t subdivision = 0; subdivision < subdivisions; subdivision++) {
        size_t sideSquares = pow(2,subdivision); /* 1 << n */
        size_t sideLength = (map.getSize() - 1) / sideSquares; /* (length - 1) >> n */

        #pragma omp single
        {
        rng.reinitialize();
        }

        t1 = omp_get_wtime();

        #pragma omp for
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                squareStep(map, rng, row, column, sideLength, sideSquares, RANGE, subdivision);
            }
        }
        #pragma omp for
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                diamondStep(map, rng, row, column, sideLength, sideSquares, RANGE, subdivision);
            }
        }

        #pragma omp single
        {
            t2 = omp_get_wtime();
            std::cout
                << std::string("Squares: ") + std::to_string(((size_t)pow(sideSquares, 2)) * 6) + std::string("\t")
                << std::to_string((pow(sideSquares, 2) * 6) / (t2 - t1)) + std::string(" squares/second\t")
                << std::string("Generators: ") + std::to_string(rng.getNumGenerators())
                << std::endl
            ;
        }
    }
    }
}