///@file

#include "DiamondSquare.hpp"

#include <omp.h>
#include <iostream>
#include <math.h>


/**
 * Returns the height, based on the supplied average height, and factored down by the current subdivision depth
 * @param rng The random number generator, to supply randomness
 * @param avg The average of the surrounding terrain, however you define it
 * @param range To specify the maximum amount by which to displace the terrain
 * @param subdivision This function reduces the range by a factor of 2 for every subdivison
*/
double randHeight(ParallelRNGSequence* rng, double avg, double range, size_t subdivision) {
    return avg + rng->getDoublePlusMinus() * (range / (1 << subdivision) /*pow(2, subdivision)*/);
}

#define RANGE 0.5

/**
 * Initializes the 4 corners on each face. All adjacent corner pixels are given the same value
 * @param map The map
 * @param subdivisions
 * @param rng The random number generator
*/
void initializeMap(CubeWorld& map, size_t subdivisions, ParallelRNGSequence* rng) {
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

/** Performs the square step. This is fairly straightforward.
 * @param map The map to manipulate
 * @param rng The random number generator sequence to use (see ParallelRNG)
 * @param row The row being manipulated, according to the current subdivision's partition scheme
 * @param column The column being manipulated, according to the current subdivision's partition scheme
 * @param sideLength The size of a single square, in pixels, according to the current subdivision's partition scheme
 * @param sideSquares The numer of vertical or horizontal squares the face has been divided into
 * @param range See randHeight()
 * @param subdivision The current subdivisions
*/
void squareStep(
    CubeWorld& map,
    ParallelRNGSequence* rng,
    size_t row,
    size_t column,
    size_t sideLength,
    size_t sideSquares,
    double range,
    size_t subdivision)
{
    //Does the same thing to all 6 faces of the cube
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

/**
 * Performs the diamond step
 * @param map The map to manipulate
 * @param rng The random number generator to use
 * @param row The row being manipulated, according to the current subdivision's partition scheme
 * @param column The column being manipulated, according to the current subdivision's partition scheme
 * @param sideLength The size of a single square, in pixels, according to the current subdivision's partition scheme
 * @param sideSquares The numer of vertical or horizontal squares the face has been divided into
 * @param range See randHeight()
 * @param subdivision The current subdivisions
*/
void diamondStep (
    CubeWorld& map,
    ParallelRNGSequence* rng,
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

        /*
         * This entire sequence is unfolded, mostly because it wouldn't look any
         * better to do this with loops. The reason why there is so much repetition
         * here is because it wouldn't look any better to use loops.
        */

        sumAltitude[0] = 0;
        sumAltitude[1] = 0;
        sumAltitude[2] = 0;
        sumAltitude[3] = 0;
        sumAltitude[4] = 0;
        sumAltitude[5] = 0;

        //if this square side is on the edge of a face, then the pixel it needs
        //to read is one pixel higher than if it were not on the edge. This
        //pattern repeats for the upper, left, right, and bottom face edges.
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

        //If the manipulated square is on an edge, then the pixel adjacent to it on
        //the adjacent face should be set to the same value. This pattern repeats
        //for all edges of the face.
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

/**
 * The serial version of the Diamond Square algorithm
*/
//For code commentary, see the parallel version
void diamondSquare(CubeWorld& map, unsigned long seed) {
    double t1, t2;

    //number of times the map will be subdivided
    size_t subdivisions = (size_t) log2(map.getSize() - 1);

    ParallelRNG rng(pow(2,subdivisions - 1), seed);

    initializeMap(map, subdivisions, rng.getSequence(0));

    for (size_t subdivision = 0; subdivision < subdivisions; subdivision++) {
        size_t sideSquares = pow(2,subdivision); /* 1 << n */
        size_t sideLength = (map.getSize() - 1) / sideSquares; /* (length - 1) >> n */

        t1 = omp_get_wtime();

        // For consistency, the parallel and serial versions of this algorithm access the random
        // number generator the same way.
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                squareStep(map, rng.getSequence(row), row, column, sideLength, sideSquares, RANGE, subdivision);
            }
        }
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                diamondStep(map, rng.getSequence(row), row, column, sideLength, sideSquares, RANGE, subdivision);
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

/**
 * The parallel version of the Diamond Square algorithm.
*/
void diamondSquareParallel(CubeWorld& map, unsigned long seed) {
    //To store timings for measuring performance
    double t1, t2;

    //number of times the map will be subdivided
    size_t subdivisions = (size_t) log2(map.getSize() - 1);

    ParallelRNG rng(pow(2,subdivisions - 1), seed);

    //initialize the corners of the faces
    initializeMap(map, subdivisions, rng.getSequence(0));

    #pragma omp parallel
    {
    /* This algorithm is run one subdivision at a time. The first subdivision divides a face into one square;
     * the second divides a face into 4 squares; the third into 9 squares; and so on.
    */
    for (size_t subdivision = 0; subdivision < subdivisions; subdivision++) {
        //the number of squares on each side of a face in the current subdivision level
        size_t sideSquares = pow(2,subdivision); /* 1 << n */
        //the length, in pixels, of each square that the faces will be divided into
        size_t sideLength = (map.getSize() - 1) / sideSquares; /* (length - 1) >> n */

        t1 = omp_get_wtime(); //timing for performance

        /* These loops run the diamond and square steps of the algorithm on the rows and colums of squares
         * in the subdivision level. Since there can be as many as one parallel thread for each row of
         * squares, the random number generators are accessed accordingly.
        */
        #pragma omp for
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                squareStep(map, rng.getSequence(row), row, column, sideLength, sideSquares, RANGE, subdivision);
            }
        }
        //these steps cannot be folded into a single loop because the diamond step depends on
        //information generated in the square step
        #pragma omp for
        for (size_t row = 0; row < sideSquares; row++) {
            for (size_t column = 0; column < sideSquares; column++) {
                diamondStep(map, rng.getSequence(row), row, column, sideLength, sideSquares, RANGE, subdivision);
            }
        }

        //Outputs some profiling information. This should probably be removed at some point.
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