#include <math.h>
#include <string.h>
#include <iostream>
#include <omp.h>

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>

#include "../lib/bitmap/bitmap_image.hpp"

#define u_int32_t_MAX 0xFFFFFFFF

template<typename T> class Array2D {
    private:
        size_t size;
        T *data;
    public:
        Array2D(size_t size) : size(size) {
            try {
                data = new T[size * size];
            }
            catch (std::bad_alloc &ba) {
                std::cerr << "bad_alloc caught: " << ba.what();
            }
        };
        ~Array2D() {
            delete[] data;
        }
        T operator() (size_t i, size_t j) {
            return data[i * size + j];
        }
        void set(size_t i, size_t j, T datum) {
            data[i * size + j] = datum;
        }
        void set(T toSet) {
            size_t i, sSquared = (size_t) pow(size, 2);
            for (i = 0; i < sSquared; i++)
            {
                data[i] = toSet;
            }
            
        }
        size_t getSize() {
            return size;
        }
};

/* The state word must be initialized to non-zero */
inline uint32_t xorshift32(uint32_t x)
{
	/* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return x;
}

//#define CRAPPY_GEN

class ParallelRNG {
    private:
        #ifndef CRAPPY_GEN
        gsl_rng **generators = NULL;
        #else
        uint32_t *generators;
        const double conversion = 1.0 / 0xffffffff;
        #endif
        size_t quantity = 0;
        bool activated = false;
    public:
        ~ParallelRNG() {
            deinitialize();
        }
        void deinitialize() {
            for (size_t i = 0; i < quantity; i++) {
                #ifndef CRAPPY_GEN
                gsl_rng_free(generators[i]);
                #endif
            }
            delete generators;
            
            activated = false;
        }
        void reinitialize() {
            if (activated) {
                deinitialize();
            }

            quantity = omp_get_num_threads();

            #ifndef CRAPPY_GEN
            generators = new gsl_rng*[quantity];
            #else
            generators = new uint32_t[quantity];
            #endif

            gsl_rng* seeder = gsl_rng_alloc(gsl_rng_taus2);

            for (size_t i = 0; i < quantity; i++) {
                #ifndef CRAPPY_GEN
                generators[i] = gsl_rng_alloc(gsl_rng_taus2);
                gsl_rng_set(generators[i], gsl_rng_get(seeder));
                #else
                generators[i] = (uint32_t) gsl_rng_get(seeder);
                #endif
            }

            gsl_rng_free(seeder);

            activated = true;
        }
        double getDouble() {
            if (activated) {
                #ifndef CRAPPY_GEN
                return gsl_rng_uniform(generators[omp_get_thread_num()]);
                #else
                uint32_t* addr = generators + omp_get_thread_num();
                uint32_t x = *addr;
                x ^= x << 13;
                x ^= x >> 17;
                x ^= x << 5;
                return (*addr = x) * conversion;
                #endif
            }
            return 0.0;
        }
        double getDoublePlusMinus() {
            return getDouble() * 2.0 - 1.0;
        }
        size_t getNumGenerators() {
            return quantity;
        }

};

template<typename T> void printTerrain(Array2D<T>& map) {
    size_t length = map.getSize();
    bitmap_image image(length, map.getSize());

    u_int32_t red, green, blue;

    u_int32_t oneThird = (u_int32_t) (((size_t) u_int32_t_MAX) / 3);
    u_int32_t twoThirds = (u_int32_t) ( 2 * ((size_t) u_int32_t_MAX) / 3);

    size_t x, y;
    double pixel; /* you could probably use ifdefs to get machine int implementation */

    for (x = 0; x < length; x++) {
        for (y = 0; y < length; y++) {
            /*
            pixel = map(x,y) * u_int32_t_MAX;
            if (pixel > oneThird) {
                if (pixel > twoThirds) {
                    red = 0;
                    green = 254;
                    blue = (pixel - twoThirds) / (oneThird / 254);
                } else {
                    red = 0;
                    green = (pixel - oneThird) / (oneThird / 254);
                    blue = 0;
                }
            } else {
                red = pixel / (oneThird / 254);
                green = 0;
                blue = 0;
            }
            image.set_pixel(
                x,
                y,
                red,
                green,
                blue
            );
            */

            pixel = map(x,y);
            image.set_pixel(
                x,
                y,
                pixel * 254,
                pixel * 254,
                pixel * 254
            );
        }
    }

    image.save_image("output.bmp");
}

#define randDouble() (double)random() / RAND_MAX
double randHeight(ParallelRNG& rng, double avg, double range, size_t subdivision) {
    #ifdef CRAPPY_GEN
    static const double conversion = 1.0 / 0xffffffff;
    static uint32_t rnd = 5;
    uint32_t x = rnd;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    rnd = x;
    return avg + x * conversion * (range / (1 << subdivision));
    #else
    return avg + rng.getDoublePlusMinus() * (range / (1 << subdivision) /*pow(2, subdivision)*/);
    #endif
}

void squareStep(Array2D<double>& map, ParallelRNG& rng, size_t row, size_t column, size_t sideLength, double range, size_t subdivision) {
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

void diamondStep (Array2D<double>& map, ParallelRNG& rng, size_t row, size_t column, size_t sideLength, size_t sideSquares, double range, size_t subdivision) {
    size_t halfEdge = sideLength / 2; /* sideLength >> 1 */
    size_t diamondX = column*sideLength + halfEdge;
    size_t diamondY = row*sideLength + halfEdge;

    size_t squareX, squareY;
    double sumAltitude;
    size_t numVerticies;
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
        squareX = diamondX + halfEdge*multiples[i][0];
        squareY = diamondY + halfEdge*multiples[i][1];

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
    rng.reinitialize();
    }
    for (size_t subdivision = 0; subdivision < subdivisions; subdivision++) {
        size_t sideSquares = pow(2,subdivision); /* 1 << n */
        size_t sideLength = (length - 1) / sideSquares; /* (length - 1) >> n */

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

double median(double data[], size_t size) {
    double swp;
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            if (data[i] > data[j]) {
                swp = data[i];
                data[i] = data[j];
                data[j] = swp;
            }
        }
    }
    return data[size / 2];
}

void printStats(double timings[], size_t numTimings) {
    std::cout << "Mean:\t" << gsl_stats_mean(timings, 1, numTimings) << std::endl;
    // std::cout << "Median:\t" << gsl_stats_median(timings, 1, numTimings) << std::endl;
    std::cout << "Standard Deviation:\t" << gsl_stats_sd(timings, 1, numTimings) << std::endl;
}

#define N 15
#define THREADS 8
#define NUM_TIMINGS 5
int main() {
    size_t length = pow(2, N) + 1;
    Array2D<double> map(length);
    map.set(0);

    double t1;
    double serialTimings[NUM_TIMINGS], parallelTimings[NUM_TIMINGS];

    //serial timings
    for (size_t i = 0; i < NUM_TIMINGS; i++) {
        t1 = omp_get_wtime();
        diamondSquare(map, length);
        serialTimings[i] = omp_get_wtime() - t1;
        std::cout << serialTimings[i] << std::endl;
        // printTerrain(map);
    }
    std::cout << "====SERIAL TIMINGS====" << std::endl;
    printStats(serialTimings, NUM_TIMINGS);
    std::cout << std::endl;
    printTerrain(map);

    //parallel timings
    for (size_t i = 0; i < NUM_TIMINGS; i++) {
        t1 = omp_get_wtime();
        diamondSquareParallel(map, length);
        parallelTimings[i] = omp_get_wtime() - t1;
        std::cout << parallelTimings[i] << std::endl;
        // printTerrain(map);
    }
    std::cout << "====PARALLEL TIMINGS====" << std::endl;
    printStats(parallelTimings, NUM_TIMINGS);
    std::cout << std::endl;

    double meanSerial = gsl_stats_mean(serialTimings, 1, NUM_TIMINGS);
    double meanParallel = gsl_stats_mean(parallelTimings, 1, NUM_TIMINGS);
    int maxThreads = omp_get_num_procs();

    std::cout << "====PARALLEL SPEEDUP====" << std::endl <<
    std::string("Total Speedup: ") + std::to_string(meanSerial / meanParallel) << std::endl <<
    std::string("Maximum Thread Count: ") + std::to_string(maxThreads) << std::endl <<
    std::string("Speedup per core: ") + std::to_string(meanSerial / (meanParallel * maxThreads)) << std::endl;

    printTerrain(map);

}
