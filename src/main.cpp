#include <string.h>
#include <iostream>

#include <omp.h>
#include <gsl/gsl_statistics_double.h>

#include "../lib/bitmap/bitmap_image.hpp"

#include "DiamondSquare.hpp"

#define u_int32_t_MAX 0xFFFFFFFF
// Outputs the terrain as a grayscale bitmap
void printTerrain(Array2D<double>& map) {
    size_t length = map.getSize();
    bitmap_image image(length, map.getSize());

    u_int32_t red, green, blue;

    u_int32_t oneThird = (u_int32_t) (((size_t) u_int32_t_MAX) / 3);
    u_int32_t twoThirds = (u_int32_t) ( 2 * ((size_t) u_int32_t_MAX) / 3);

    size_t x, y;
    double pixel;

    for (x = 0; x < length; x++) {
        for (y = 0; y < length; y++) {
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

void printStats(double timings[], size_t numTimings) {
    std::cout << "Mean:\t" << gsl_stats_mean(timings, 1, numTimings) << std::endl;
    std::cout << "Standard Deviation:\t" << gsl_stats_sd(timings, 1, numTimings) << std::endl;
}

#define N 9
#define NUM_TIMINGS 5
int main() {
    size_t length = pow(2, N) + 1;
    Array2D<double> map(length);
    map.set(0); // sets all cells in the map to 0

    double t1;
    double serialTimings[NUM_TIMINGS], parallelTimings[NUM_TIMINGS];

    //serial timings
    for (size_t i = 0; i < NUM_TIMINGS; i++) {
        t1 = omp_get_wtime();
        diamondSquare(map, length);
        serialTimings[i] = omp_get_wtime() - t1;
        std::cout << serialTimings[i] << std::endl;
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
    }
    std::cout << "====PARALLEL TIMINGS====" << std::endl;
    printStats(parallelTimings, NUM_TIMINGS);
    std::cout << std::endl;

    // Print speedup statistics

    double meanSerial = gsl_stats_mean(serialTimings, 1, NUM_TIMINGS);
    double meanParallel = gsl_stats_mean(parallelTimings, 1, NUM_TIMINGS);
    int maxThreads = omp_get_num_procs();

    std::cout << "====PARALLEL SPEEDUP====" << std::endl <<
    std::string("Total Speedup: ") + std::to_string(meanSerial / meanParallel) << std::endl <<
    std::string("Maximum Thread Count: ") + std::to_string(maxThreads) << std::endl <<
    std::string("Speedup per core: ") + std::to_string(meanSerial / (meanParallel * maxThreads)) << std::endl;

    printTerrain(map);

}
