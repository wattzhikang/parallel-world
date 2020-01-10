#include <string.h>
#include <iostream>

#include <omp.h>
#include <gsl/gsl_statistics_double.h>

#include "DiamondSquare.hpp"
#include "CubeWorld.hpp"

void printStats(double timings[], size_t numTimings) {
    std::cout << "Mean:\t" << gsl_stats_mean(timings, 1, numTimings) << std::endl;
    std::cout << "Standard Deviation:\t" << gsl_stats_sd(timings, 1, numTimings) << std::endl;
}

#define N 10
#define NUM_TIMINGS 5
int main() {
    size_t length = pow(2, N) + 1;
    CubeWorld map(length);

    double t1;
    double serialTimings[NUM_TIMINGS]; //, parallelTimings[NUM_TIMINGS];

    //serial timings
    for (size_t i = 0; i < NUM_TIMINGS; i++) {
        map.setAll(0.1);
        t1 = omp_get_wtime();
        diamondSquare(map);
        serialTimings[i] = omp_get_wtime() - t1;
        std::cout << serialTimings[i] << std::endl;
    }
    std::cout << "====SERIAL TIMINGS====" << std::endl;
    printStats(serialTimings, NUM_TIMINGS);
    std::cout << std::endl;
    printTerrain(map);

    //parallel timings
    // for (size_t i = 0; i < NUM_TIMINGS; i++) {
    //     t1 = omp_get_wtime();
    //     diamondSquareParallel(map, length);
    //     parallelTimings[i] = omp_get_wtime() - t1;
    //     std::cout << parallelTimings[i] << std::endl;
    // }
    // std::cout << "====PARALLEL TIMINGS====" << std::endl;
    // printStats(parallelTimings, NUM_TIMINGS);
    // std::cout << std::endl;

    // Print speedup statistics

    // double meanSerial = gsl_stats_mean(serialTimings, 1, NUM_TIMINGS);
    // double meanParallel = gsl_stats_mean(parallelTimings, 1, NUM_TIMINGS);
    // int maxThreads = omp_get_num_procs();

    // std::cout << "====PARALLEL SPEEDUP====" << std::endl <<
    // std::string("Total Speedup: ") + std::to_string(meanSerial / meanParallel) << std::endl <<
    // std::string("Maximum Thread Count: ") + std::to_string(maxThreads) << std::endl <<
    // std::string("Speedup per core: ") + std::to_string(meanSerial / (meanParallel * maxThreads)) << std::endl;

    // printTerrain(map);

}
