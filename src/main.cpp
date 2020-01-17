#include <string.h>
#include <iostream>

#include <omp.h>
#include <gsl/gsl_statistics_double.h>

#include "DiamondSquare.hpp"
#include "CubeWorld.hpp"

struct options
{
    bool failure;
    bool readInput;
    char *inputFileName;
    size_t mapSize;
    bool writeOutput;
    char *outputFileName;
};

int argInputFile(int argIndex, char *argv[], struct options *opts) {
    opts->readInput = true;
    opts->inputFileName = argv[argIndex];
    return argIndex + 1;
}

int argMapSize(int argIndex, char *argv[], struct options *opts) {
    char *endptr;
    long size = strtol(argv[argIndex], &endptr, 0);
    opts->mapSize = size;
    if (endptr == argv[argIndex]) {
        opts->failure = true;
    }
    return argIndex + 1;
}

int argOutputFile(int argIndex, char *argv[], struct options *opts) {
    opts->writeOutput = true;
    opts->outputFileName = argv[argIndex];
    return argIndex + 1;
}

struct options parseArguments(int argc, char *argv[]) {
    int argIndex = 1;
    struct options opts = {
        .failure =          false,
        .readInput =        false,
        .inputFileName =    NULL,
        .mapSize =          5,
        .writeOutput =      false,
        .outputFileName =   NULL
    };
    while (argIndex < argc) {
        if (!strncmp(argv[argIndex], "-i", 3)) {
            argIndex = argInputFile(argIndex + 1, argv, &opts);
        } else if (!strncmp(argv[argIndex], "-s", 3)) {
            argIndex = argMapSize(argIndex + 1, argv, &opts);
        } else if (!strncmp(argv[argIndex], "-o", 3)) {
            argIndex = argOutputFile(argIndex + 1, argv, &opts);
        } else {
            opts.failure = true;
        }
        if (opts.failure) {
            argIndex = argc;
        }
    }
    return opts;
}

void printStats(double timings[], size_t numTimings) {
    std::cout << "Mean:\t" << gsl_stats_mean(timings, 1, numTimings) << std::endl;
    std::cout << "Standard Deviation:\t" << gsl_stats_sd(timings, 1, numTimings) << std::endl;
}

const char *help[] = {
    "Usage:",
    "\tparallel-world [-i <input_file>] [-s <map_size>] [-o <output_file>]",
    "\t",
    "\tMap size is the length of the side of a map. It must be 2^n + 1",
    "\t\twhere n is any positive integer."
};
const int helpLines = 5;
void printHelp() {
    for (size_t i = 0; i < helpLines; i++) {
        std::cout << help[i] << std::endl;
    }
}

#define N 9
#define NUM_TIMINGS 5
int main(int argc, char *argv[]) {
    struct options opts = parseArguments(argc, argv);

    if (opts.failure) {
        printHelp();
        return 0;
    }

    if (opts.readInput) {
        std::cout << "Input is not yet supported." << std::endl;
    }

    CubeWorld map(opts.mapSize);

    diamondSquareParallel(map);

    if (opts.writeOutput)
    {
        printTerrain(opts.outputFileName, map);
    }
    
    
    return 0;
}
