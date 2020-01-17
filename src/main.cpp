#include <string.h>
#include <iostream>

#include <omp.h>
#include <gsl/gsl_statistics_double.h>

#include "DiamondSquare.hpp"
#include "CubeWorld.hpp"

struct options
{
    bool failure;
    const char* failureMsg;
    bool readInput;
    const char *inputFileName;
    size_t mapSize;
    bool writeOutput;
    const char *outputFileName;
    bool benchmark;
    size_t benchmarkIntervals;
};

int argInputFile(int argIndex, int argc, char *argv[], struct options *opts) {
    opts->readInput = true;
    
    if (argIndex < argc) {
        opts->inputFileName = argv[argIndex];
    } else {
        opts->failure = true;
        opts->failureMsg = "You must specify a file name.";
    }

    return argIndex + 1;
}

int argMapSize(int argIndex, int argc, char *argv[], struct options *opts) {
    if (argIndex < argc) {
        char *endptr;
        long size = strtol(argv[argIndex], &endptr, 0);
        opts->mapSize = size;
        if (endptr == argv[argIndex]) {
            opts->failure = true;
            opts->failureMsg = "You must specify an integer 2^n + 1 for the map size";
        }
    } else {
        opts->failure = true;
        opts->failureMsg = "You must specify an integer 2^n + 1 for the map size";
    }
    return argIndex + 1;
}

int argOutputFile(int argIndex, int argc, char *argv[], struct options *opts) {
    opts->writeOutput = true;
    
    if (argIndex < argc) {
        opts->outputFileName = argv[argIndex];
    } else {
        opts->failure = true;
        opts->failureMsg = "You must specify a file name.";
    }
    return argIndex + 1;
}

int argBenchmark(int argIndex, int argc, char *argv[], struct options *opts) {
    opts->benchmark = true;

    if (argIndex < argc) {
        char *endptr;
        long intervals = strtol(argv[argIndex], &endptr, 0);
        if (endptr == argv[argIndex]) {
            opts->failure = true;
            opts->failureMsg = "You must specify the number of intervals";
        }
        opts->benchmarkIntervals = (int) intervals;
    } else {
        opts->failure = true;
        opts->failureMsg = "You must specify the number of intervals";
    }

    return argIndex + 1;
}

int argPower(int argIndex, int argc, char *argv[], struct options *opts) {
    if (argIndex < argc) {
        char *endptr;
        long size = strtol(argv[argIndex], &endptr, 0);
        opts->mapSize = (1 << size) + 1;
        if (endptr == argv[argIndex]) {
            opts->failure = true;
            opts->failureMsg = "You must specify an integer for the map power size";
        }
    } else {
        opts->failure = true;
        opts->failureMsg = "You must specify an integer for the map power size";
    }
    return argIndex + 1;
}

struct options parseArguments(int argc, char *argv[]) {
    int argIndex = 1;
    struct options opts = {
        .failure =              false,
        .readInput =            false,
        .inputFileName =        NULL,
        .mapSize =              5,
        .writeOutput =          false,
        .outputFileName =       NULL,
        .benchmark =            false,
        .benchmarkIntervals =   0
    };
    while (argIndex < argc) {
        if (!strncmp(argv[argIndex], "-i", 3))
        {
            argIndex = argInputFile(argIndex + 1, argc, argv, &opts);
        }
        else if (!strncmp(argv[argIndex], "-s", 3))
        {
            argIndex = argMapSize(argIndex + 1, argc, argv, &opts);
        }
        else if (!strncmp(argv[argIndex], "-o", 3))
        {
            argIndex = argOutputFile(argIndex + 1, argc, argv, &opts);
        }
        else if
            (!strncmp(argv[argIndex], "-b", 3) ||
             !strncmp(argv[argIndex], "--benchmark", 11))
        {
            argIndex = argBenchmark(argIndex + 1, argc, argv, &opts);
        }
        else if (!strncmp(argv[argIndex], "--power", 7))
        {
            argIndex = argPower(argIndex + 1, argc, argv, &opts);
        }
        else
        {
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
    "\tparallel-world [-i <input_file>] [-s <map_size> | -p <map_power>] [-o <output_file>]",
    "\t\t[-b|--benchmark <intervals>]",
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

int main(int argc, char *argv[]) {
    struct options opts = parseArguments(argc, argv);

    if (opts.failure) {
        if (opts.failureMsg != NULL) {
            std::cout << opts.failureMsg << "\n" << std::endl;
        }
        printHelp();
        return 0;
    }

    if (opts.readInput) {
        std::cout << "Input is not yet supported." << std::endl;
    }

    CubeWorld map(opts.mapSize);

    if (opts.benchmark) {
        double   serialTimes[opts.benchmarkIntervals];
        double parallelTimes[opts.benchmarkIntervals];

        for (size_t i = 0; i < opts.benchmarkIntervals; i++) {
            double t1 = omp_get_wtime();
            diamondSquare(map);
            std::cout
                << "Generation Time: "
                << (serialTimes[i] = omp_get_wtime() - t1)
                << "s"
                << std::endl
            ;
        }
        std::cout << "\nSERIAL TIMINGS" << std::endl;
        printStats(serialTimes, opts.benchmarkIntervals);
        std::cout << std::endl;

        for (size_t i = 0; i < opts.benchmarkIntervals; i++) {
            double t1 = omp_get_wtime();
            diamondSquareParallel(map);
            std::cout
                << "Generation Time: "
                << (parallelTimes[i] = omp_get_wtime() - t1)
                << "s"
                << std::endl
            ;
        }
        std::cout << "\nPARALLEL TIMINGS" << std::endl;
        printStats(parallelTimes, opts.benchmarkIntervals);
        std::cout << std::endl;
    } else {
        diamondSquareParallel(map);
    }

    if (opts.writeOutput)
    {
        printTerrain(opts.outputFileName, map);
    }
    
    
    return 0;
}
