/// @file

#include <string.h>
#include <iostream>

#include <omp.h>
#include <gsl/gsl_statistics_double.h>

#include "DiamondSquare.hpp"
#include "CubeWorld.hpp"

/** \brief Stores all information parsed from the command line
 * 
*/
struct options
{
    bool failure; ///< Indicates whether the command line arguments were successfully parsed or not
    const char* failureMsg; ///< Describes the nature of the error encountered, if any, in parsing the command line arguments
    bool readInput; ///< True if the user specified an input file
    const char *inputFileName; ///< The name of the input file
    size_t mapSize; ///< Size of the map. Defaults to 5
    bool writeOutput; ///< True if the user specified an output file
    const char *outputFileName; ///< Name of the output file
    bool benchmark; ///< True if the user asked for benchmarks
    size_t benchmarkIntervals; ///< The number of benchmark intervals the user asked for
    bool help; ///< True if the user asked for the command line options
    unsigned long seed; ///< The random seed. Defaults to 0
    bool writeRectangle; ///< True if the user asked to write an equirectangular projection
    const char *outputEquirectangularFileName; ///< Name of the output file (equirectangular)
};

/**
 * Gets the name of the input file
 * @param argIndex The index of the next token after the command switch
 * @param argc The total number of tokens pulled from the command line
 * @param argv The raw tokens pulled from the command line
 * @param opts Stores the information parsed in this function
*/ 
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

/** \brief Reads a long from a string
 *
*/
bool readLong(char *str, long *destination) { //TODO: use this for map size and other
    char *endptr;                             //functions that need numbers
    long size = strtol(str, &endptr, 0);
    *destination = size;
    if (endptr == str) {
        return false;
    }
    return true;
}

/** \brief Gets the map size that the user specified
 * 
*/ 
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

/** \brief Gets the name of the user-specified output file
 * 
*/
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

/** \brief Gets the name of the user-specified equirectangular output file
 * 
*/
int argOutputEquirectangular(int argIndex, int argc, char *argv[], struct options *opts) {
    opts->writeRectangle = true;
    
    if (argIndex < argc) {
        opts->outputEquirectangularFileName = argv[argIndex];
    } else {
        opts->failure = true;
        opts->failureMsg = "You must specify a file name.";
    }
    return argIndex + 1;
}

/** \brief Gets the number of benchmark intervals the user specified
 * 
*/
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

/** \brief Get the map size as the power of 2 that the user specified
 * 
*/
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

/** \brief Gets the specified random seed
 * 
*/
int argSeed(int argIndex, int argc, char *argv[], struct options *opts) {
    long seed;
    if (argIndex < argc && readLong(argv[argIndex], &seed)) {
        opts->seed = seed;
    } else {
        opts->failure = true;
        opts->failureMsg = "You must input a valid number for a random seed";
    }
    return argIndex + 1;
}

/** \brief Indicates that the user has requested the help options
 * 
*/ 
int argHelp(int argIndex, int argc, char *argv[], struct options *opts) {
    opts->help = true;
    return argIndex;
}

/** \brief Parses options from the command line
 * Finds command switches and calls the appropriate parser methods to handle them
 * @param argc The number of tokens from the command line
 * @param argv The raw tokens from the command line
*/
struct options parseArguments(int argc, char *argv[]) {
    int argIndex = 1;

    //sets the defaults
    struct options opts = {
        .failure =              false,
        .readInput =            false,
        .inputFileName =        NULL,
        .mapSize =              5,
        .writeOutput =          false,
        .outputFileName =       NULL,
        .benchmark =            false,
        .benchmarkIntervals =   0,
        .help =                 false,
        .seed =                 0l
    };

    // Runs through the command line tokens, finds the command switches, and
    // calls the appropriate methods to parse user specifications, if applicable.
    while (argIndex < argc) {
        if (!strncmp(argv[argIndex], "-i", 3)) // The ! is necessary because 0 is false in a conditional
        {
            // The output from all of these functions is the index of the next switch, after the user
            // specifications, if any. The point is to set the index so that the loop is ready to parse
            // the next swtich.
            argIndex = argInputFile(argIndex + 1, argc, argv, &opts);
        }
        if (!strncmp(argv[argIndex], "-s", 3)) {
            argIndex = argSeed(argIndex + 1, argc, argv, &opts);
        }
        else if (!strncmp(argv[argIndex], "-r", 3))
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
        else if (!strncmp(argv[argIndex], "-h", 2))
        {
            argIndex = argHelp(argIndex + 1, argc, argv, &opts);
        }
        else if (!strncmp(argv[argIndex], "--power", 7))
        {
            argIndex = argPower(argIndex + 1, argc, argv, &opts);
        }
        else if (!strncmp(argv[argIndex], "-e", 3)) {
            argIndex = argOutputEquirectangular(argIndex + 1, argc, argv, &opts);
        }
        else
        {
            opts.failure = true;
        }

        // If any failure occurs in parsing, the argIndex will be set such that the loop will
        // immediately exit.
        if (opts.failure) {
            argIndex = argc;
        }
    }
    return opts;
}

/** \brief When benchmarks are run, this function will print the mean and standard deviations of the interval times
 * @param timings An array of the duration of each interval
 * @param numTimings The number of elements in timings
*/
void printStats(double timings[], size_t numTimings) {
    std::cout << "Mean:\t" << gsl_stats_mean(timings, 1, numTimings) << std::endl;
    std::cout << "Standard Deviation:\t" << gsl_stats_sd(timings, 1, numTimings) << std::endl;
}

/** \brief Each line of the usage information for this program
 * 
*/
const char *help[] = {
    "Usage:",
    "\tparallel-world [-i <input_file> | -s <random_seed>] [-r <map_resolution> | --power <map_power>] [-o <output_file>]",
    "\t\t[-b|--benchmark <intervals>] [-e <output_file>] [ -h ]",
    "\t",
    "\tMap resolution is the length of the side of a map. It must be 2^n + 1",
    "\t\twhere n is any positive integer."
};
const int helpLines = 5; ///< Number of lines in the usage

/*** \brief Print the usage information
 * 
*/
void printHelp() {
    for (size_t i = 0; i < helpLines; i++) {
        std::cout << help[i] << std::endl;
    }
}

int main(int argc, char *argv[]) {
    struct options opts = parseArguments(argc, argv);

    if (opts.failure || opts.help) {
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
        // The time durations of the serial intervals and parallel intervals
        double   serialTimes[opts.benchmarkIntervals];
        double parallelTimes[opts.benchmarkIntervals];

        // Run the serial intervals
        for (size_t i = 0; i < opts.benchmarkIntervals; i++) {
            double t1 = omp_get_wtime();

            diamondSquare(map, opts.seed);

            std::cout
                << "Generation Time: "
                << (serialTimes[i] = omp_get_wtime() - t1)
                << "s"
                << std::endl
            ;
            serialTimes[i] = omp_get_wtime() - t1;
        }
        std::cout << "\nSERIAL TIMINGS" << std::endl;
        printStats(serialTimes, opts.benchmarkIntervals);
        std::cout << std::endl;

        for (size_t i = 0; i < opts.benchmarkIntervals; i++) {
            double t1 = omp_get_wtime();

            diamondSquareParallel(map, opts.seed);
            
            std::cout
                << "Generation Time: "
                << (parallelTimes[i] = omp_get_wtime() - t1)
                << "s"
                << std::endl
            ;
            serialTimes[i] = omp_get_wtime() - t1;
        }
        std::cout << "\nPARALLEL TIMINGS" << std::endl;
        printStats(parallelTimes, opts.benchmarkIntervals);
        std::cout << std::endl;

    } else {
        diamondSquareParallel(map, opts.seed);
    }

    if (opts.writeOutput)
    {
        printTerrain(opts.outputFileName, map);
    }

    if (opts.writeRectangle) {
        printEquirectangular(opts.outputEquirectangularFileName, map);
    }

    return 0;
}
