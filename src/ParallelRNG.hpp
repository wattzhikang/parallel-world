#ifndef PARALLEL_RNG
#define PARALLEL_RNG

#include <gsl/gsl_rng.h>

/**
 * Provides access to a single random number generator
*/
class ParallelRNGSequence {
    private:
        gsl_rng *generator;
    public:
        /**
         * This method is not meant to be used by anything other than a ParallelRNG instance, though
         * it would cause no errors otherwise. This method constructs a Random Number Generator with
         * a specified seed.
         * @param seed The seed
        */
        ParallelRNGSequence(unsigned long seed) : generator(generator) {
            generator = gsl_rng_alloc(gsl_rng_taus2);
            gsl_rng_set(generator, seed);
        };
        /**
         * Cleanly destroys an instance of this class
        */
        ~ParallelRNGSequence() {
            gsl_rng_free(generator);
        }
        /**
         * @return a double between 0 and 1
        */
        float getFloat() {
            return (float) gsl_rng_uniform(generator);
        }
        /**
         * @return a double between -1 and 1
        */
        float getFloatPlusMinus() {
            return getFloat() * 2.0 - 1.0;
        }
        /**
         * Provides raw access to the GNU Scientific Library random number generator struct
         * @return the underlying random number generator
        */
        gsl_rng *getRNG() {
            return generator;
        }
};

/**
Random number generators require state, which means that they cannot be used in parallel.
This class facilitates parallel optimization by instantiating an arbitrary number of
random number generators. To execute in parallel, simply give each thread a different
random number generator. For consistency between serial and parallel areas, it is
important that these random number generators are used the same way in both---even
though they do not need to be.
*/
class ParallelRNG {
    private:
        ParallelRNGSequence **sequences;
        size_t quantity = 0;
    public:
        /**
         * Creates a number of random number generators.
         * @param numSequences The number of random number generators to create
        */
        ParallelRNG(size_t numSequences, unsigned long int seed) : quantity(numSequences) {
            sequences = new ParallelRNGSequence*[quantity];

            gsl_rng* seeder = gsl_rng_alloc(gsl_rng_taus2);
            gsl_rng_set(seeder, seed);

            for (size_t i = 0; i < quantity; i++) {
                sequences[i] = new ParallelRNGSequence(gsl_rng_get(seeder));
            }

            gsl_rng_free(seeder);
        }
        /**
         * Cleanly destroys all the random number generators
        */
        ~ParallelRNG() {
            for (size_t i = 0; i < quantity; i++) {
                delete sequences[i];
            }
            delete[] sequences;
        }
        /**
         * Gets a particular random number generator for use
         * @param index Specifies which random number generator you want. This can be used to ensure consistency.
         * @return A pointer to a ParallelRNGSequence, by which you can get random numbers
        */
        ParallelRNGSequence* getSequence(size_t index) {
            return sequences[index];
        }
        /**
         * @return The number of random number generators that can be used in parallel
        */
        size_t getNumGenerators() {
            return quantity;
        }
};

#endif