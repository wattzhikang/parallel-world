#ifndef PARALLEL_RNG
#define PARALLEL_RNG

#include <gsl/gsl_rng.h>

/*
Because every random number generator requires its own state, when using the
standard library random number generator in parallel, all threads have to block
when using it. This results in extremely poor performance. This class solves
that problem by creating a GNU Scientific Library random number generator for
each thread. Moreover, it is completely transparent---simply initialize it in
a parallel region and it will allocate a RNG for every thread, and every call
will call the appropriate RNG for that thread.
*/
class ParallelRNG {
    private:
        gsl_rng **generators = NULL;
        size_t quantity = 0;
        bool activated = false;
    public:
        ~ParallelRNG() {
            if (activated) {
                deinitialize();
            }
        }
        void deinitialize() {
            for (size_t i = 0; i < quantity; i++) {
                gsl_rng_free(generators[i]);
            }
            delete[] generators;
            
            activated = false;
        }
        void reinitialize() {
            if (activated) {
                deinitialize();
            }

            quantity = omp_get_num_threads();

            generators = new gsl_rng*[quantity];

            gsl_rng* seeder = gsl_rng_alloc(gsl_rng_taus2);

            for (size_t i = 0; i < quantity; i++) {
                generators[i] = gsl_rng_alloc(gsl_rng_taus2);
                gsl_rng_set(generators[i], gsl_rng_get(seeder));
            }

            gsl_rng_free(seeder);

            activated = true;
        }
        double getDouble() {
            if (activated) {
                return gsl_rng_uniform(generators[omp_get_thread_num()]);
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

#endif