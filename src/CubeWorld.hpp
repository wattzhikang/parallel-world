#ifndef CUBE_WORLD
#define CUBE_WORLD

#include <gsl/gsl_linalg.h>

#include "../lib/bitmap/bitmap_image.hpp"

struct transform {
    char face;
    double matrix[9];
};

class CubeWorld {
    private:
        struct transform transformMap[6][4][4];
        gsl_matrix *data[6];
        size_t size;

        void initializeMatricies();
        long getDisplacedFaces(long position);
        void navigate(char *face, long *x, long *y);
    public:
        CubeWorld(size_t size);
        ~CubeWorld();
        size_t getSize();
        double get(char face, long x, long y);
        void set(char face, long x, long y, double value);
        void setAll(double value);
};

void printTerrain(std::string title, CubeWorld& map);
void printTerrain(CubeWorld& map);

#endif