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
        gsl_matrix *data[6];

        size_t size;

        long getDisplacedFaces(long position);
        char displacementIndex(long facesDisplaced);

        #define COS0     1.0
        #define COSPID2  0.0
        #define COSPI    -1.0
        #define COS3PID2 0.0
        #define SIN0     0.0
        #define SINPID2  1.0
        #define SINPI    0.0
        #define SIN3PID2 -1.0
        struct transform face4TMosaic[4] = {
            {4, {COS0,     -(SIN0),     0.0, SIN0,     COS0,     0.0, 0.0, 0.0, 1.0}}, // face 4 rotation 0
            {1, {COS3PID2, -(SIN3PID2), 0.0, SIN3PID2, COS3PID2, 0.0, 0.0, 0.0, 1.0}}, // face 1 rotation 3pi/2
            {5, {COSPI,    -(SINPI),    0.0, SINPI,    COSPI,    0.0, 0.0, 0.0, 1.0}}, // face 5 rotation pi
            {3, {COSPID2,  -(SINPID2),  0.0, SINPID2,  COSPID2,  0.0, 0.0, 0.0, 1.0}}  // face 3 rotation pi/2
        };
        struct transform face5TMosaic[4] = {
            {5, {COS0,     -(SIN0),     0.0, SIN0,     COS0,     0.0, 0.0, 0.0, 1.0}}, // face 5 rotation 0
            {1, {COSPID2,  -(SINPID2),  0.0, SINPID2,  COSPID2,  0.0, 0.0, 0.0, 1.0}}, // face 1 rotation pi/2
            {4, {COSPI,    -(SINPI),    0.0, SINPI,    COSPI,    0.0, 0.0, 0.0, 1.0}}, // face 4 rotation pi
            {3, {COS3PID2, -(SIN3PID2), 0.0, SIN3PID2, COS3PID2, 0.0, 0.0, 0.0, 1.0}} // face 3 rotation 3pi/2
        };
        struct transform horizontalTMosaic[4][4] = {
            { // up from face 0
                {0, {COS0,  -(SIN0),  0.0, SIN0,  COS0,  0.0, 0.0, 0.0, 1.0}}, // face 0 rotation 0
                {4, {COS0,  -(SIN0),  0.0, SIN0,  COS0,  0.0, 0.0, 0.0, 1.0}}, // face 4 rotation 0
                {2, {COSPI, -(SINPI), 0.0, SINPI, COSPI, 0.0, 0.0, 0.0, 1.0}}, // face 2 rotation pi
                {5, {COS0,  -(SIN0),  0.0, SIN0,  COS0,  0.0, 0.0, 0.0, 1.0}}  // face 5 rotation 0
            },
            { // up from face 1
                {1, {COS0,     -(SIN0),     0.0, SIN0,     COS0,     0.0, 0.0, 0.0, 1.0}}, // face 1 rotation 0
                {4, {COSPID2,  -(SINPID2),  0.0, SINPID2,  COSPID2,  0.0, 0.0, 0.0, 1.0}}, // face 4 rotation pi/2
                {3, {COSPI,    -(SINPI),    0.0, SINPI,    COSPI,    0.0, 0.0, 0.0, 1.0}}, // face 3 rotation pi
                {5, {COS3PID2, -(SIN3PID2), 0.0, SIN3PID2, COS3PID2, 0.0, 0.0, 0.0, 1.0}}  // face 5 rotation 3pi/2
            },
            { // up from face 2
                {2, {COS0,  -(SIN0),  0.0, SIN0,  COS0,  0.0, 0.0, 0.0, 1.0}}, // face 2 rotation 0
                {4, {COSPI, -(SINPI), 0.0, SINPI, COSPI, 0.0, 0.0, 0.0, 1.0}}, // face 4 rotation pi
                {0, {COSPI, -(SINPI), 0.0, SINPI, COSPI, 0.0, 0.0, 0.0, 1.0}}, // face 0 rotation pi
                {5, {COSPI, -(SINPI), 0.0, SINPI, COSPI, 0.0, 0.0, 0.0, 1.0}}  // face 5 rotation pi
            },
            { // up from face 3
                {3, {COS0,     -(SIN0),     0.0, SIN0,     COS0,     0.0, 0.0, 0.0, 1.0}}, // face 3 rotation 0
                {4, {COS3PID2, -(SIN3PID2), 0.0, SIN3PID2, COS3PID2, 0.0, 0.0, 0.0, 1.0}}, // face 4 rotation 3pi/2
                {1, {COSPI,    -(SINPI),    0.0, SINPI,    COSPI,    0.0, 0.0, 0.0, 1.0}}, // face 1 rotation pi
                {5, {COSPID2,  -(SINPID2),  0.0, SINPID2,  COSPID2,  0.0, 0.0, 0.0, 1.0}}  // face 5 rotation pi/2
            }
        };
        double translate[9] = {1.0, 0.0, -((double)(size / 2)), 0.0, 1.0, -((double)(size / 2)), 0.0, 0.0, 1.0};
        double backTranslate[9] = {1.0, 0.0, ((double)(size / 2)), 0.0, 1.0, ((double)(size / 2)), 0.0, 0.0, 1.0};

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