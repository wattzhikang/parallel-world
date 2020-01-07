#ifndef CUBE_WORLD
#define CUBE_WORLD

#include <gsl/gsl_linalg.h>

#include "../lib/bitmap/bitmap_image.hpp"

class CubeWorld {
    private:
        gsl_matrix *data[6];

        size_t size;
        
        char horizontalNavigation[6][4] = {
            {0, 1, 2, 3}, // Face 0
            {1, 2, 3, 0}, // Face 1
            {2, 3, 0, 1}, // Face 2
            {3, 0, 1, 2}, // Face 3
            {4, 1, 5, 3}, // Face 4
            {5, 1, 4, 3}  // Face 5
        };
        char verticalNavigation[6][4] = {
            {0, 4, 2, 5}, // Face 0
            {1, 4, 3, 5}, // Face 1
            {2, 4, 0, 5}, // Face 2
            {3, 4, 1, 5}, // Face 3
            {4, 2, 5, 0}, // Face 4
            {5, 0, 4, 2}  // Face 5
        };
        #define SIN 1
        #define COS 0
        char transformMosaic[4][4][2] = {
            { // up from face 0
                {1, 0}, // face 0 rotation 0
                {1, 0}, // face 4 rotation 0
                {-1,0}, // face 2 rotation pi
                {1, 0}  //face 5 rotation 0
            },
            { // up from face 1
                {1, 0}, // face 1 rotation 0
                {0, 1}, // face 4 rotation pi/2
                {-1,0}, // face 3 rotation pi
                {0,-1}  // face 5 rotation 3pi/2
            },
            { // up from face 2
                {1, 0}, // face 2 rotation 0
                {-1,0}, // face 4 rotation pi
                {-1,0}, // face 0 rotation pi
                {-1,0}  // face 5 rotation pi
            },
            { // up from face 3
                {1, 0}, // face 3 rotation 0
                {0,-1}, // face 4 rotation 3pi/2
                {-1,0}, // face 1 rotation pi
                {0, 1}  // face 5 rotation pi/2
            }
        };
        long getDisplacedFaces(long position) {
            if (position < 0) {
                return ((position + 1) / ((long)size)) - 1;
            }
            if (position >= (long)size) {
                return position / ((long)size);
            }
            return 0;
        }
        char displacementIndex(long facesDisplaced) {
            if (facesDisplaced < 0) {
                return (facesDisplaced % 4) + 4;
            } else {
                return facesDisplaced % 4;
            }
        }
        struct transform {
            char face;
            double matrix[9];
        };
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
        void navigate(char *face, long *x, long *y) {
            double xy[3] = {
                *x - getDisplacedFaces(*x) * ((long)size),
                *y - getDisplacedFaces(*y) * ((long)size),
                1.0
            };
            double result[3];
            double *transformations[4];
            char multiplication = 0;
            char multiplications = 3;

            transformations[multiplication] = translate;
            multiplication++;

            long hIndex, vIndex;

            // first handle horizontal face displacement, including the special cases
            // of faces 4 and 5

            switch (*face)
            {
            case 5:
                hIndex = getDisplacedFaces(*x) % 4;
                if (hIndex < 0) {
                    hIndex += 4;
                }
                transformations[multiplication] = face5TMosaic[hIndex].matrix;
                multiplication++;
                multiplications = 4;
                switch (face5TMosaic[hIndex].face) { //set hIndex and vIndex for the vertical displacement phase
                    case 5:
                        hIndex = 0;
                        vIndex = 3;
                        break;
                    case 4:
                        hIndex = 0;
                        vIndex = 1;
                        break;
                    default:
                        hIndex = face5TMosaic[hIndex].face;
                        vIndex = 0;
                        break;
                }
                break;
            case 4:
                hIndex = getDisplacedFaces(*x) % 4;
                if (hIndex < 0) {
                    hIndex += 4;
                }
                transformations[multiplication] = face4TMosaic[hIndex].matrix;
                multiplication++;
                multiplications = 4;
                switch (face4TMosaic[hIndex].face) {
                    case 5:
                        hIndex = 0;
                        vIndex = 3;
                        break;
                    case 4:
                        hIndex = 0;
                        vIndex = 1;
                        break;
                    default:
                        hIndex = face4TMosaic[hIndex].face;
                        vIndex = 0;
                        break;
                }
                break;
            default:
                //hIndex should be the resulting face
                hIndex = ( *face + getDisplacedFaces(*x)) % 4;
                if (hIndex < 0) {
                    hIndex += 4;
                }
                vIndex = 0;
                break;
            }

            // then handle vertical face displacement

            vIndex += getDisplacedFaces(*y);
            vIndex %= 4;
            if (vIndex < 0) {
                vIndex += 4;
            }

            transformations[multiplication] = horizontalTMosaic[hIndex][vIndex].matrix;
            multiplication++;

            transformations[multiplication] = backTranslate;
            multiplication++;

            for (multiplication = 0; multiplication < multiplications; multiplication++) {
                cblas_dgemv(
                    CblasRowMajor,
                    CblasNoTrans,
                    3, 3,
                    1.0,
                    transformations[multiplication],
                    3, 
                    xy,
                    1,
                    0.0,
                    result,
                    1
                );
                memcpy(xy, result, 3 * sizeof(double));
            }

            *face = horizontalTMosaic[hIndex][vIndex].face;
            *x = (long) result[0];
            *y = (long) result[1];
        }
        // void navigate(char *face, long *x, long *y) {
        //     char newFace;
        //     long newX, newY;

        //     long horizontalFacesDisplaced = getDisplacedFaces(*x);
        //     newFace = horizontalNavigation[*face][displacementIndex(horizontalFacesDisplaced)];

        //     long verticalFacesDisplaced = getDisplacedFaces(*y);
        //     newFace = verticalNavigation[newFace][displacementIndex(verticalFacesDisplaced)];

        //     newX = *x - horizontalFacesDisplaced * ((long)size);
        //     newY = *y - verticalFacesDisplaced * ((long)size);

        //     //rotate coordinates
        //     rotate(*face, horizontalFacesDisplaced, verticalFacesDisplaced, &newX, &newY);

        //     //output
        //     *face = newFace;
        //     *x = newX;
        //     *y = newY;
        // }
    public:
        CubeWorld(size_t size) : size(size) {
            for (size_t i = 0; i < 6; i++) {
                data[i] = gsl_matrix_alloc(size, size);
            }
        }
        ~CubeWorld() {
            for (size_t i = 0; i < 6; i++) {
                gsl_matrix_free(data[i]);
            }
        }
        size_t getSize() {
            return size;
        }
        double get(char face, long x, long y) {
            navigate(&face, &x, &y);
            return gsl_matrix_get(data[face], x, y);
        }
        void set(char face, long x, long y, double value) {
            navigate(&face, &x, &y);
            gsl_matrix_set(data[face], x, y, value);
        }
        void setAll(double value) {
            for (char face = 0; face < 6; face++) {
                for (size_t x = 0; x < size; x++) {
                    for (size_t y = 0; y < size; y++) {
                        gsl_matrix_set(data[face], x, y, value);
                    }
                }
            }
        }
};

static void printTerrain(std::string title, CubeWorld& map) {
    size_t length = map.getSize();
    bitmap_image image(map.getSize() * 4, map.getSize() * 3);

    #define X 0
    #define Y 1
    size_t faceOffsets[6][2] = {
        {map.getSize(), map.getSize()+1},      // face 0
        {2 * map.getSize(), map.getSize()+1},  // face 1
        {3 * map.getSize(), map.getSize()+1},  // face 2
        {0, map.getSize()+1},                  // face 3
        {map.getSize(),2*map.getSize()+1},     // face 4
        {map.getSize(),0+1}                    // face 5
    };
    size_t x, y;
    double pixel;

    for (char face = 0; face < 6; face++) {
        for (x = 0; x < length; x++) {
            for (y = 0; y < length; y++) {
                pixel = map.get(face, x, y);
                image.set_pixel(
                    x + faceOffsets[face][X],
                    map.getSize() * 3 - (y + faceOffsets[face][Y]),
                    pixel * 254,
                    pixel * 254,
                    pixel * 254
                );
            }
        }
    }

    image.save_image(title);
};

static void printTerrain(CubeWorld& map) {
    printTerrain(std::string("output.bmp"), map);
};

#endif