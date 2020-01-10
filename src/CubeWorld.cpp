#include "CubeWorld.hpp"

void CubeWorld::initializeMatricies() {
    //                COS      SIN    0 ;  COS   SIN   0 ;  0    0    1
    #define ROT0     {1.0,  -( 0.0), 0.0,  0.0,  1.0, 0.0, 0.0, 0.0, 1.0}
    #define ROTPID2  {0.0,  -( 1.0), 0.0,  1.0,  0.0, 0.0, 0.0, 0.0, 1.0}
    #define ROTPI    {-1.0, -( 0.0), 0.0,  0.0, -1.0, 0.0, 0.0, 0.0, 1.0}
    #define ROT3PID2 {0.0,  -(-1.0), 0.0, -1.0,  0.0, 0.0, 0.0, 0.0, 1.0}
    struct transform horizontalTransforms[6][4] =
    {
        { // right from face 0
            {0, ROT0},
            {1, ROT0},
            {2, ROT0},
            {3, ROT0}
        },
        { // right from face 1
            {1, ROT0},
            {2, ROT0},
            {3, ROT0},
            {0, ROT0}
        },
        { // right from face 2
            {2, ROT0},
            {3, ROT0},
            {0, ROT0},
            {1, ROT0}
        },
        { // right from face 3
            {3, ROT0},
            {0, ROT0},
            {1, ROT0},
            {2, ROT0}
        },
        { // right from face 4
            {4, ROT0},
            {1, ROT3PID2},
            {5, ROTPI},
            {3, ROTPID2}
        },
        { // right from face 5
            {5, ROT0},
            {1, ROTPID2},
            {4, ROTPI},
            {3, ROT3PID2}
        }
    };
    struct transform verticalTransforms[6][4] =
    {
        { // up from face 0
            {0, ROT0},
            {4, ROT0},
            {2, ROTPI},
            {5, ROT0}
        },
        { // up from face 1
            {1, ROT0},
            {4, ROTPID2},
            {3, ROTPI},
            {5, ROT3PID2}
        },
        { // up from face 2
            {2, ROT0},
            {4, ROTPI},
            {0, ROTPI},
            {5, ROTPI}
        },
        { // up from face 3
            {3, ROT0},
            {4, ROT3PID2},
            {1, ROTPI},
            {5, ROTPID2}
        },
        { // up from face 4
            {4, ROT0},
            {2, ROTPI},
            {5, ROT0},
            {0, ROT0}
        },
        { // up from face 5
            {5, ROT0},
            {0, ROT0},
            {4, ROT0},
            {2, ROTPI}
        }
    };

    double translate[9] = {1.0, 0.0, -((double)(size / 2)), 0.0, 1.0, -((double)(size / 2)), 0.0, 0.0, 1.0};
    double backTranslate[9] = {1.0, 0.0, ((double)(size / 2)), 0.0, 1.0, ((double)(size / 2)), 0.0, 0.0, 1.0};

    double C[9];

    for (size_t face = 0; face < 6; face++) {
        for (size_t hIndex = 0; hIndex < 4; hIndex++) {
            for (size_t vIndex = 0; vIndex < 4; vIndex++)
            {
                cblas_dgemm(
                    CblasRowMajor,
                    CblasNoTrans,
                    CblasNoTrans,
                    3, 3,
                    3,
                    1.0,
                    horizontalTransforms[face][hIndex].matrix,
                    3,
                    translate,
                    3,
                    0.0,
                    C,
                    3
                );
                memcpy(transformMap[face][hIndex][vIndex].matrix, C, 9 * sizeof(double));
                cblas_dgemm(
                    CblasRowMajor,
                    CblasNoTrans,
                    CblasNoTrans,
                    3, 3,
                    3,
                    1.0,
                    verticalTransforms[face][vIndex].matrix,
                    3,
                    transformMap[face][hIndex][vIndex].matrix,
                    3,
                    0.0,
                    C,
                    3
                );
                memcpy(transformMap[face][hIndex][vIndex].matrix, C, 9 * sizeof(double));
                cblas_dgemm(
                    CblasRowMajor,
                    CblasNoTrans,
                    CblasNoTrans,
                    3, 3,
                    3,
                    1.0,
                    backTranslate,
                    3,
                    transformMap[face][hIndex][vIndex].matrix,
                    3,
                    0.0,
                    C,
                    3
                );
                memcpy(transformMap[face][hIndex][vIndex].matrix, C, 9 * sizeof(double));

                transformMap[face][hIndex][vIndex].face = verticalTransforms[horizontalTransforms[face][hIndex].face][vIndex].face;
            }
        }
    }
}

long CubeWorld::getDisplacedFaces(long position) {
    if (position < 0) {
        return ((position + 1) / ((long)size)) - 1;
    }
    if (position >= (long)size) {
        return position / ((long)size);
    }
    return 0;
}

void CubeWorld::navigate(char *face, long *x, long *y) {
    double xy[3] = {
        *x - getDisplacedFaces(*x) * ((long)size),
        *y - getDisplacedFaces(*y) * ((long)size),
        1.0
    };
    double result[3];

    long hIndex = getDisplacedFaces(*x) % 4;
    if (hIndex < 0) {
        hIndex += 4;
    }
    long vIndex = getDisplacedFaces(*y) % 4;
    if (vIndex < 0) {
        vIndex += 4;
    }

    cblas_dgemv(
        CblasRowMajor,
        CblasNoTrans,
        3, 3,
        1.0,
        transformMap[*face][hIndex][vIndex].matrix,
        3,
        xy,
        1,
        0.0,
        result,
        1
    );

    *y = result[1];
    *x = result[0];
    *face = transformMap[*face][hIndex][vIndex].face;
}

CubeWorld::CubeWorld(size_t size) : size(size) {
    for (size_t i = 0; i < 6; i++) {
        data[i] = gsl_matrix_alloc(size, size);
    }
    initializeMatricies();
}

CubeWorld::~CubeWorld() {
    for (size_t i = 0; i < 6; i++) {
        gsl_matrix_free(data[i]);
    }
}

size_t CubeWorld::getSize() {
    return size;
}

double CubeWorld::get(char face, long x, long y) {
    navigate(&face, &x, &y);
    return gsl_matrix_get(data[face], x, y);
}

void CubeWorld::set(char face, long x, long y, double value) {
    navigate(&face, &x, &y);
    gsl_matrix_set(data[face], x, y, value);
}

void CubeWorld::setAll(double value) {
    for (char face = 0; face < 6; face++) {
        for (size_t x = 0; x < size; x++) {
            for (size_t y = 0; y < size; y++) {
                gsl_matrix_set(data[face], x, y, value);
            }
        }
    }
}

void printTerrain(std::string title, CubeWorld& map) {
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

void printTerrain(CubeWorld& map) {
    printTerrain(std::string("output.bmp"), map);
};