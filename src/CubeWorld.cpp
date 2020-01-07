#include "CubeWorld.hpp"

long CubeWorld::getDisplacedFaces(long position) {
    if (position < 0) {
        return ((position + 1) / ((long)size)) - 1;
    }
    if (position >= (long)size) {
        return position / ((long)size);
    }
    return 0;
}

char CubeWorld::displacementIndex(long facesDisplaced) {
    if (facesDisplaced < 0) {
        return (facesDisplaced % 4) + 4;
    } else {
        return facesDisplaced % 4;
    }
}

void CubeWorld::navigate(char *face, long *x, long *y) {
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

CubeWorld::CubeWorld(size_t size) : size(size) {
    for (size_t i = 0; i < 6; i++) {
        data[i] = gsl_matrix_alloc(size, size);
    }
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