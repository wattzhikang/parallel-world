#include <math.h>

#include "projector.hpp"

namespace projector
{
    /*
        All formulas from Wolfram MathWorld
        https://mathworld.wolfram.com/GnomonicProjection.html
    */

    #define LAT 0
    #define LON 1
    float faceCenters[6][2] = {
        // LAT     LON
        {  0.0f,   0.0f},//face0
        {  0.0f,  M_PI_2},//face1
        {  0.0f, M_PI},//face2
        {  0.0f, -M_PI_2},//face3
        {-M_PI_2,   0.0f},//face4
        { M_PI_2,   0.0f} //face5
    };

    float longitudeZones[4] = {
        M_PI_4, 3 * M_PI_4, 5 * M_PI_4, 7 * M_PI_4
    };

    float p(float x, float y) {
        return sqrt(x * x + y * y);
    }

    float c(float p) {
        return atan(p);
    }

    float getLongitude(char face, float xf, float yf, float sinC, float cosCTmp, float tmpP) {
        return faceCenters[face][LAT] + atan2(
            ((float)xf * sinC),
            (tmpP * cos(faceCenters[face][LON]) * cosCTmp - yf * sin(faceCenters[face][LON]) * sinC)
        );
    }

    float getLatitude(char face, float yf, float sinC, float cosCTmp, float tmpP) {
        return asin(
            cosCTmp * sin(faceCenters[face][LON]) +
            ((yf != 0) ? ((float)yf * sinC * cos(faceCenters[face][LON])) / tmpP : 0)
        );
    }

    float cosC(char face, float lat, float lon) {
        return 
            sin(faceCenters[face][LAT]) * sin(lat) +
            cos(faceCenters[face][LAT]) * cos(lat) * cos(lon - faceCenters[face][LON])
        ;
    }

    // This is easy because we already know what face it is
    void getSpherical(CubeWorld& map, char face, size_t x, size_t y, float *lat, float *lon) {
        float size = (float)map.getSize();
        float xf = (float)x / size - 0.5f;
        float yf = (float)y / size - 0.5f;

        float tmpP = p((float)xf, (float)yf);
        float cSto = c(tmpP);

        float sinC = sin(cSto);
        float cosCTmp = cos(cSto);

        *lon = getLongitude(face, xf, yf, sinC, cosCTmp, tmpP);

        *lat = getLatitude(face, yf, sinC, cosCTmp, tmpP);
    }

    void getCartesian(CubeWorld& map, float lat, float lon, char *face, size_t *x, size_t *y) {
        //first, determine the face
        char tmpFace;
        //determine longitude zone
        //first normalize longitude
        while (lon < 0) {
            lon += 2 * M_PI;
        }
        while (lon > 2 * M_PI) {
            lon -= 2 * M_PI;
        }
        //zone determination
        if (lon < M_PI_4) {
            tmpFace = 0;
        } else if (lon < 3 * M_PI_4) {
            tmpFace = 1;
        } else if (lon < 5 * M_PI_4) {
            tmpFace = 2;
        } else if (lon < 7 * M_PI_4) {
            tmpFace = 3;
        } else {
            tmpFace = 0;
        }
        // //determine if it is an upper, lower, or middle face based on latitude
        // if (lat < -0.463647f) {
        //     //if the latitude is lower than the highest possible latitude for the south polar face (face 4), then
        //     //it could be in the south polar face. Check to see if the point's latitude, at the point's longitude,
        //     //fall within the south polar face
        //     float cSto = c(1.0f);
        //     if (lat < getLatitude(0, -1.0f, sin(cSto), cos(cSto), 1.0f)) {
        //         tmpFace = 4;
        //     }
        // } else if ( lat > 0.463647f) {
        //     //if the latitude is higher than the lowest possible latitude for the north polar face (face 5), do
        //     //basically the same thing
        //     float cSto = c(1.0f);
        //     if (lat > getLatitude(0, 1.0f, sin(cSto), cos(cSto), 1.0f)) {
        //         tmpFace = 5;
        //     }
        // }

        float cosCTmp = cosC(tmpFace, lat, lon);
        float tmpy =
            (cos(faceCenters[tmpFace][LAT]) * sin(lat) - sin(faceCenters[tmpFace][LAT]) * cos(lat) * cos(lon - faceCenters[tmpFace][LON]))
            / cosCTmp
        ;

        if (tmpy < -1 || lat < -M_PI_4) {
            tmpFace = 4;
            cosCTmp = cosC(tmpFace, lat, lon);
            tmpy =
                (cos(faceCenters[tmpFace][LAT]) * sin(lat) - sin(faceCenters[tmpFace][LAT]) * cos(lat) * cos(lon - faceCenters[tmpFace][LON]))
                / cosCTmp
            ;
        } else if (tmpy > 1 || lat > M_PI_4) {
            tmpFace = 5;
            cosCTmp = cosC(tmpFace, lat, lon);
            tmpy =
                (cos(faceCenters[tmpFace][LAT]) * sin(lat) - sin(faceCenters[tmpFace][LAT]) * cos(lat) * cos(lon - faceCenters[tmpFace][LON]))
                / cosCTmp
            ;
        }

        float tmpx = ( cos(lat) * sin(lon - faceCenters[tmpFace][LON]) ) / cosCTmp;

        *x = (size_t)((((-tmpx) + 1.0f)*0.5) * map.getSize());
        *y = (size_t)((((-tmpy) + 1.0f)*0.5) * map.getSize());
        if (tmpFace == 1) {
            *face = 3;
        } else if (tmpFace == 3) {
            *face = 1;
        } else {
            *face = tmpFace;
        }
    }
} // namespace projector
