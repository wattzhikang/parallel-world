#include "Unit.hpp"

template<class T> class ParallelBlock {
    public:
        /**
         * This method gets the cell that this latitude and longitude
         * vector falls within.
         * 
         * @param latitude The latitude, in radians
         * @param longitude The longitude, in radians
        */ 
        virtual Unit<T>& getCell(double latitude, double longitude);
        /**
         * This method gets the cell where x and y form a vector
         * pointing from the global origin. Note that this gets
         * a cell that may be outside this ParallelBlock, and
         * therefore access is not guaranteed. For best results
         * use getLocalCell() when possible.
         * 
         * @param x The x coordinate, from the global grid origin.
         * @param y The y coordinate, from the global grid origin.
        */
        virtual Unit<T>& getCell(size_t x, size_t y);

        /**
         * This method accesses the cell pointed to by the vector
         * UV, centered at the origin from this ParallelBlock.
         * 
         * @param u The u coordinate, from the local block origin
         * @param v The v coordinate, from the local block origin
        */ 
        virtual Unit<T>& getLocalCell(size_t u, size_t v);
};