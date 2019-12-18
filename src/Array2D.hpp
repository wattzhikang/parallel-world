#ifndef ARR_2D
#define ARR2D

#include <iostream>
#include <math.h>

/*
2D array syntax is difficult in C++, especially when having to pass
arrays between functions. So this class basically emulates the
behavior of a 2D array.
*/
template<typename T> class Array2D {
    private:
        size_t size;
        T *data;
    public:
        Array2D(size_t size) : size(size) {
            try {
                data = new T[size * size];
            }
            catch (std::bad_alloc &ba) {
                std::cerr << "bad_alloc caught: " << ba.what();
            }
        };
        ~Array2D() {
            delete[] data;
        }
        T operator() (size_t i, size_t j) {
            return data[i * size + j];
        }
        void set(size_t i, size_t j, T datum) {
            data[i * size + j] = datum;
        }
        void set(T toSet) {
            size_t i, sSquared = (size_t) pow(size, 2);
            for (i = 0; i < sSquared; i++)
            {
                data[i] = toSet;
            }
            
        }
        size_t getSize() {
            return size;
        }
};

#endif