/// @brief Contains auxiliar functions (most of them come from numpy librairy)
#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>

template <typename T>
/**
 * Numpy linspace function
 * @param  a first T element
 * @param  b last T element
 * @param  N size of the vector
 * @return   vector of N T type elements bet a and b
 */
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

#endif