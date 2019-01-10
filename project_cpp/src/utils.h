/// @brief Contains auxiliar functions (most of them come from numpy librairy)
#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <math.h>

# define M_PI           3.14159265358979323846


using namespace cv;
using namespace std;

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


template <typename T>
/**
 * Polar coordinates grid computation
 * @param  xx line vector of meshgrid
 * @param  yy column vector of meshgrid
 * @param  nx size of xx vector
 * @param  ny size of yy vector
 * @return opencv Matrix with polar coordinates
 */
cv::Mat* polar_coordinates(const vector<T> &xx, const vector<T> &yy, size_t nx, size_t ny)
    cv::Mat* res = new Mat(nx, ny, CV_32F);
    for (int i=0; i<nx; i++) {
        for (int j=0; i<ny; j++) {
            res->at<T>(i,j) = sqrt(xx[i]**2 + yy[j]**2);
        }
    }
    return res;

template <typename T>
/**
 * Angular coordinates grid computation
 * @param  wx line vector of meshgrid
 * @param  wy column vector of meshgrid
 * @param  nx size of wx vector
 * @param  ny size of wy vector
 * @return    OpenCV Matrix with angular coordinates
 */
cv::Mat *angular_coordinates(const vector<T> &wx, const vector<T> &wy, size_t nx, size_t ny) {
    cv::Mat *res = new Mat(nx, ny, CV_32F);
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            if (wy[j]==0.0f && wx[i]<0) {
                res->at<T>(i,j) = M_PI;
            }
            else {
                res->at<T>(i,j) = atan2(wy[j], wx[i]);
            }
        }
    }
    return res;
}

#endif
