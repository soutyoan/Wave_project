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

/*
Our own implementation of the ft Shift
Inversing blocks in the Matrix
Matrix :     1 2
             3 4

fftshift :   4 3
             2 1

WARNING : TODO : odd sizes
*/
void ft_shift(Mat &in, Mat &out){
    // We need to copy each block
    int middle_x = (int)(in.rows/2); 
    int middle_y = (int)(in.cols/2); 

    // Block 1 to 4
    for (int i = 0; i < middle_x; i++){
        for (int j = 0; j < middle_y; j++){
            out.at<float>(i + middle_x, j + middle_y) = in.at<float>(i, j); 
        }
    }

    // Block 4 to 1
    for (int i = 0; i < middle_x; i++){
        for (int j = 0; j < middle_y; j++){
            out.at<float>(i, j) = in.at<float>(i + middle_x, j + middle_y); 
        }
    }

    // Block 3 to 2
    for (int i = 0; i < middle_x; i++){
        for (int j = 0; j < middle_y; j++){
            out.at<float>(i + middle_x, j) = in.at<float>(i, j + middle_y); 
        }
    }

    // Block 2 to 3
    for (int i = 0; i < middle_x; i++){
        for (int j = 0; j < middle_y; j++){
            out.at<float>(i, j + middle_y) = in.at<float>(i + middle_x, j); 
        }
    }
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
Mat* polar_coordinates(const vector<T> &xx, const vector<T> &yy, size_t nx, size_t ny){
    Mat* res = new Mat(nx, ny, CV_64FC1);
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            // cout << i << " " << j << " " << xx[i] << "\n";
            res->at<T>(i,j) = sqrt(pow(xx[i], 2) + pow(yy[j], 2));
        }
    }
    return res;
}

template <typename T>
/**
 * Angular coordinates grid computation
 * @param  wx line vector of meshgrid
 * @param  wy column vector of meshgrid
 * @param  nx size of wx vector
 * @param  ny size of wy vector
 * @return    OpenCV Matrix with angular coordinates
 */
Mat *angular_coordinates(const vector<T> &wx, const vector<T> &wy, size_t nx, size_t ny) {
    Mat *res = new Mat(nx, ny, CV_64FC1);
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
