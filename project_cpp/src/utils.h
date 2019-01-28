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
template <typename T>
void ft_shift(Mat &in, Mat &out){
    // We need to copy each block
    int middle_x = (int)(in.rows/2);
    int middle_y = (int)(in.cols/2);

    // Block 1 to 4
    for (int i = 0; i < middle_x; i++){
        for (int j = 0; j < middle_y; j++){
            out.at<T>(i + middle_x, j + middle_y) = in.at<T>(i, j);
        }
    }

    // Block 4 to 1
    for (int i = 0; i < middle_x; i++){
        for (int j = 0; j < middle_y; j++){
            out.at<T>(i, j) = in.at<T>(i + middle_x, j + middle_y);
        }
    }

    // Block 3 to 2
    for (int i = 0; i < middle_x; i++){
        for (int j = 0; j < middle_y; j++){
            out.at<T>(i + middle_x, j) = in.at<T>(i, j + middle_y);
        }
    }

    // Block 2 to 3
    for (int i = 0; i < middle_x; i++){
        for (int j = 0; j < middle_y; j++){
            out.at<std::complex<double> >(i, j + middle_y) = in.at<std::complex<double> >(i + middle_x, j);
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
void polar_coordinates(const vector<T> &xx, const vector<T> &yy,
    size_t nx, size_t ny, Mat res){
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            // cout << i << " " << j << " " << xx[i] << "\n";
            res.at<T>(i,j) = sqrt(pow(xx[i], 2) + pow(yy[j], 2));
        }
    }
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
void angular_coordinates(const vector<T> &wx, const vector<T> &wy,
    size_t nx, size_t ny, Mat res) {
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            if (wy[j]==0.0f && wx[i]<0) {
                res.at<T>(i,j) = M_PI;
            }
            else {
                res.at<T>(i,j) = atan2(wy[j], wx[i]);
            }
        }
    }
}

Mat mul_complex(Mat &m1, Mat &m2){
    assert (m1.rows == m2.rows);
    assert (m1.cols == m2.cols);

    cout << "ASSERTIONS OK" << endl;

    Mat m = Mat_<std::complex<double> >(Size(m1.rows, m1.cols), CV_64FC1);

    for (int i = 0; i < m1.rows; i++){
        for (int j = 0; j < m1.cols; j++){
            m.at<complex<double> >(i, j) = m1.at<double>(i, j) * m2.at<complex<double> >(i, j);
        }
    }

    return m;
}

#endif
