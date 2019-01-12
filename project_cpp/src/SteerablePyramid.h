#ifndef STEERABLE_PYRAMID_H
#define STEERABLE_PYRAMID_H

#include <vector>
#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>
#include "utils.h"

using namespace cv;
using namespace std;

class SteerablePyramid{

public:
    // Matrix of the image
    Mat image;
    // Horizontal resolution
    int xRes;
    // Vertical resolution
    int yRes;
    // Depth, number of scales
    int n;
    // Number of orientations (must be in (4, 6, 8, 10, 12, 15, 18, 20, 30, 60))
    int k;
    // Name of the chosen image
    string image_name;
    // Name of the ouput path of the image
    string output_path;
    // Can set the program as verbose
    bool verbose;

    SteerablePyramid(Mat image,
                    int xRes,
                    int yRes,
                    int n,
                    int k,
                    string image_name,
                    string output_path,
                    bool verbose);

    ~SteerablePyramid();

    void caliculate_one_polar(Mat RS, Mat AT, int i); 

    void caliculate_polar(vector<Mat > &RS, vector<Mat > &AT);

    // caliculate H0 values on the grid.
    void calicurate_h0_filter(Mat &fil, vector<Mat > &RS);

    // caliculate L0 values on the grid.
    void calicurate_l0_filter(Mat &fil, vector<Mat > &RS);

    // caliculate L filter values on the grid.
    void calicurate_l_filter(int i, Mat &f, vector<Mat > &RS);

    // caliculate H filter values on the grid.
    void calicurate_h_filter(vector<Mat> &f, vector<Mat> &RS);

    // caliculate B filter values on the grid.
    void calicurate_b_filter(int i, int j, Mat &fil, const vector<Mat> &AT);

    // create steerable pyramid
    void createPyramids();

    // image reconstruction from steerable pyramid in Fourier domain.
    void collapsePyramids(Mat &f);

    // clear the steerable pyramid
    void clearPyramids(Mat &f);

private:
    float alphak;


};

#endif
