#ifndef STEERABLE_PYRAMID_H
#define STEERABLE_PYRAMID_H

#include <vector>
#include <opencv2/opencv.hpp>
#include <string>
#include <iostream>

using namespace cv;
using namespace std;

class SteerablePyramid{

public:
    // Matrix of the image
    Mat *image;
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

    SteerablePyramid(Mat *image,
                    int xRes,
                    int yRes,
                    int n,
                    int k,
                    string image_name,
                    string output_path,
                    bool verbose);

    ~SteerablePyramid();

    vector<vector<float> >* calicurate_polar();

    // caliculate H0 values on the grid.
    Mat* calicurate_h0_filter();

    // caliculate L0 values on the grid.
    Mat* calicurate_l0_filter();

    // caliculate L filter values on the grid.
    vector<Mat*> * calicurate_l_filter();

    // caliculate H filter values on the grid.
    vector<Mat*> * calicurate_h_filter();

    // caliculate B filter values on the grid.
    Mat calicurate_b_filter();

    // create steerable pyramid
    void createPyramids();

    // image reconstruction from steerable pyramid in Fourier domain.
    Mat* collapsePyramids();

    // clear the steerable pyramid
    Mat* clearPyramids();

private:
    float alphak;


};

#endif
