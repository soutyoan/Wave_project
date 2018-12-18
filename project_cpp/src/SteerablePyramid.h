#ifndef STEERABLE_PYRAMID_H
#define STEERABLE_PYRAMID_H

#include <stdio>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace cv;

class SteerablePyramid{

public:
    // Matrix of the image
    Mat *image;
    // Horizontal resolution
    int XRes;
    // Vertical resolution
    int YRes;
    // Depth, number of scales
    int n;
    // Number of orientations (must be in (4, 6, 8, 10, 12, 15, 18, 20, 30, 60))
    int k;
    // Name of the chosen image
    string image_name:
    // Name of the ouput path of the image
    string out_path;

    SteerablePyramid(Mat *image,
                    int xRes,
                    int yRes,
                    int n,
                    int k,
                    string image_name,
                    string output_path,
                    bool verbose);

    ~SteerablePyramid();


};

#endif
