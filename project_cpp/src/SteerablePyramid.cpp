#include "SteerablePyramid.h"

#define _USE_MATH_DEFINES
#include <math.h>

SteerablePyramid::SteerablePyramid(Mat *_image,
                int _xRes,
                int _yRes,
                int _n,
                int _k,
                string _image_name,
                string _output_path,
                bool _verbose){
    this->image = _image;
    this->xRes = _xRes;
    this->yRes = _yRes;
    this->n = _n;
    this->k = _k;
    this->image_name = _image_name;
    this->output_path = _output_path;
    this->verbose = _verbose;

}

vector<vector<float> >* calicurate_polar(){

}

Mat* calicurate_h0_filter(){
    vector<vector<Mat> > * polar = calicurate_polar();
    Mat RS = polar->at(0)[0];
    Mat *fil = new Mat(size(xRes, yRes), CV_64FC1, 0);
    // Possible parallelisation?
    for (int i = 0; i < xRes; i++){
        for (int j = 0; j < yRes; j++){
            if (RS.at<float>(i, j) >= M_PI){
                fil->at<float>(i, j) = 1;
            } else if (RS.at<float>(i, j) >= M_PI/2) {
                fil->at<float>(i, j) = cos(M_PI/2 * log2(RS.at<float>(i, j)/M_PI));
            }
        }
    }
    return fil;
}

Mat* calicurate_l0_filter(){
    vector<vector<Mat> > * polar = calicurate_polar();
    Mat RS = polar->at(0)[0];
    Mat *fil = new Mat(size(xRes, yRes), CV_64FC1, 0);
    // Possible parallelisation?
    for (int i = 0; i < xRes; i++){
        for (int j = 0; j < yRes; j++){
            if (RS.at<float>(i, j) >= 0){
                fil->at<float>(i, j) = 1;
            } else if (RS.at<float>(i, j) >= M_PI/2) {
                fil->at<float>(i, j) = cos(M_PI/2 * log2(RS.at<float>(i, j)/M_PI));
            }
        }
    }
    return fil;
}

vector<Mat* > * calicurate_l_filter(){
    vector<vector<Mat> > * polar = calicurate_polar();
    vector<Mat *> * f = new vector<Mat* >(n);
    for (int i = 0; i < n; i++){
        Mat* m = new Mat(size(xRes, yRes), CV_64FC1, 0);
        for (int x = 0; x < xRes; x++){
            for(int y = 0; y < yRes; y++){
                if (RS.at<float>(i, j) >= M_PI/2){
                    m->at<float>(i, j) = 0;
                } else if (RS.at<float>(i, j) <= M_PI/4){
                    m->at<float>(i, j) = 1;
                } else {
                    m->at<float>(i, j) = cos(M_PI/2 * log2(4 * RS.at<float>(i, j)/M_PI));
                }
            }
        }
        f[n] = m;
    }
    return f;
}

vector<Mat* > * calicurate_h_filter(){
    vector<vector<Mat> > * polar = calicurate_polar();
    vector<Mat *> * f = new vector<Mat* >(n);
    for (int i = 0; i < n; i++){
        Mat* m = new Mat(size(xRes, yRes), CV_64FC1, 0);
        for (int x = 0; x < xRes; x++){
            for(int y = 0; y < yRes; y++){
                if (RS.at<float>(i, j) >= M_PI/2){
                    m->at<float>(i, j) = 1;
                } else if (RS.at<float>(i, j) <= M_PI/4){
                    m->at<float>(i, j) = 0;
                } else {
                    m->at<float>(i, j) = cos(M_PI/2 * log2(2 * RS.at<float>(i, j)/M_PI));
                }
            }
        }
        f[n] = m;
    }
    return f;
}

Mat* calicurate_b_filter(){

}

void createPyramids(){

}

Mat* collapsePyramids(){

}

Mat* clearPyramids(){

}

SteerablePyramid::~SteerablePyramid(){
    this->image->release();

}
