#include "SteerablePyramid.h"

#define _USE_MATH_DEFINES
#include <math.h>

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

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

    this->alphak = 	pow(2, (k-1)) * factorial(k-1)/sqrt(k * factorial(2*(k-1)));
}

void SteerablePyramid::caliculate_polar(vector<Mat*> &RS, vector<Mat *> &AT){

    for (int i=0; i<this->n; i++) {
        // Computation polar coordinates (radius) on the grid
        float _tmp = pow(2.0, i);
        size_t grid_x = this->image->rows / _tmp;
        size_t grid_y = this->image->rows / _tmp;
        vector<float> _wx = linspace<float>(-M_PI, M_PI, grid_x);
        vector<float> _wy = linspace<float>(-M_PI, M_PI, grid_y);
        Mat *_rs = new Mat(grid_x, grid_y, CV_32F);
    }

}

void SteerablePyramid::calicurate_h0_filter(Mat &fil, vector<Mat*> &RS){
    Mat* RS_0 = RS[0];
    // Mat* fil = new Mat(Size(xRes, yRes), CV_64FC1);
    // Possible parallelisation?
    for (int i = 0; i < xRes; i++){
        for (int j = 0; j < yRes; j++){
            if (RS_0->at<float>(i, j) >= M_PI){
                fil.at<float>(i, j) = 1;
            } else if (RS_0->at<float>(i, j) >= M_PI/2) {
                fil.at<float>(i, j) = cos(M_PI/2 * log2(RS_0->at<float>(i, j)/M_PI));
            }
        }
    }
}

void SteerablePyramid::calicurate_l0_filter(Mat &fil, vector<Mat*> &RS){
    Mat* RS_0 = RS[0];
    // Possible parallelisation?
    for (int i = 0; i < xRes; i++){
        for (int j = 0; j < yRes; j++){
            if (RS_0->at<float>(i, j) >= 0){
                fil.at<float>(i, j) = 1;
            } else if (RS_0->at<float>(i, j) >= M_PI/2) {
                fil.at<float>(i, j) = cos(M_PI/2 * log2(RS_0->at<float>(i, j)/M_PI));
            }
        }
    }
}

void SteerablePyramid::calicurate_l_filter(int i, Mat &f, vector<Mat*> &RS){
    Mat* RS_0 = RS[0];
    for (int x = 0; x < xRes; x++){
        for(int y = 0; y < yRes; y++){
            if (RS_0->at<float>(x, y) >= M_PI/2){
                f.at<float>(x, y) = 0;
            } else if (RS_0->at<float>(x, y) <= M_PI/4){
                f.at<float>(x, y) = 1;
            } else {
                f.at<float>(x, y) = cos(M_PI/2 * log2(4 * RS_0->at<float>(x, y)/M_PI));
            }
        }
    }
}

void SteerablePyramid::calicurate_h_filter(vector<Mat *> &f, vector<Mat*> &RS){
    Mat* RS_0 = RS[0];
    for (int i = 0; i < n; i++){
        Mat* m = new Mat(Size(xRes, yRes), CV_64FC1, 0);
        for (int x = 0; x < xRes; x++){
            for(int y = 0; y < yRes; y++){
                if (RS_0->at<float>(x, y) >= M_PI/2){
                    m->at<float>(x, y) = 1;
                } else if (RS_0->at<float>(x, y) <= M_PI/4){
                    m->at<float>(x, y) = 0;
                } else {
                    m->at<float>(x, y) = cos(M_PI/2 * log2(2 * RS_0->at<float>(x, y)/M_PI));
                }
            }
        }
    }
}

// On calcule les b filters un à un
void SteerablePyramid::calicurate_b_filter(int i, int j, Mat &fil, const vector<Mat *> &AT){

    for (int x = 0; x < xRes; x++){
        for(int y = 0; y < yRes; y++){
            float currentCoefficient = AT[i]->at<float>(x, y);
            if (AT[i]->at<float>(x, y) - j * M_PI/k < -M_PI){
                currentCoefficient += 2*M_PI;
            } else if (AT[i]->at<float>(x, y) - j * M_PI/k > M_PI){
                currentCoefficient -= 2*M_PI;
            } else if (abs(AT[i]->at<float>(x, y) - j*M_PI/k)< M_PI/2){
                currentCoefficient = alphak * pow(cos(AT[i]->at<float>(x, y)
                - j * M_PI/k), k-1);
            }
        }
    }
}

void fillDownfillMatrix(Mat &down, int x, int y){
    for (int i = x; i < 3 * x; i++){
        for (int j = y; j < 3 * y; j++){
            down.at<float>(i, j) = 1;
        }
    }
}

void SteerablePyramid::createPyramids(){

    vector<Mat *> RS; // calculate RS
    vector<Mat *> AT; // calculate AT

    // Create all the matrices used during the collapse function
    Mat h0f;
    Mat h0s;
    Mat l0f;
    Mat l0s;

    // Find library for fast fourier transform
    Mat ft;
    dct((*image), ft);
    Mat ft_shift;

    ///// THREAD 1 in openmp use thread ID //////

    // Create HO filter
    Mat h0(Size(xRes, yRes), CV_64FC1);
    calicurate_h0_filter(h0, RS);

    h0 = h0.mul(ft_shift);  // calculation of h0
    Mat f_ishift;  //ishift


    //// THREAD 2 in openmp use thread ID /////

    Mat l0(Size(xRes, yRes), CV_64FC1);
    calicurate_l0_filter(l0, RS);

    Mat lastImage = Mat_<std::complex<float> >(l0);

    vector<Mat *> BND;

    // Parallelize here (use openmp for)
    for (int i = 0; i < n; i ++){

        for (int j = 0; j < yRes; j++){

            Mat b_filter(Size(xRes, yRes), CV_64FC1);
            calicurate_b_filter(i, j, b_filter, AT);

            Mat lb = lastImage.mul(b_filter);
            Mat f_ishift; // Shift lb
            Mat img_back; // ishift2 f_ishift

            BND.push_back(&lb);
            BND.push_back(&img_back);

        }

        // Apply low pas filter to image downsampled
        int img_x = lastImage.size[0];
        int img_y = lastImage.size[1];
        int quantification_x = (int)(img_x/4);
        int quantification_y = (int)(img_y/4);

        Mat down_fil(Size(img_x, img_y), CV_64FC1);
        fillDownfillMatrix(down_fil, quantification_x, quantification_y);

        Mat l(Size(xRes, yRes), CV_64FC1);
        calicurate_l_filter(i, l, RS);

        down_fil = lastImage.mul(l).mul(down_fil);

        Mat down_image = Mat_<std::complex<double> >(2 * quantification_x, 2 * quantification_y);
        for (int i = quantification_x; i < 3 * quantification_x; i ++){
            for (int j = quantification_y; j < 3* quantification_y; j++){
                down_image.at<complex<double> >(i, j) = down_fil.at<complex<double> >(i, j);
            }
        }

        // Par rapport à la version python, on ne sauvegarde pas low.
        // Ca ne semble pas avoir d'interet.

        lastImage = down_image;
    }

    Mat LRf = lastImage;

}

void collapsePyramids(Mat &f){

}

void SteerablePyramid::clearPyramids(Mat &f){

}

SteerablePyramid::~SteerablePyramid(){
    this->image->release();

}
