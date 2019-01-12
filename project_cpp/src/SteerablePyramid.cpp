#include "SteerablePyramid.h"

#define _USE_MATH_DEFINES
#include <math.h>

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

SteerablePyramid::SteerablePyramid(Mat _image,
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

    cout << xRes << "\n"; 
    cout << image.rows << "\n"; 

    assert(xRes%2 == 0); // Some function like ft shift not implemented 
    // for odd numbers
    assert(yRes%2 == 0); 
}

/*
in order to do a good parallelization we need to be able to calculate these
matrices one at a time
*/
void SteerablePyramid::caliculate_one_polar(Mat RS, Mat AT, int i){
    float _tmp = pow(2.0, i);
    size_t nx = this->image.rows / _tmp;
    size_t ny = this->image.cols / _tmp;
    vector<float> _wx = linspace<float>(-M_PI, M_PI, nx);
    vector<float> _wy = linspace<float>(-M_PI, M_PI, ny);
    cout << _wx.size() << " " << _wy.size() << " " << nx << " " << ny << "\n"; 
    polar_coordinates<float>(_wx, _wy, nx, ny, RS);
    angular_coordinates<float>(_wx, _wy, nx, ny, AT);
}

void SteerablePyramid::caliculate_polar(vector<Mat> &RS, vector<Mat> &AT){

    for (int i=0; i < this->n; i++) {
        // Computation polar coordinates (radius) on the grid
        float _tmp = pow(2.0, i);
        Mat _rs(Size(this->image.rows / _tmp, this->image.cols / _tmp), CV_64FC1); 
        Mat _at(Size(this->image.rows / _tmp, this->image.cols / _tmp), CV_64FC1); 
        caliculate_one_polar(_rs, _at, i); 
        cout << "RS " << _rs.rows << " " << _rs.cols << "\n"; 
        RS.push_back(_rs);
        AT.push_back(_at);
    }
}

void SteerablePyramid::calicurate_h0_filter(Mat &fil, vector<Mat> &RS){
    Mat RS_0 = RS[0];
    cout << "size" << RS_0.rows << "\n";  
    // Mat fil = new Mat(Size(xRes, yRes), CV_64FC1);
    // Possible parallelisation?
    for (int i = 0; i < xRes; i++){
        for (int j = 0; j < yRes; j++){
            if (RS_0.at<float>(i, j) >= M_PI){
                fil.at<float>(i, j) = 1;
            } else if (RS_0.at<float>(i, j) >= M_PI/2) {
                fil.at<float>(i, j) = cos(M_PI/2 * log2(RS_0.at<float>(i, j)/M_PI));
            }
        }
    }
}

void SteerablePyramid::calicurate_l0_filter(Mat &fil, vector<Mat> &RS){
    Mat RS_0 = RS[0];
    // Possible parallelisation?
    for (int i = 0; i < xRes; i++){
        for (int j = 0; j < yRes; j++){
            if (RS_0.at<float>(i, j) >= 0){
                fil.at<float>(i, j) = 1;
            } else if (RS_0.at<float>(i, j) >= M_PI/2) {
                fil.at<float>(i, j) = cos(M_PI/2 * log2(RS_0.at<float>(i, j)/M_PI));
            }
        }
    }
}

void SteerablePyramid::calicurate_l_filter(int i, Mat &f, vector<Mat> &RS){
    Mat RS_0 = RS[0];
    for (int x = 0; x < xRes; x++){
        for(int y = 0; y < yRes; y++){
            if (RS_0.at<float>(x, y) >= M_PI/2){
                f.at<float>(x, y) = 0;
            } else if (RS_0.at<float>(x, y) <= M_PI/4){
                f.at<float>(x, y) = 1;
            } else {
                f.at<float>(x, y) = cos(M_PI/2 * log2(4 * RS_0.at<float>(x, y)/M_PI));
            }
        }
    }
}

void SteerablePyramid::calicurate_h_filter(vector<Mat> &f, vector<Mat> &RS){
    Mat RS_0 = RS[0];
    for (int i = 0; i < n; i++){
        Mat m(Size(xRes, yRes), CV_64FC1, 0);
        for (int x = 0; x < xRes; x++){
            for(int y = 0; y < yRes; y++){
                if (RS_0.at<float>(x, y) >= M_PI/2){
                    m.at<float>(x, y) = 1;
                } else if (RS_0.at<float>(x, y) <= M_PI/4){
                    m.at<float>(x, y) = 0;
                } else {
                    m.at<float>(x, y) = cos(M_PI/2 * log2(2 * RS_0.at<float>(x, y)/M_PI));
                }
            }
        }
    }
}

// On calcule les b filters un à un
void SteerablePyramid::calicurate_b_filter(int i, int j, Mat &fil, const vector<Mat> &AT){

    for (int x = 0; x < xRes; x++){
        for(int y = 0; y < yRes; y++){
            float currentCoefficient = AT[i].at<float>(x, y);
            if (AT[i].at<float>(x, y) - j * M_PI/k < -M_PI){
                currentCoefficient += 2*M_PI;
            } else if (AT[i].at<float>(x, y) - j * M_PI/k > M_PI){
                currentCoefficient -= 2*M_PI;
            } else if (abs(AT[i].at<float>(x, y) - j*M_PI/k)< M_PI/2){
                currentCoefficient = alphak * pow(cos(AT[i].at<float>(x, y)
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

    vector<Mat> RS; // calculate RS
    vector<Mat> AT; // calculate AT

    this->caliculate_polar(RS, AT);

    // Create all the matrices used during the collapse function
    Mat h0f(Size(xRes, yRes), CV_64FC1);
    Mat h0s(Size(xRes, yRes), CV_64FC1);
    Mat l0f(Size(xRes, yRes), CV_64FC1);
    Mat l0s(Size(xRes, yRes), CV_64FC1);

    // Find library for fast fourier transform
    Mat ft(Size(xRes, yRes), CV_64FC1); 
    Mat _ft(Size(xRes, yRes), CV_64FC1); 

    ///// THREAD 1 in openmp use thread ID //////

    Mat imgBack(Size(xRes, yRes), CV_64FC1); 
    image.convertTo(image, CV_64FC1);

    cout << "test" << image.type(); 

    cout << "hey" << (image).at<float>(5, 25) << "\n"; 

    dft((image), ft); 
    ft_shift(ft, _ft); 

    // Create HO filter
    Mat h0(Size(xRes, yRes), CV_64FC1);
    calicurate_h0_filter(h0, RS);

    h0 = h0.mul(_ft);  // calculation of h0
    Mat f_ishift(Size(xRes, yRes), CV_64FC1);  
    ft_shift(h0, f_ishift); // FFT opencv
    idft(f_ishift, imgBack); // IFFT opencv

    //// THREAD 2 in openmp use thread ID /////

    Mat l0(Size(xRes, yRes), CV_64FC1);
    calicurate_l0_filter(l0, RS);

    l0 = l0.mul(_ft);  // calculation of h0 

    dft(h0, f_ishift); // FFT opencv
    idft(f_ishift, imgBack); // IFFT opencv

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

}
