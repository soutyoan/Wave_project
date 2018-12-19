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

vector<vector<float> >* SteerablePyramid::calicurate_polar(){

}

Mat* SteerablePyramid::calicurate_h0_filter(){
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

Mat* SteerablePyramid::calicurate_l0_filter(){
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

vector<Mat* > * SteerablePyramid::calicurate_l_filter(){
    vector<vector<Mat> > * polar = calicurate_polar();
    vector<Mat *> * f = new vector<Mat* >(n);
    for (int i = 0; i < n; i++){
        Mat* m = new Mat(size(xRes, yRes), CV_64FC1, 0);
        for (int x = 0; x < xRes; x++){
            for(int y = 0; y < yRes; y++){
                if (RS.at<float>(x, y) >= M_PI/2){
                    m->at<float>(x, y) = 0;
                } else if (RS.at<float>(x, y) <= M_PI/4){
                    m->at<float>(x, y) = 1;
                } else {
                    m->at<float>(x, y) = cos(M_PI/2 * log2(4 * RS.at<float>(i, j)/M_PI));
                }
            }
        }
        f[n] = m;
    }
    return f;
}

vector<Mat* > * SteerablePyramid::calicurate_h_filter(){
    vector<vector<Mat> > * polar = calicurate_polar();
    vector<Mat *> * f = new vector<Mat* >(n);
    for (int i = 0; i < n; i++){
        Mat* m = new Mat(size(xRes, yRes), CV_64FC1, 0);
        for (int x = 0; x < xRes; x++){
            for(int y = 0; y < yRes; y++){
                if (RS.at<float>(x, y) >= M_PI/2){
                    m->at<float>(x, y) = 1;
                } else if (RS.at<float>(x, y) <= M_PI/4){
                    m->at<float>(x, y) = 0;
                } else {
                    m->at<float>(x, y) = cos(M_PI/2 * log2(2 * RS.at<float>(i, j)/M_PI));
                }
            }
        }
        f[n] = m;
    }
    return f;
}

// On calcule les b filters un à un
Mat* SteerablePyramid::calicurate_b_filter(int i, int j){

    Mat *fil = new Mat();
    // A changer, recupere direcement la matrice AT
    vector<Mat*> * AT = calicurate_polar()->at(1);

    Mat* fil = new Mat(size(xRes, yRes), CV_64FC1, 0);
    for (int x = 0; x < xRes; x++){
        for(int y = 0; y < yRes; y++){
            float currentCoefficient = AT->at(i)->at<float>(x, y);
            if (AT->at(i)->at<float>(x, y) - j * M_PI/k < -M_PI){
                currentCoefficient += 2*M_PI;
            } else if (AT->at(i)->at<float>(x, y) - j * M_PI/k > M_PI){
                currentCoefficient -= 2*M_PI;
            } else if (abs(AT->at(i)->at<float>(x, y) - j*M_PI/k)< M_PI/2){
                currentCoefficient = alphak * pow(cos(AT->at(i)->at<float>(x, y)
                - j * M_PI/k), k-1);
            }
        }
    }

    return fil;
}

void SteerablePyramid::createPyramids(){

}

Mat* SteerablePyramid::collapsePyramids(){

}

Mat* SteerablePyramid::clearPyramids(){

}

SteerablePyramid::~SteerablePyramid(){
    this->image->release();

}
