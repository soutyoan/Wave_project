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

vector<float>* calicurate_h0_filter(){
    vector<vector<float> > * polar = calicurate_polar();
    vector<float> RS = polar->at(0);
    vector<float> * fil = new vector<float>();
    // Possible parallelisation?
    for (int i = 0; i < RS.size(); i++){
        if (RS[i] >= M_PI/2){
            fil->push_back(1);
        } else if (RS[i] < 0){
            fil->push_back(0);
        } else {
            fil->push_back(cos(M_PI/2 * log2(RS[i]/M_PI)));
        }
    }
    return fil;
}

vector<float>* calicurate_l0_filter(){
    vector<vector<float> > * polar = calicurate_polar();
    vector<float> RS = polar->at(0);
    vector<float> * fil = new vector<float>();
    // Possible parallelisation
    for (int i = 0; i < RS.size(); i++){
        if (RS[i] <= M_PI/2){
            fil->push_back(1);
        } else if (RS[i] >= 0){
            fil->push_back(0);
        } else {
            fil->push_back(cos(M_PI/2 * log2(2  *RS[i]/M_PI)));
        }
    }
    return fil;
}

vector<float>* calicurate_l_filter(){

}

vector<float>* calicurate_h_filter(){

}

vector<float>* calicurate_b_filter(){

}

void createPyramids(){

}

vector<float>* collapsePyramids(){

}

vector<float>* clearPyramids(){

}

SteerablePyramid::~SteerablePyramid(){
    this->image->release();

}
