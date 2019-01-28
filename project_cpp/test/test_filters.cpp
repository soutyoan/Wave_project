#include "SteerablePyramid.h"

int main(){
    
    Mat m = imread("../../images/lena.jpg", 0);

    SteerablePyramid s(m,
        m.rows,
        m.cols,
        3,
        2,
        "lena",
        ".",
        false,
        2,
        0.1);

    vector<Mat> RS; // calculate RS
    vector<Mat> AT; // calculate AT
    Mat filt_h(Size(m.rows, m.cols), CV_64FC1); 

    s.caliculate_polar(RS, AT);
    s.calicurate_h0_filter(filt_h, RS); 

    cout << filt_h.at<double>(511, 495) << endl;

    imwrite("../../output/high_filt.png", filt_h);

    Mat filt_u(Size(m.rows, m.cols), CV_64FC1); 

    s.calicurate_l0_filter(filt_u, RS); 

    cout << filt_u.at<double>(511, 495) << endl;

    imwrite("../../output/low_filt.png", filt_u);
    imwrite("../../output/rs.png", RS[0]);

}