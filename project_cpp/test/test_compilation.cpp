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
        false);

    s.collapsePyramids(); 
    
}
