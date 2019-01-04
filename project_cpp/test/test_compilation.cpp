#include "SteerablePyramid.h"

int main(){

    Mat m = imread("../../images/lena.jpg", 0); 

    SteerablePyramid s(&m, 
        m.size[0], 
        m.size[1], 
        3, 
        2, 
        "lena", 
        ".", 
        false);

    s.createPyramids(); 
    
}