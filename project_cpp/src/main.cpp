#include "SteerablePyramid.h"

int main(int argc, char *argv[]){
    Mat image = imread("../../images/lena.jpg", CV_LOAD_IMAGE_COLOR);
    SteerablePyramid h(
        image,
        image.rows, 
        image.cols, 
        5,
        3,
        "lena",
        "lena.png",
        false
    );
    return 0;
}
