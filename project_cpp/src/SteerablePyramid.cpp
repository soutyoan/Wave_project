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
                bool _verbose,
                double C,
                double beta){

    if ( _image.empty() ){
        cerr << "Image is empty" << endl;
        exit(1);
    }
    this->image = _image;
    this->xRes = _xRes;
    this->yRes = _yRes;
    this->n = _n;
    this->k = _k;
    this->image_name = _image_name;
    this->output_path = _output_path;
    this->verbose = _verbose;
    this->C = C;
    this->beta = beta;

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
void SteerablePyramid::caliculate_one_polar(Mat &RS, Mat &AT, int i){
    double _tmp = pow(2.0, i);
    size_t nx = this->image.rows / _tmp;
    size_t ny = this->image.cols / _tmp;
    vector<double> _wx = linspace<double>(-M_PI, M_PI, nx);
    vector<double> _wy = linspace<double>(-M_PI, M_PI, ny);
    cout << _wx.size() << " " << _wy.size() << " " << nx << " " << ny << "\n";
    polar_coordinates<double>(_wx, _wy, nx, ny, RS);
    angular_coordinates<double>(_wx, _wy, nx, ny, AT);
}

void SteerablePyramid::caliculate_polar(vector<Mat> &RS, vector<Mat> &AT){

    for (int i=0; i < this->n; i++) {
        // Computation polar coordinates (radius) on the grid
        double _tmp = pow(2.0, i);
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
    // Possible parallelisation?
    for (int i = 0; i < xRes; i++){
        for (int j = 0; j < yRes; j++){
            if (RS_0.at<double>(i, j) >= M_PI){
                fil.at<double>(i, j) = 1;
            } else if (RS_0.at<double>(i, j) >= M_PI/2) {
                fil.at<double>(i, j) = cos(M_PI/2 * log2(RS_0.at<double>(i, j)/M_PI));
            }
        }
    }
}

void SteerablePyramid::calicurate_l0_filter(Mat &fil, vector<Mat> &RS){
    Mat RS_0 = RS[0];
    // Possible parallelisation?
    for (int i = 0; i < xRes; i++){
        for (int j = 0; j < yRes; j++){
            if (RS_0.at<double>(i, j) >= M_PI){
                fil.at<double>(i, j) = 0;
            } else if (RS_0.at<double>(i, j) <= M_PI/2){
                fil.at<double>(i, j) = 1;
            } else {
                fil.at<double>(i, j) = cos(M_PI/2 * log2(2 * RS_0.at<double>(i, j)/M_PI));
            }
        }
    }
}

void SteerablePyramid::calicurate_l_filter(int i, Mat &f, vector<Mat> &RS){
    Mat RS_0 = RS[0];
    double _tmp = pow(2.0, i);
    for (int x = 0; x < xRes/_tmp; x++){
        for(int y = 0; y < yRes/_tmp; y++){
            if (RS_0.at<double>(x, y) >= M_PI/2){
                f.at<double>(x, y) = 0;
            } else if (RS_0.at<double>(x, y) <= M_PI/4){
                f.at<double>(x, y) = 1;
            } else {
                f.at<double>(x, y) = cos(M_PI/2 * log2(4 * RS_0.at<double>(x, y)/M_PI));
            }
        }
    }
}

void SteerablePyramid::calicurate_h_filter(vector<Mat> &f, vector<Mat> &RS){
    Mat RS_0 = RS[0];
    for (int i = 0; i < n; i++){
        double _tmp = pow(2.0, i);
        Mat m(Size(xRes/_tmp, yRes/_tmp), CV_64FC1, 0);
        for (int x = 0; x < xRes/_tmp; x++){
            for(int y = 0; y < yRes/_tmp; y++){
                if (RS_0.at<double>(x, y) >= M_PI/2){
                    m.at<double>(x, y) = 1;
                } else if (RS_0.at<double>(x, y) <= M_PI/4){
                    m.at<double>(x, y) = 0;
                } else {
                    m.at<double>(x, y) = cos(M_PI/2 * log2(2 * RS_0.at<double>(x, y)/M_PI));
                }
            }
        }
    }
}

// On calcule les b filters un à un
void SteerablePyramid::calicurate_b_filter(int i, int j, Mat &fil, const vector<Mat> &AT){

    double _tmp = pow(2.0, i);

    for (int x = 0; x < xRes/_tmp; x++){
        for(int y = 0; y < yRes/_tmp; y++){
            double currentCoefficient = AT[i].at<double>(x, y);
            if (AT[i].at<double>(x, y) - j * M_PI/k < -M_PI){
                currentCoefficient += 2*M_PI;
            } else if (AT[i].at<double>(x, y) - j * M_PI/k > M_PI){
                currentCoefficient -= 2*M_PI;
            } else if (abs(AT[i].at<double>(x, y) - j*M_PI/k)< M_PI/2){
                currentCoefficient = alphak * pow(cos(AT[i].at<double>(x, y)
                - j * M_PI/k), k-1);
            }
            fil.at<double>(x, y) = currentCoefficient;
        }
    }
}

void fillDownfillMatrix(Mat &down, int x, int y){
    for (int i = x; i < 3 * x; i++){
        for (int j = y; j < 3 * y; j++){
            down.at<double>(i, j) = 1;
        }
    }
}

Mat SteerablePyramid::createPyramids(vector<Mat> &RS, vector<Mat> &AT, vector<Mat> &BND){

    // Create all the matrices used during the collapse function
    Mat h0f(Size(xRes, yRes), CV_64FC1);
    Mat h0s(Size(xRes, yRes), CV_64FC1);
    Mat l0f(Size(xRes, yRes), CV_64FC1);
    Mat l0s(Size(xRes, yRes), CV_64FC1);

    // Find library for fast fourier transform
    Mat ft = Mat_<std::complex<double> > (Size(xRes, yRes), CV_64FC1);
    Mat _ft = Mat_<std::complex<double> > (Size(xRes, yRes), CV_64FC1);

    ///// THREAD 1 in openmp use thread ID //////

    Mat imgBack(Size(xRes, yRes), CV_64FC1);
    image.convertTo(image, CV_64FC1);

    dft((image), ft);
    ft_shift<std::complex<double> >(ft, _ft);

    // Create HO filter
    Mat h0(Size(xRes, yRes), CV_64FC1);
    calicurate_h0_filter(h0, RS);

    Mat m = mul_complex(h0, _ft);  // calculation of h0
    Mat f_ishift = Mat_<std::complex<double> > (Size(xRes, yRes), CV_64FC1);
    ft_shift<std::complex<double> >(m, f_ishift); // FFT opencv
    idft(f_ishift, imgBack); // IFFT opencv

    //// THREAD 2 in openmp use thread ID /////

    Mat l0(Size(xRes, yRes), CV_64FC1);
    calicurate_l0_filter(l0, RS);

    _ft = mul_complex(l0, _ft);  // calculation of h0

    ft_shift<std::complex<double> >(_ft, f_ishift); // FFT opencv
    idft(f_ishift, imgBack); // IFFT opencv

    // Mat lastImage = l0;

    double _tmp = 1;

    Mat lastImage = Mat_<std::complex<double> >(xRes/_tmp, yRes/_tmp, CV_64FC1);

    // PROBLEME FLOAT ET COMPLEXE

    // Parallelize here (use openmp for)
    for (int i = 0; i < n; i ++){

        _tmp = pow(2, i);

        for (int j = 0; j < k; j++){

            Mat b_filter(Size(xRes/_tmp, yRes/_tmp), CV_64FC1);
            calicurate_b_filter(i, j, b_filter, AT);

            Mat lb = mul_complex(b_filter, lastImage);

            Mat f_ishift = Mat_<std::complex<double> >(Size(xRes/_tmp, yRes/_tmp), CV_64FC1); // Shift lb
            Mat img_back = Mat_<std::complex<double> >(Size(xRes/_tmp, yRes/_tmp), CV_64FC1); // ishift2 f_ishift

            ft_shift<std::complex<double> >(lb, img_back);
            idft(lb, f_ishift);

            BND.push_back(lb);
            BND.push_back(img_back);

        }

        cout << xRes/_tmp << endl;

        // Apply low pas filter to image downsampled
        int img_x = lastImage.size[0];
        int img_y = lastImage.size[1];
        int quantification_x = (int)(img_x/4);
        int quantification_y = (int)(img_y/4);

        Mat_<std::complex<double> > down_fil(Size(img_x, img_y), CV_64FC1);
        fillDownfillMatrix(down_fil, quantification_x, quantification_y);

        Mat_<std::complex<double> > l(Size(img_x, img_y), CV_64FC1);
        calicurate_l_filter(i, l, RS);

        cout << "OK\n";

        l = l.mul(down_fil);

        cout << "ICI\n";

        down_fil = mul_complex(l, lastImage);

        Mat down_image = Mat_<std::complex<double> >(2 * quantification_x, 2 * quantification_y);
        for (int i = quantification_x; i < 3 * quantification_x; i ++){
            for (int j = quantification_y; j < 3* quantification_y; j++){
                down_image.at<complex<double> >(i - quantification_x,
                    j - quantification_y) = down_fil.at<complex<double> >(i, j);
            }
        }

        // Par rapport à la version python, on ne sauvegarde pas low.
        // Ca ne semble pas avoir d'interet.

        lastImage.create(2 * quantification_x, 2 * quantification_y, CV_64FC1);

        lastImage = down_image;
    }

    Mat LRf = lastImage;

    return LRf;

}

Mat SteerablePyramid::collapsePyramids(){

    vector<Mat> RS; // calculate RS
    vector<Mat> AT; // calculate AT
    vector<Mat> BND; // calculate AT

    this->caliculate_polar(RS, AT);

    //Vecteur retour de la valeur de BND
    Mat resid = createPyramids(RS, AT, BND);

    cout << "Pyramid creation OK" << endl;

    for (int i = 0; i < BND.size(); i++){
        cout << i << " " << BND[i].rows << " " << BND[i].cols << endl;
    }

    for (int i = n-1; i > -1; i--){

        Mat tmp = Mat_<std::complex<double> > (Size(2*resid.rows, 2 * resid.rows), CV_64FC1);

        int quant4x = (int)(resid.rows/2);
        int quant4y = (int)(resid.cols/2);

        // Recopy the matrix
        for (int i = quant4x; i < 3 * quant4x; i++){
            for (int j = quant4y; j < 3 * quant4y; j++){
                tmp.at<complex<double> >(i, j) =
                    resid.at<complex<double> >(i - quant4x, j - quant4y);
            }
        }

        Mat filt(Size(2*resid.rows, 2 * resid.rows), CV_64FC1);
        calicurate_l_filter(i, filt, AT);

        tmp = mul_complex(filt, tmp);

        for (int j = 0; j < k; j++){
            Mat filt = Mat_<std::complex<double> >(Size(2*resid.rows, 2 * resid.rows), CV_64FC1);
            calicurate_b_filter(i, j, filt, AT);
            cout << BND[2*(i * k + j)].rows << " " << BND[2*(i * k + j)].cols << endl;
            cout << filt.rows << " " << filt.cols << endl;
            cout << tmp.rows << " " << tmp.cols << endl;
            tmp = tmp + BND[2*(i * k + j)].mul(filt);
        }

        resid.create(2*resid.rows, 2*resid.cols, CV_64FC1);
        resid = tmp;
    }

    Mat h0(Size(xRes, yRes), CV_64FC1);
    calicurate_h0_filter(h0, RS);

    Mat l0(Size(xRes, yRes), CV_64FC1);
    calicurate_l0_filter(l0, RS);

    Mat ft = Mat_<std::complex<double> >(Size(xRes, yRes), CV_64FC1);
    Mat _ft = Mat_<std::complex<double> >(Size(xRes, yRes), CV_64FC1);

    Mat imgBack(Size(xRes, yRes), CV_64FC1);
    image.convertTo(image, CV_64FC1);

    dft((image), ft);
    ft_shift<std::complex<double> >(ft, _ft);

    Mat recon(Size(xRes, yRes), CV_64FC1);
    Mat mul = mul_complex(h0, _ft);
    recon = mul_complex(l0, resid) + mul_complex(h0, mul);

    Mat recon_real(Size(xRes, yRes), CV_64FC1);
    idft(recon, recon_real);

    cout << recon_real.channels() << endl;

    Mat ch1;
    // "channels" is a vector of 3 Mat arrays:
    vector<Mat> channels(2);
    // split img:
    split(recon_real, channels);
    // get the channels (dont forget they follow BGR order in OpenCV)
    ch1 = channels[0];

    return ch1;

}

/**
 * Computation of w weights for LPC strength computation
 * @param w resulting vector of double w
 */
void SteerablePyramid::get_w(vector<double>& w)
{
    // TODO
}

/**
 * Computation of Local Phase Coherence strength
 * @param  w vector of weights
 * @param  j chosen orientation
 * @param  k spatial location
 * @return   S_{LPC}^{j,k}
 */
double SteerablePyramid::LPCStrength(vector<double> w, int j, int k)
{
    // TODO
}

/**
 * Computation of Local Phase Coherence coefficient
 * @param  k spatial location
 * @return   S_{LPC}^{k}
 */
double SteerablePyramid::SpatialLPC(int k)
{
    // TODO
}

/**
 * Computation of LPC map
 * @param map output LPC maps
 */
void SteerablePyramid::LPCMap(Mat& map)
{
    // TODO
}

/**
 * Computation of the sharpness index
 * @return LPC-based sharpness index
 */
double SteerablePyramid::LPCSharpnessIndex()
{
    // TODO
}

void SteerablePyramid::clearPyramids(Mat &f){

}

SteerablePyramid::~SteerablePyramid(){

}
