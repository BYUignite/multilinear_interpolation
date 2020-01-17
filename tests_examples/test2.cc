#include <iostream>
#include <vector>
#include <cmath>
#include "multilinear_interpolation.h"
//#include "multilinear_interpolation_vector.h"
//#include "multilinear_interpolation_pseudorecursive.h"

using namespace std;

/////////////////////////////////////////

int main(){

    vector<double> x{1,2,5,9};
    vector<double> y{1,2,5,9,10};
    vector<vector<double>> f(4, vector<double>(5));
    vector<double> ff(4*5);

    for(int i=0; i<x.size(); ++i)
        for(int j=0; j<y.size(); ++j){
            f[i][j] = sqrt(i+j);
            ff[i*5+j] = sqrt(i+j);
        }

    double xP = 7.0;
    double yP = 2.2;

    double fP = LI_2D(4, 5, &x[0], &y[0], &ff[0], xP, yP);   // multilinear_interpolation.h
    //double fP = LI_2D(x, y, f, xP, yP);                        // multilinear_interpolation_vector.h, _pseudorecursive.h

    cout << endl << fP << endl;

    return 0;
}

