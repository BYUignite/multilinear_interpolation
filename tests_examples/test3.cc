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
    vector<double> z{1,2,5,9,10,11};
    vector<vector<vector<double>>> f(4, vector<vector<double>>(5, vector<double>(6)));
    vector<double> ff(4*5*6);

    for(int i=0; i<x.size(); ++i)
        for(int j=0; j<y.size(); ++j)
            for(int k=0; k<z.size(); ++k){
                f[i][j][k] = sqrt(i+j+k);
                ff[i*6*5+j*6+k] = sqrt(i+j+k);
            }

    double xP = 7.0;
    double yP = 2.2;
    double zP = 5.5;

    double fP = LI_3D(4, 5, 6, &x[0], &y[0], &z[0], &ff[0], xP, yP, zP);   // multilinear_interpolation.h
    //double fP = LI_3D(x, y, z, f, xP, yP, zP);                        // multilinear_interpolation_vector.h, _pseudorecursive.h

    cout << endl << fP << endl;

    return 0;
}

