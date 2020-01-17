#include <iostream>
#include <vector>
#include "multilinear_interpolation.h"
//#include "multilinear_interpolation_vector.h"
//#include "multilinear_interpolation_pseudorecursive.h"

using namespace std;

/////////////////////////////////////////

int main(){

    vector<double> x{1,2,5,9};
    vector<double> f{1,2,3,4};

    double xP = 7.0;

    double fP = LI_1D(4, &x[0], &f[0], xP);   // multilinear_interpolation.h
    //double fP = LI_1D(x, f, xP);            // multilinear_interpolation_vector and _pseudorecursive.h

    cout << endl << fP << endl;

    return 0;
}

