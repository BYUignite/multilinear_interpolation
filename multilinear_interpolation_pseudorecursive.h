/* 
 Multilinear interpolation.
 Up to five dimensions.
 These are hard coded. This is relatively easy to read, but not as general as, 
 e.g., https://github.com/parsiad/mlinterp
 This is an older version that effectively uses a pseudo-recursive approach. stl vectors are used.
 Pass in the grid values for each dimension,
    then the vector of function values at each point,
    then the locations of the desired points in each dimension.
 See LI_2D for specific description of values passed.
 LI_3D and LI_4D are tested by comparison of the rcSLW code to the python version.
*/

#include <vector>
#include <tuple>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

tuple<int, int> get_bounding_points_uniform(const vector<double> &x, const double xP){
    double dx = x[1] - x[0];
    int ilo = int((xP-x[0])/dx);
    if (ilo < 0)
        ilo = 0;
    else if (ilo >= x.size() - 1)
        ilo = x.size() - 2;
    int ihi = ilo + 1;
    return make_tuple(ilo, ihi);
}
//-------------------------------------------------------------------------------

tuple<int, int> get_bounding_points(const vector<double> &x, const double xP){
    int ilo, ihi;
    if(xP <= x[0])
        ihi = 1;
    else if(xP >= x.back())
        ihi = x.size()-1;
    else {
        vector<const double>::iterator itHi = lower_bound(x.begin(), x.end(), xP); // lower_bound gives values >= xP
        ihi = itHi - x.begin();
    }
    ilo = ihi-1;
    return make_tuple(ilo, ihi);
}

////////////////////////////////////////////////////////////////////////////////

double LI_1D(const vector<double> &x,
             const vector<double> &f, 
             const double xP){

    int ilo, ihi;
    tie(ilo, ihi) = get_bounding_points(x, xP);

    return f[ilo] + (xP - x[ilo]) * (f[ihi]-f[ilo])/(x[ihi]-x[ilo]);

}

////////////////////////////////////////////////////////////////////////////////

double LI_2D(const vector<double> &x,                   // x grid values
             const vector<double> &y,                   // y grid values
             const vector<vector<double> > &f,          // field of function values f[i][j] with i-->x and j-->y
             const double xP, const double yP){         // interpolation point in x
                                                        // interpolation point in y
    int ilo, ihi;
    tie(ilo, ihi) = get_bounding_points(x, xP);

    ////////////// interpolate grid to yP ////////////////////////

    double flo = LI_1D(y, f[ilo], yP);
    double fhi = LI_1D(y, f[ihi], yP);

    ////////////// interpolate final x direction /////////////////

    return LI_1D( vector<double> {x[ilo], x[ihi]}, vector<double> {flo, fhi}, xP);
}

////////////////////////////////////////////////////////////////////////////////

double LI_3D(const vector<double> &x, 
             const vector<double> &y, 
             const vector<double> &z, 
             const vector<vector<vector<double> > > &f, 
             const double xP, const double yP, const double zP){

    int ilo, ihi;
    tie(ilo, ihi) = get_bounding_points(x, xP);

    ///////////////////// interpolate grid to yP, zP /////////////////

    double flo = LI_2D(y, z, f[ilo], yP, zP);
    double fhi = LI_2D(y, z, f[ihi], yP, zP);

    //////////////////// interpolate final x direction //////////////

    return LI_1D(vector<double> {x[ilo], x[ihi]}, vector<double> {flo, fhi}, xP);

}

////////////////////////////////////////////////////////////////////////////////

double LI_4D(const vector<double> &x, 
             const vector<double> &y, 
             const vector<double> &z, 
             const vector<double> &w, 
             const vector<vector<vector<vector<double> > > > &f, 
             const double xP, const double yP, const double zP, const double wP){

    int ilo, ihi;
    tie(ilo, ihi) = get_bounding_points(x, xP);

    ////////////////// interpolate grid to yP, zP ////////////////////

    double flo = LI_3D(y, z, w, f[ilo], yP, zP, wP);
    double fhi = LI_3D(y, z, w, f[ihi], yP, zP, wP);
    
    //////////////// interpolate final x direction //////////////////

    return LI_1D( vector<double> {x[ilo], x[ihi]}, vector<double> {flo, fhi}, xP);

}

////////////////////////////////////////////////////////////////////////////////

double LI_5D(const vector<double> &x, 
             const vector<double> &y, 
             const vector<double> &z, 
             const vector<double> &w, 
             const vector<double> &a, 
             const vector<vector<vector<vector<vector<double> > > > > &f, 
             const double xP, const double yP, const double zP, const double wP, const double aP){

    int ilo, ihi;
    tie(ilo, ihi) = get_bounding_points(x, xP);


    ////////////////// interpolate grid to yP, zP ////////////////////

    double flo = LI_4D(y, z, w, a,  f[ilo], yP, zP, wP, aP);
    double fhi = LI_4D(y, z, w, a,  f[ihi], yP, zP, wP, aP);
    
    //////////////// interpolate final x direction //////////////////

    return LI_1D( vector<double> {x[ilo], x[ihi]}, vector<double> {flo, fhi}, xP);

}

