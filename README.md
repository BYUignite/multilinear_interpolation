## Multilinear Interpolation in C++.

* This is a header-only code for multi-dimensional linear interpolation.
* Up to five dimensions are supported.
* Three versions are provided:
    1. multilinear_interpolation.h: this uses simple arrays, with the grid mapped to a one-D array in row-major ordering.
    2. multilinear_interpolation_vector.h: this uses STL vectors with the grid as vectors of vectors: f[i][j][k].
    3. multilinear_interpolation_pseudorecursive.h: this is an older version that effectively uses a pseudo-recursive approach.

* Grid locations are computed using the std::lower_bound function, which uses a binary search.

* There are examples and comparison to python for one to three dimensions.
* The code has been tested for three and four dimensions using another production code.

* Note that values out of bounds are linearly extrapolated using the nearest two grid points in a given dimension.
