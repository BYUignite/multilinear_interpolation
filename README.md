## Multilinear Interpolation in C++.

This is a header-only code for multi-dimensional linear interpolation.
Code supports up to five dimensions.
Three versions are provided:
1. multilinear_interpolation.h: this uses simple arrays, with the grid mapped to a one-D array in row-major ordering.
2. multilinear_interpolation_vector.h: this uses STL vectors with the grid as vectors of vectors: f[i][j][k].
3. multilinear_interpolation_pseudorecursive.h: this is an older version that effectively uses a pseudo-recursive approach.
Grid locations are computed using the std::lower_bound function, which uses a binary search.

