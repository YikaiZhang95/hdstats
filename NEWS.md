# hdststs 0.1.0

* Initial release.
* Unifies the `hdhuber`, `hdsvm`, and `hdqr` solvers into a single package.
* Coordinate-descent kernels re-implemented in C++ via Rcpp (translated from the
  original Fortran), with the Huber kernel restructured to enable SIMD
  vectorization (~3x faster than the Fortran version).
* User-facing API unchanged from the original packages: `hdhuber()`, `hdsvm()`,
  `hdqr()`, their `cv.*` cross-validation wrappers, `nc.hdsvm()`/`nc.hdqr()`
  non-convex (SCAD/MCP) variants, and `predict()`/`coef()` methods.
