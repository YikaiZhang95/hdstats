# hdstats

<!-- badges: start -->
[![R-CMD-check](https://github.com/YikaiZhang95/hdstats/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YikaiZhang95/hdstats/actions/workflows/R-CMD-check.yaml)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- badges: end -->

**High-dimensional penalized Huber, quantile, and rank regression and
support vector machines in one package.**

`hdstats` fits the entire regularization path of four elastic-net penalized
models for high-dimensional data behind a single, `glmnet`-like interface:

| Function    | Loss                    | Response         | Model-specific argument |
|-------------|-------------------------|------------------|-------------------------|
| `hdhuber()` | Huber                   | continuous       | `delta` (transition point) |
| `hdqr()`    | quantile (check) loss   | continuous       | `tau` (quantile level) |
| `hdrr()`    | Wilcoxon rank (Jaeckel) | continuous       | none (built on `hdqr()`) |
| `hdsvm()`   | hinge (linear SVM)      | binary           | none |

All solvers use coordinate descent with warm starts and active-set or
strong-rule screening. The non-smooth hinge and check losses are handled by
the *finite smoothing algorithm* (uniform-kernel convolution smoothing with
a decreasing bandwidth) of Tang, Zhang, and Wang (2024). The compute
kernels are written in C++ via `Rcpp`.

Each model comes with

- `cv.*()` for k-fold cross-validation, with `plot()` and `print()` methods;
- `coef()` and `predict()` methods that interpolate between `lambda` values;
- `nc.hdsvm()` / `nc.hdqr()` (and `cv.nc.*()`) for the nonconvex SCAD and
  MCP penalties via the local linear approximation;
- adaptive penalty factors (`pf`, `pf2`), variable exclusion, and
  path-truncation controls (`dfmax`, `pmax`).

This package combines and re-implements the previously separate
[`hdqr`](https://CRAN.R-project.org/package=hdqr) and
[`hdsvm`](https://CRAN.R-project.org/package=hdsvm) CRAN packages and the
[`hdhuber`](https://github.com/YikaiZhang95/hdhuber) package, whose Fortran
kernels were translated to C++ (see [Correctness](#correctness)).

## Installation

```r
# install.packages("remotes")
remotes::install_github("YikaiZhang95/hdstats")
```

A C++ compiler is required (Rtools on Windows, Xcode Command Line Tools on
macOS, `r-base-dev` on Linux).

## Usage

```r
library(hdstats)
set.seed(1)
n <- 100; p <- 200
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1] * 2 - x[, 2] * 1.5 + rt(n, df = 3)
```

### Huber regression

```r
fit <- hdhuber(x, y, delta = 1, lam2 = 0.01)
fit                                  # df and lambda along the path
cvfit <- cv.hdhuber(x, y, delta = 1, lam2 = 0.01)
plot(cvfit)
coef(cvfit, s = "lambda.min")
predict(fit, newx = x[1:5, ], s = 0.1)
```

### Quantile regression

```r
fit <- hdqr(x, y, tau = 0.5, lam2 = 0.01)
cvfit <- cv.hdqr(x, y, tau = 0.5, lam2 = 0.01)
coef(cvfit, s = "lambda.1se")

# nonconvex MCP penalty (lambda sequence taken from the lasso fit if omitted)
ncfit <- nc.hdqr(x, y, tau = 0.5, lam2 = 0.01, pen = "mcp")
```

### Wilcoxon rank regression

Canonical (Wilcoxon) rank regression minimizes Jaeckel's dispersion of the
residuals, which equals the median check loss of all pairwise differences of
the data. `hdrr()` builds the differenced design and calls `hdqr()`; the
returned object inherits the `hdqr` methods.

```r
fit <- hdrr(x, y, lam2 = 0.01)
cvfit <- cv.hdrr(x, y, lam2 = 0.01)
coef(cvfit, s = "lambda.min")
```

### Support vector machine

```r
yb <- factor(ifelse(plogis(x[, 1] - x[, 2]) > runif(n), "case", "control"))
fit <- hdsvm(x, yb, lam2 = 0.01)
cvfit <- cv.hdsvm(x, yb, lam2 = 0.01)
predict(cvfit, newx = x[1:5, ], s = "lambda.min")   # class labels (-1/1)

# nonconvex SCAD penalty
ncfit <- nc.hdsvm(x, yb, lam2 = 0.01, pen = "scad")
```

See the vignette (`vignette("hdstats")`) for a longer tour.

## Correctness

The C++ kernels are validated against the original Fortran implementations
on identical inputs (`lambda` sequences, coefficients, intercepts, and
objective values along 50-value paths, for several `n`/`p` configurations):

- `hdhuber()`: identical algorithm; coefficients agree to about 1e-8.
- `hdqr()` and `hdsvm()`: both implementations converge to a KKT tolerance,
  so coefficients agree to that tolerance (1e-4 to 1e-2) while the values
  of the penalized objective agree to 1e-7 to 5e-5.

The package also ships a `testthat` suite that checks the KKT conditions of
the Huber solutions, compares unpenalized and lasso-penalized quantile
regression with `quantreg::rq()`, compares `hdsvm(is_exact = TRUE)` with the
linear SVM of `e1071::svm()`, and checks rank regression against a direct
minimization of Jaeckel's dispersion.

## Performance

Run times depend on the compiler and platform. The benchmark in the
[`hdhuber`](https://github.com/YikaiZhang95/hdhuber) repository reports a
median 2.5x speedup of the restructured C++ Huber kernel over the Fortran
original (100-lambda paths, `-O2`). On Linux with `gcc`/`gfortran` 13 at
`-O2` and `-O3`, however, we measured parity: 0.85-0.9x for `hdhuber()` and
1.0-1.2x for `hdqr()` and `hdsvm()` (n = 200-500, p = 200-500). The
`is_exact = TRUE` option of `hdsvm()` and `hdqr()` adds a projection step
and is correspondingly slower.

## Reference

The finite smoothing algorithm is described in

> Tang, Q., Zhang, Y., and Wang, B. (2024). Finite smoothing algorithm for
> high-dimensional support vector machines and quantile regression.
> *Proceedings of the 41st International Conference on Machine Learning*.
> <https://openreview.net/forum?id=RvwMTDYTOb>

To cite the package, see `citation("hdstats")`.

## License

GPL (>= 2). See [LICENSE.md](LICENSE.md).
