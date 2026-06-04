# hdstats

<!-- badges: start -->
[![R-CMD-check](https://github.com/YikaiZhang95/hdstats/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/YikaiZhang95/hdstats/actions/workflows/R-CMD-check.yaml)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- badges: end -->

**High-dimensional smoothed Huber, SVM, and quantile regression in one package.**

`hdstats` unifies three penalized high-dimensional regression solvers — Huber
regression, support vector machines, and quantile regression — behind a single
R interface. All three use the *finite smoothing algorithm* (uniform-density
convolution smoothing) with coordinate descent, and an elastic-net penalty
(weighted L1 + L2). The compute kernels are implemented in **C++ via Rcpp**.

This package combines and re-implements the previously separate
[`hdhuber`](https://CRAN.R-project.org/package=hdhuber),
[`hdsvm`](https://CRAN.R-project.org/package=hdsvm), and
[`hdqr`](https://CRAN.R-project.org/package=hdqr) Fortran packages. The numerical
results are equivalent to the originals (see [Correctness](#correctness)).

## Installation

```r
# install.packages("remotes")
remotes::install_github("YikaiZhang95/hdstats")
```

A C++ compiler is required (Rtools on Windows, Xcode CLT on macOS, `r-base-dev`
on Linux).

## Models

| Function    | Loss              | Response       | Key extra argument |
|-------------|-------------------|----------------|--------------------|
| `hdhuber()` | Huber             | continuous     | `delta` (transition point) |
| `hdsvm()`   | hinge (SVM)       | binary `{-1,1}`| `hval` (smoothing) |
| `hdqr()`    | quantile / check  | continuous     | `tau` (quantile), `hval` |
| `hdrr()`    | Wilcoxon rank     | continuous     | — (built on `hdqr`) |

Each model ships with:

- `cv.*()` — k-fold cross-validation
- `predict()` / `coef()` — S3 methods
- `nc.hdsvm()` / `nc.hdqr()` — non-convex (SCAD / MCP) penalties via the local
  linear approximation
- `hdrr()` — canonical Wilcoxon rank regression, solved as median quantile
  regression on pairwise differences (reuses the `hdqr` `coef`/`predict` methods)

## Usage

### Huber regression

```r
library(hdstats)
set.seed(1)
n <- 100; p <- 200
x <- matrix(rnorm(n * p), n, p)
y <- x[, 1] * 2 - x[, 2] * 1.5 + rnorm(n)

fit <- hdhuber(x, y, lam2 = 0.01, delta = 1)
cvfit <- cv.hdhuber(x, y, delta = 1)
coef(cvfit, s = cvfit$lambda.min)
predict(fit, newx = x[1:5, ], s = fit$lambda[10])
```

### Support vector machine

```r
prob <- plogis(c(x %*% rnorm(p) * 0.1))
yb <- 2 * rbinom(n, 1, prob) - 1

fit <- hdsvm(x, yb, lam2 = 0.01)
cvfit <- cv.hdsvm(x, yb)
predict(fit, newx = x[1:5, ], s = cvfit$lambda.min)   # class labels

# non-convex SCAD penalty
ncfit <- nc.hdsvm(x, yb, lambda = fit$lambda[1:10], lam2 = 0.01, pen = "scad")
```

### Quantile regression

```r
fit <- hdqr(x, y, tau = 0.5, lam2 = 0.01)
cvfit <- cv.hdqr(x, y, tau = 0.5)
coef(cvfit, s = cvfit$lambda.1se)

# non-convex MCP penalty
ncfit <- nc.hdqr(x, y, tau = 0.5, lambda = fit$lambda[1:10], lam2 = 0.01, pen = "mcp")
```

### Wilcoxon rank regression

Canonical (Wilcoxon) rank regression is the median quantile-regression objective
on all pairwise differences of the data, so `hdrr()` builds the differenced
design and calls `hdqr()` internally. The returned object inherits the `hdqr`
`coef`/`predict` methods.

```r
fit <- hdrr(x, y, lam2 = 0.01)
coef(fit, s = fit$lambda[5])
predict(fit, newx = x[1:5, ], s = fit$lambda[5])
```

## Performance

The C++ kernels are benchmarked against the original Fortran implementations
(both compiled at `-O2` on the same machine; 100-lambda paths). Speedup = Fortran
time / C++ time.

| Model     | Speedup vs Fortran            | Notes |
|-----------|-------------------------------|-------|
| `hdhuber` | **2.6–3.4×** (median 2.9×)    | removes a loop-carried dependency that blocked SIMD vectorization |
| `hdsvm`   | ~parity (0.9–1.0×)            | Fortran inner loop was already well vectorized |
| `hdqr`    | up to **1.7×**, grows with *n*| faster as the sample size increases |

Solutions are numerically identical to the Fortran packages: Huber coefficients
agree to ~1e-6; SVM/QR objectives agree to ~1e-6 (these solvers converge to a KKT
tolerance, so coefficients agree at that level rather than to machine precision).

## Correctness

The package is validated against the original Fortran packages on shared inputs
and on the simulation designs from the reference paper (Wang et al. 2006 for SVM;
Friedman et al. 2010 for quantile regression). Lambda paths are identical and
objective values match to within solver tolerance.

## Reference

The finite smoothing algorithm is described in:

> *Finite smoothing algorithm for high-dimensional support vector machines and quantile regression*.

## License

GPL (>= 2). See [LICENSE.md](LICENSE.md).
