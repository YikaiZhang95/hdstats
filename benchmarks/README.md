# Benchmarks: hdstats versus mainstream packages

Reproducible comparison of **hdstats** with the mainstream R packages that fit
the same or closely related models. All scripts run from the package root and
write their results to `results/` and figures to `figures/`:

| Script | What it measures |
|---|---|
| `01-speed.R` | Time to compute a 100-value solution path on identical data, single-threaded, median of 3 runs. |
| `02-quality.R` | Optimization accuracy: the value of the penalized objective reached by each solver at the same `lambda` values (no standardization), relative to the exact optimum. |
| `03-accuracy.R` | Statistical performance at the cross-validation-selected `lambda` of each package: estimation error, variable selection, test error, over replicated simulations. |
| `04-figures.R`, `05-tables.R` | Figures and markdown tables from the CSV results. |

Run with `HDSTATS_BENCH_QUICK=1` for a smoke test; `HDSTATS_BENCH_REPS` sets
the number of timing repetitions (speed) or simulation replications (accuracy).

## Competitors and how the problems were matched

| Model | hdstats | Competitors | Notes on matching |
|---|---|---|---|
| Huber regression | `hdhuber()` | **hqreg** (semismooth Newton coordinate descent) | hqreg's Huber loss is ours divided by `gamma`, so `lambda_hqreg = lambda / delta` reproduces the same problem. |
| Quantile regression | `hdqr()` | **hqreg** (`method = "quantile"`), **conquer** (convolution-smoothed QR, LAMM), **quantreg** (`rq(method = "lasso")`, interior-point LP) | Same check loss and `(1/n) sum + lambda * L1` scaling in hqreg and conquer. quantreg's lasso penalizes the summed loss with `lambda/2 * L1`, so `lambda_rq = 2 n lambda`. |
| Linear SVM | `hdsvm()` | **sparseSVM** (huberized hinge, semismooth Newton CD), **gcdnet** (`hhsvm`, Huberized *squared* hinge), **LiblineaR** (`type = 5`, L1-regularized *squared* hinge, one cost value), **quantreg** LP (hinge = check loss at `tau = 1`) | Only sparseSVM targets the same hinge-loss objective; gcdnet and LiblineaR fit different losses and are reported for speed only. sparseSVM codes the first level of `y` as +1, so its weights are sign-flipped before evaluating the objective. |
| Wilcoxon rank regression | `hdrr()` | **Rfit** (unpenalized), **quantreg** on all pairwise differences (unpenalized LP) | No mainstream package fits penalized rank regression; `hdrr()` computes a 100-value path while the competitors compute one unpenalized fit. |
| Least squares (reference) | | **glmnet** | Reference for the statistical comparison. |

`rqPen` was not included: its `alg = "huber"` is hqreg and its `alg = "br"` is
quantreg, both covered above.

## Results

### 1. Speed

Single-threaded elapsed time for a 100-value `lambda` path on the same
simulated data (sparse signal with 10 nonzero coefficients, AR(1) design,
t3 errors; for the SVM a linear signal with 10% label noise), median of 3
runs. `lambda.min / lambda.max` is 0.01 for p > n and 1e-4 for n >= p in all
packages that expose it. quantreg's linear-programming solver has no warm
starts and is run on 10 `lambda` values only, for the smallest size.

![speed](figures/speed.png)

**Huber**: median elapsed seconds of 3 runs; the ratio is time / time of hdstats' default solver at that size.

| method (seconds; ratio to hdstats default) | 200 x 1000 | 500 x 2000 | 1000 x 5000 | 2000 x 2000 |
|---|---|---|---|---|
| hdhuber() | 0.09 s (1.0x) | 0.80 s (1.0x) | 2.43 s (1.0x) | 34.31 s (1.0x) |
| hqreg(method = 'huber') | 0.11 s (1.3x) | 0.90 s (1.1x) | 2.54 s (1.0x) | 22.44 s (0.7x) |

**Quantile**: median elapsed seconds of 3 runs; the ratio is time / time of hdstats' default solver at that size.

| method (seconds; ratio to hdstats default) | 200 x 1000 | 500 x 2000 | 1000 x 5000 | 2000 x 2000 |
|---|---|---|---|---|
| hdqr() | 0.60 s (1.0x) | 3.11 s (1.0x) | 10.65 s (1.0x) | 15.87 s (1.0x) |
| hdqr(is_exact = TRUE) | 1.23 s (2.1x) | 5.44 s (1.7x) |  |  |
| hqreg(method = 'quantile') | 9.28 s (15.6x) | 26.99 s (8.7x) | 80.72 s (7.6x) | 259.00 s (16.3x) |
| conquer.reg() on hdqr's lambda sequence | 1.85 s (3.1x) | 7.16 s (2.3x) | 57.94 s (5.4x) | 27.50 s (1.7x) |
| rq(method = 'lasso'), 10 lambda values | 38.00 s (63.8x) |  |  |  |

**SVM**: median elapsed seconds of 3 runs; the ratio is time / time of hdstats' default solver at that size.

| method (seconds; ratio to hdstats default) | 200 x 1000 | 500 x 2000 | 1000 x 5000 | 2000 x 2000 |
|---|---|---|---|---|
| hdsvm() | 0.10 s (1.0x) | 0.41 s (1.0x) | 1.78 s (1.0x) | 6.75 s (1.0x) |
| hdsvm(is_exact = TRUE) | 2.69 s (27.1x) | 10.36 s (25.0x) |  |  |
| sparseSVM() | 0.09 s (0.9x) | 0.37 s (0.9x) | 2.11 s (1.2x) | 2.13 s (0.3x) |
| gcdnet(method = 'hhsvm'), Huberized squared hinge | 0.15 s (1.5x) | 0.69 s (1.7x) | 3.02 s (1.7x) | 5.79 s (0.9x) |
| LiblineaR(type = 5), one cost value | 0.06 s (0.6x) | 0.17 s (0.4x) | 0.79 s (0.4x) | 0.90 s (0.1x) |

**Rank**: median elapsed seconds of 3 runs; the ratio is time / time of hdstats' default solver at that size.

| method (seconds; ratio to hdstats default) | 100 x 10 | 200 x 20 | 400 x 20 |
|---|---|---|---|
| hdrr(), 100-value path | 0.26 s (1.0x) | 1.17 s (1.0x) | 5.01 s (1.0x) |
| rfit(), one unpenalized fit | 0.02 s (0.1x) | 0.03 s (0.0x) | 0.03 s (0.0x) |
| rq() on all pairwise differences, one fit | 0.04 s (0.2x) | 0.63 s (0.5x) | 10.06 s (2.0x) |


What the timings show:

- **Huber regression.** `hdhuber()` is 1.05 to 1.3 times faster than
  hqreg's semismooth Newton coordinate descent when p > n and 1.5 times
  slower at n = p = 2000, where the path ends in the nearly unpenalized,
  ill-conditioned regime (lambda.min = 1e-4 lambda.max) in which the
  majorization steps of the finite smoothing algorithm are most damped.
- **Quantile regression.** `hdqr()` is 7.6 to 16 times faster than hqreg's
  quantile solver and 1.7 to 5.4 times faster than conquer, while reaching
  a better objective than either (next section). The exact LP solver in
  quantreg needs 38 s for 10 `lambda` values at the smallest size, against
  0.6 s for the 100-value `hdqr()` path. `is_exact = TRUE` roughly doubles
  the run time.
- **SVM.** sparseSVM and `hdsvm()` are on par when p > n (0.9 to 1.2 times
  `hdsvm()`'s time); sparseSVM is 3 times faster at n = p = 2000, where its
  advantage grows with n. gcdnet's Huberized squared hinge takes 1.5 to
  1.7 times as long for p > n, and a single LiblineaR fit costs 13 to 60%
  of a whole `hdsvm()` path. `hdsvm(is_exact = TRUE)` is about 25 times
  slower than the default.
- **Rank regression.** No mainstream package fits penalized rank
  regression. A 100-value `hdrr()` path costs 0.3 to 5 s for n = 100 to
  400, against 0.02 to 0.03 s for one unpenalized `Rfit::rfit()` fit; the
  exact LP on all pairwise differences (quantreg) takes twice as long as
  the whole `hdrr()` path at n = 400.

These timings are for the current C++ kernels, which keep the derivative
of the smoothed loss for every observation and use unrolled, branch-free
inner loops; they are 1.4 to 2.8 times faster than the first C++ version
and 2 to 2.8 (Huber, quantile) or 1.4 to 1.6 (SVM) times faster than the
Fortran originals on the same machine.

**A Newton-type coordinate step was evaluated and not adopted.** The
majorization step of the finite smoothing algorithm uses the curvature
bound mean(x_k^2) for every coordinate; hqreg's semismooth Newton step uses
instead the curvature of the smoothed loss along the coordinate (the mean
of x_ik^2 over the residuals inside the band), with a continuation term when
few residuals are in band. We implemented this step in the Huber kernel
together with a safeguard that replaces a step by the majorization step
whenever the directional derivative of the coordinate objective changes
sign (this restores monotone descent; without it the step cycles on some
of the package's own test problems). On the speed-benchmark data with
delta = 1 the Newton-type step reduced the number of coordinate-descent
passes by 10-20% only, because most residuals lie inside the band once the
fit is good, and by about 45% with delta = 0.1; the extra pass that
computes the curvature costs about as much as the whole majorization
update, so the wall time increased in every setting (500 x 2000, t3
errors: 2.1 s against 0.8 s). The majorization step is therefore kept.

### 2. Optimization quality

For the same `lambda` values (five positions along the 100-value path) and
with no standardization, each solver's objective value is compared with the
exact optimum: the Huber problem is smooth and all solvers converge to the
same value; the quantile and hinge problems are linear programs solved
exactly with `quantreg::rq(method = "lasso")`. The full per-`lambda` tables
(with the number of nonzero coefficients and the time for the five fits)
are in [`results/tables.md`](results/tables.md).

Maximum relative objective gap over the five lambda values, by solver (0 = the best solution found; the LP references are exact).

| model | solver | max gap, n = 200, p = 500 | max gap, n = 400, p = 200 |
|---|---|---|---|
| Huber | hdstats::hdhuber (eps = 1e-12) | 6.3e-13 |   0 |
| Huber | hqreg::hqreg_raw (eps = 1e-10) | 1.8e-06 | 3.6e-08 |
| Huber | hdstats::hdhuber (default) | 6.7e-05 | 2.1e-06 |
| Huber | hqreg::hqreg_raw (default) | 0.00066 | 3.1e-05 |
| QR tau = 0.3 | quantreg::rq lasso (LP, exact) |   0 |   0 |
| QR tau = 0.3 | hdstats::hdqr (is_exact = TRUE) | 0.0063 | 0.0053 |
| QR tau = 0.3 | hqreg::hqreg_raw (eps = 1e-10) | 0.032 | 0.016 |
| QR tau = 0.3 | hdstats::hdqr (default) | 0.048 | 0.019 |
| QR tau = 0.3 | hqreg::hqreg_raw (default) | 0.078 | 0.016 |
| QR tau = 0.3 | conquer::conquer.reg (uniform kernel, epsilon = 1e-6) | 0.26 | 0.026 |
| QR tau = 0.3 | conquer::conquer.reg (Gaussian kernel, default) | 0.39 | 0.062 |
| QR tau = 0.5 | quantreg::rq lasso (LP, exact) |   0 |   0 |
| QR tau = 0.5 | hdstats::hdqr (is_exact = TRUE) | 0.0064 | 0.0048 |
| QR tau = 0.5 | hqreg::hqreg_raw (eps = 1e-10) | 0.041 | 0.01 |
| QR tau = 0.5 | hdstats::hdqr (default) | 0.055 | 0.015 |
| QR tau = 0.5 | conquer::conquer.reg (uniform kernel, epsilon = 1e-6) | 0.09 | 0.021 |
| QR tau = 0.5 | hqreg::hqreg_raw (default) | 0.1 | 0.01 |
| QR tau = 0.5 | conquer::conquer.reg (Gaussian kernel, default) | 0.23 | 0.04 |
| SVM | quantreg::rq lasso (LP, exact) |   0 |   0 |
| SVM | hdstats::hdsvm (default, hval = 1) | 1.9 | 2.1 |
| SVM | hdstats::hdsvm (hval = 0.1) |  11 |  42 |
| SVM | hdstats::hdsvm (is_exact = TRUE) |  12 | 6.2 |
| SVM | sparseSVM::sparseSVM (gamma = 0.01, eps = 1e-8) |  98 | 1.1e+02 |
| SVM | sparseSVM::sparseSVM (gamma = 0.1, default) | 1e+02 | 1.1e+02 |


What the objective gaps show:

- **Huber.** `hdhuber()` and hqreg both reach the optimum; at default
  tolerances `hdhuber()` is within 7e-5 and hqreg within 7e-4, and both go
  to machine precision with tighter tolerances.
- **Quantile regression.** With default settings `hdqr()` is within 2%
  (n > p) and 5.5% (p > n) of the LP optimum, `is_exact = TRUE` within
  0.65%. hqreg's defaults are within 1.6% (n > p) and 10% (p > n), and
  conquer within 6% (n > p) and 39% (p > n). conquer targets the
  convolution-smoothed estimator with a statistically chosen bandwidth,
  so its gap to the non-smoothed optimum is by design; `hdqr()` and hqreg
  target the non-smoothed problem.
- **SVM.** No coordinate-descent solver comes close to the hinge-loss
  optimum on these problems. With the default bandwidth (`hval = 1`)
  `hdsvm()` is only solving the h = 1 smoothed problem, and its hinge
  objective is 1.1 to 3 times the optimum (gap 0.09 to 2.1). Decreasing the
  bandwidth (`hval = 0.1`) or using `is_exact = TRUE` makes the objective
  *worse* (gap up to 42), because coordinate descent on the smoothed hinge
  slows down as the bandwidth shrinks and stops far from the optimum even
  with 1e6 passes. sparseSVM (huberized hinge) is further away still
  (gap 0.4 to 110). With the ridge penalty alone (`lambda = 0`),
  `hdsvm(is_exact = TRUE)` does reproduce the SVM quadratic-program
  solution of e1071 (see `tests/testthat/test-hdsvm.R`). Cross-validated
  prediction accuracy is much less sensitive than the objective value
  (section 3), but the paper should not claim exact solutions of the
  L1-penalized SVM.

### 3. Statistical performance

Each package selects its own `lambda` by 5-fold cross-validation on the
same folds (`lambda.min`); the CV time includes the full-data fit.
Regression: n = 200, p = 500, 10 nonzero coefficients, AR(1) design with
correlation 0.5, 10 replications, errors N(0,1), t3, or 90% N(0,1) + 10%
N(0,100); test error is the mean absolute error on 1000 new observations.
Classification: the same design with two blocks of five coefficients of
opposite sign, Gaussian noise before taking the sign, and 5% of the labels
flipped; test error is the misclassification rate on 1000 new observations.
`hdqr()` and `hdrr()` default to an elastic-net ridge term `lam2 = 0.01`; the
`lam2 = 0` rows use the lasso penalty like the competitors.

![regression](figures/accuracy-regression.png)

![classification](figures/accuracy-classification.png)

**classification, flip 5%** (n = 200, p = 500, s = 10, 10 replications; mean (s.e.); test error = misclassification rate on 1000 test observations).

| method | test error | TPR | FPR | df | CV seconds |
|---|---|---|---|---|---|
| hdstats::hdsvm | 0.193 (0.007) | 0.810 (0.035) | 0.032 (0.006) | 23.9 | 0.6 |
| glmnet (logistic lasso) | 0.197 (0.005) | 0.830 (0.026) | 0.045 (0.007) | 30.2 | 0.1 |
| hdstats::hdsvm (is_exact = TRUE) | 0.216 (0.008) | 0.790 (0.055) | 0.074 (0.029) | 44.0 | 18.7 |
| sparseSVM | 0.221 (0.006) | 0.830 (0.026) | 0.049 (0.013) | 32.3 | 0.2 |
| gcdnet (Huberized squared hinge) | 0.242 (0.010) | 0.860 (0.037) | 0.162 (0.028) | 88.0 | 0.5 |

**regression, normal** (n = 200, p = 500, s = 10, 10 replications; mean (s.e.); test error = mean absolute error on 1000 test observations).

| method | L2 error | test error | TPR | FPR | df | CV seconds |
|---|---|---|---|---|---|---|
| glmnet (least squares lasso) | 1.136 (0.048) | 1.056 (0.021) | 1.000 (0.000) | 0.159 (0.012) | 87.8 | 0.1 |
| hqreg (huber, gamma = 1) | 1.140 (0.052) | 1.065 (0.021) | 1.000 (0.000) | 0.165 (0.015) | 91.0 | 0.3 |
| hdstats::hdhuber (delta = 1) | 1.147 (0.051) | 1.065 (0.021) | 1.000 (0.000) | 0.160 (0.015) | 88.6 | 0.8 |
| conquer (lasso, tau = 0.5) | 1.257 (0.070) | 1.099 (0.024) | 1.000 (0.000) | 0.155 (0.012) | 86.1 | 2.0 |
| hdstats::hdqr (tau = 0.5, lam2 = 0) | 1.270 (0.043) | 1.147 (0.030) | 1.000 (0.000) | 0.214 (0.028) | 115.0 | 2.6 |
| hqreg (quantile, tau = 0.5) | 1.277 (0.052) | 1.157 (0.028) | 1.000 (0.000) | 0.286 (0.037) | 150.1 | 26.5 |
| hdstats::hdqr (tau = 0.5) | 1.558 (0.055) | 1.238 (0.034) | 1.000 (0.000) | 0.207 (0.029) | 111.5 | 4.4 |
| hdstats::hdrr (rank) | 1.627 (0.077) | 1.222 (0.034) | 1.000 (0.000) | 0.154 (0.030) | 85.3 | 136.3 |

**regression, t3** (n = 200, p = 500, s = 10, 10 replications; mean (s.e.); test error = mean absolute error on 1000 test observations).

| method | L2 error | test error | TPR | FPR | df | CV seconds |
|---|---|---|---|---|---|---|
| hdstats::hdhuber (delta = 1) | 1.660 (0.119) | 1.472 (0.040) | 0.950 (0.022) | 0.118 (0.010) | 67.4 | 1.4 |
| hqreg (huber, gamma = 1) | 1.663 (0.118) | 1.472 (0.040) | 0.950 (0.022) | 0.117 (0.009) | 66.9 | 0.4 |
| hqreg (quantile, tau = 0.5) | 1.727 (0.132) | 1.524 (0.046) | 0.920 (0.039) | 0.179 (0.014) | 97.0 | 34.8 |
| hdstats::hdqr (tau = 0.5, lam2 = 0) | 1.730 (0.140) | 1.515 (0.046) | 0.910 (0.038) | 0.139 (0.011) | 77.4 | 4.2 |
| conquer (lasso, tau = 0.5) | 1.773 (0.148) | 1.515 (0.047) | 0.890 (0.043) | 0.117 (0.013) | 66.0 | 2.2 |
| glmnet (least squares lasso) | 1.819 (0.127) | 1.522 (0.034) | 0.940 (0.031) | 0.114 (0.012) | 65.1 | 0.1 |
| hdstats::hdqr (tau = 0.5) | 1.987 (0.112) | 1.603 (0.045) | 0.910 (0.031) | 0.150 (0.013) | 82.7 | 5.0 |
| hdstats::hdrr (rank) | 2.075 (0.109) | 1.597 (0.037) | 0.910 (0.038) | 0.106 (0.011) | 60.9 | 176.3 |

**regression, contaminated** (n = 200, p = 500, s = 10, 10 replications; mean (s.e.); test error = mean absolute error on 1000 test observations).

| method | L2 error | test error | TPR | FPR | df | CV seconds |
|---|---|---|---|---|---|---|
| hdstats::hdqr (tau = 0.5, lam2 = 0) | 1.799 (0.123) | 1.973 (0.046) | 0.910 (0.035) | 0.122 (0.012) | 69.0 | 9.3 |
| hqreg (huber, gamma = 1) | 1.799 (0.126) | 1.954 (0.041) | 0.910 (0.046) | 0.098 (0.013) | 57.0 | 0.9 |
| hdstats::hdhuber (delta = 1) | 1.807 (0.125) | 1.956 (0.041) | 0.920 (0.039) | 0.095 (0.012) | 55.7 | 3.9 |
| hqreg (quantile, tau = 0.5) | 1.868 (0.151) | 2.004 (0.055) | 0.920 (0.047) | 0.151 (0.020) | 83.0 | 47.5 |
| conquer (lasso, tau = 0.5) | 1.938 (0.149) | 2.009 (0.056) | 0.880 (0.053) | 0.098 (0.013) | 56.7 | 2.6 |
| hdstats::hdqr (tau = 0.5) | 2.114 (0.130) | 2.092 (0.046) | 0.880 (0.053) | 0.133 (0.020) | 74.2 | 6.4 |
| hdstats::hdrr (rank) | 2.416 (0.139) | 2.182 (0.036) | 0.830 (0.052) | 0.091 (0.019) | 52.9 | 269.1 |
| glmnet (least squares lasso) | 3.175 (0.068) | 2.505 (0.024) | 0.370 (0.052) | 0.038 (0.008) | 22.2 | 0.1 |


What the statistical comparison shows:

- **Huber regression.** `hdhuber()` and hqreg's Huber solver give the same
  estimator and the same numbers: on par with the least-squares lasso under
  normal errors and clearly better under t3 errors (L2 error 1.66 versus
  1.82) and under contamination (1.81 versus 3.18, where the least-squares
  lasso recovers only 37% of the true variables).
- **Quantile regression.** With `lam2 = 0`, `hdqr()` matches hqreg's
  quantile solver under normal and t3 errors and is slightly better under
  contamination (1.80 versus 1.87), at 5 to 18 times lower CV time; it is
  on par with conquer. With the default `lam2 = 0.01` its estimation error
  is 15 to 25% higher, because the ridge term shrinks the large true
  coefficients and moves the CV choice to a larger `lambda`. The same
  default applies to `hdrr()`, whose results here were obtained with it
  (`HDSTATS_BENCH_PART=extra2` reruns the rank regression with `lam2 = 0`;
  it takes about two hours and was not completed for this write-up).
- **Rank regression.** `hdrr()` is robust (it degrades far less than the
  least-squares lasso under contamination) but is the least accurate of the
  robust methods here and by far the most expensive, since its design has
  n(n - 1)/2 rows.
- **SVM.** At the CV-selected `lambda`, `hdsvm()` has the lowest test error
  (0.193) and the lowest false-positive rate (0.032, 24 selected variables),
  ahead of the logistic lasso (0.197), sparseSVM (0.221) and gcdnet (0.242);
  `is_exact = TRUE` is worse (0.216) and 30 times slower. So the large
  optimization gaps of section 2 do not translate into worse predictions
  in this setting, which is consistent with the h = 1 smoothed hinge being
  a sensible loss in its own right.

### 4. Platform

**Platform**:  R: R version 4.3.3 (2024-02-29); platform: x86_64-pc-linux-gnu; os: Linux 6.18.44-fc-v33; cpu: Intel(R) Xeon(R) Processor @ 2.80GHz; BLAS: libblas.so.3.12.0; compiler: g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0 

**Package versions**

| package | version |
|---|---|
| hdstats | 0.2.0 |
| hqreg | 1.4.1 |
| sparseSVM | 1.1.7 |
| conquer | 1.3.3 |
| quantreg | 5.97 |
| gcdnet | 1.0.6 |
| LiblineaR | 2.10.25 |
| Rfit | 0.27.0 |
| glmnet | 4.1.8 |


## Implications for the JSS manuscript

1. **Lead with quantile regression and robustness; the speed story is
   now solid.** `hdqr()` is 8 to 16 times faster than hqreg's quantile
   solver and 2 to 5 times faster than conquer while reaching a
   closer-to-exact objective than either; `hdhuber()` matches or beats
   hqreg for p > n, and `hdsvm()` matches sparseSVM for p > n. The one
   regime where the majorization steps lag (n = p, nearly unpenalized end
   of the path) is worth stating rather than hiding.
2. **Reconsider the default `lam2 = 0.01` of `hdqr()` and `hdrr()`** (a
   default inherited from the original packages): it costs 15 to 25% in
   estimation error in these sparse settings and is inconsistent with
   `hdhuber()` and `hdsvm()`, whose default is `lam2 = 0`.
3. **Do not describe the L1-penalized SVM solutions as exact.** Both modes
   return smoothed-problem solutions whose hinge objective can be several
   times the LP optimum; the `is_exact` projection step does not fix this
   and coordinate descent stalls for small bandwidths. If exactness is a
   claim of the paper, the algorithm needs a different termination
   criterion or a final LP/active-set polish. sparseSVM has the same
   problem, so a comparison on objective values against an LP reference is
   a fair and interesting table for the paper.
4. **Report optimization accuracy alongside speed** (the LP-reference
   protocol used here, with the `lambda_rq = 2 n lambda` conversion for
   quantreg's lasso) so that reviewers cannot attribute the speed to loose
   convergence.
