**Huber**: median elapsed seconds of 3 runs; the ratio is time / time of hdstats' default solver at that size.

| method (seconds; ratio to hdstats default) | 200 x 1000 | 500 x 2000 | 1000 x 5000 | 2000 x 2000 |
|---|---|---|---|---|
| hdhuber() | 0.25 s (1.0x) | 2.29 s (1.0x) | 6.64 s (1.0x) | 84.94 s (1.0x) |
| hqreg(method = 'huber') | 0.11 s (0.4x) | 0.87 s (0.4x) | 2.37 s (0.4x) | 22.18 s (0.3x) |

**Quantile**: median elapsed seconds of 3 runs; the ratio is time / time of hdstats' default solver at that size.

| method (seconds; ratio to hdstats default) | 200 x 1000 | 500 x 2000 | 1000 x 5000 | 2000 x 2000 |
|---|---|---|---|---|
| hdqr() | 1.32 s (1.0x) | 6.41 s (1.0x) | 21.84 s (1.0x) | 27.07 s (1.0x) |
| hdqr(is_exact = TRUE) | 3.03 s (2.3x) | 13.19 s (2.1x) |  |  |
| hqreg(method = 'quantile') | 9.30 s (7.1x) | 26.79 s (4.2x) | 78.69 s (3.6x) | 249.01 s (9.2x) |
| conquer.reg() on hdqr's lambda sequence | 1.86 s (1.4x) | 7.30 s (1.1x) | 56.50 s (2.6x) | 24.65 s (0.9x) |
| rq(method = 'lasso'), 10 lambda values | 37.95 s (28.8x) |  |  |  |

**SVM**: median elapsed seconds of 3 runs; the ratio is time / time of hdstats' default solver at that size.

| method (seconds; ratio to hdstats default) | 200 x 1000 | 500 x 2000 | 1000 x 5000 | 2000 x 2000 |
|---|---|---|---|---|
| hdsvm() | 0.18 s (1.0x) | 0.93 s (1.0x) | 3.65 s (1.0x) | 13.70 s (1.0x) |
| hdsvm(is_exact = TRUE) | 6.92 s (37.4x) | 34.30 s (36.9x) |  |  |
| sparseSVM() | 0.10 s (0.5x) | 0.44 s (0.5x) | 2.82 s (0.8x) | 1.60 s (0.1x) |
| gcdnet(method = 'hhsvm'), Huberized squared hinge | 0.19 s (1.0x) | 0.74 s (0.8x) | 3.32 s (0.9x) | 9.15 s (0.7x) |
| LiblineaR(type = 5), one cost value | 0.06 s (0.3x) | 0.21 s (0.2x) | 0.92 s (0.3x) | 1.10 s (0.1x) |

**Rank**: median elapsed seconds of 3 runs; the ratio is time / time of hdstats' default solver at that size.

| method (seconds; ratio to hdstats default) | 100 x 10 | 200 x 20 | 400 x 20 |
|---|---|---|---|
| hdrr(), 100-value path | 0.45 s (1.0x) | 2.04 s (1.0x) | 9.02 s (1.0x) |
| rfit(), one unpenalized fit | 0.02 s (0.0x) | 0.03 s (0.0x) | 0.03 s (0.0x) |
| rq() on all pairwise differences, one fit | 0.04 s (0.1x) | 0.60 s (0.3x) | 9.97 s (1.1x) |

**Huber, n = 400, p = 200**: relative objective gap (objective - best) / best at five lambda values; df = number of nonzero coefficients; seconds = time for the five fits.

| solver | lambda = 0.158 | lambda = 0.0245 | lambda = 0.00381 | lambda = 0.000593 | lambda = 9.23e-05 | seconds |
|---|---|---|---|---|---|---|
| hdstats::hdhuber (default) | 1.9e-08 (df 5) | 3.9e-08 (df 76) | 8.5e-07 (df 176) | 1.6e-06 (df 194) | 2.1e-06 (df 199) | 0.05 |
| hdstats::hdhuber (eps = 1e-12) | 0 (df 5) | 0 (df 76) | 0 (df 176) | 0 (df 194) | 0 (df 200) | 0.10 |
| hqreg::hqreg_raw (default) | 1.2e-08 (df 5) | 3.3e-07 (df 77) | 6.3e-06 (df 176) | 3.1e-05 (df 194) | 2.5e-05 (df 200) | 0.03 |
| hqreg::hqreg_raw (eps = 1e-10) | 2.7e-12 (df 5) | 3.7e-10 (df 76) | 9.2e-09 (df 176) | 2.5e-08 (df 194) | 3.6e-08 (df 200) | 0.06 |

**Quantile (tau = 0.5), n = 400, p = 200**: relative objective gap (objective - best) / best at five lambda values; df = number of nonzero coefficients; seconds = time for the five fits.

| solver | lambda = 0.0799 | lambda = 0.0124 | lambda = 0.00193 | lambda = 0.000301 | lambda = 4.68e-05 | seconds |
|---|---|---|---|---|---|---|
| quantreg::rq lasso (LP, exact) | 0 (df 5) | 0 (df 107) | 0 (df 177) | 0 (df 197) | 0 (df 200) | 0.52 |
| hdstats::hdqr (default) | 6e-05 (df 5) | 0.0047 (df 106) | 0.011 (df 175) | 0.015 (df 198) | 0.015 (df 200) | 0.30 |
| hdstats::hdqr (is_exact = TRUE) | 5.7e-05 (df 5) | 0.0012 (df 122) | 0.0033 (df 197) | 0.0046 (df 200) | 0.0048 (df 200) | 0.32 |
| hqreg::hqreg_raw (default) | 0.00077 (df 5) | 0.01 (df 95) | 0.0071 (df 177) | 0.0029 (df 200) | 0.0032 (df 200) | 0.35 |
| hqreg::hqreg_raw (eps = 1e-10) | 0.00077 (df 5) | 0.01 (df 95) | 0.007 (df 176) | 0.0016 (df 200) | 0.0006 (df 200) | 3.94 |
| conquer::conquer.reg (Gaussian kernel, default) | 0.00037 (df 5) | 0.0089 (df 105) | 0.029 (df 189) | 0.037 (df 200) | 0.04 (df 200) | 0.10 |
| conquer::conquer.reg (uniform kernel, epsilon = 1e-6) | 0.00015 (df 5) | 0.0065 (df 99) | 0.016 (df 174) | 0.02 (df 193) | 0.021 (df 199) | 3.27 |

**Quantile (tau = 0.3), n = 400, p = 200**: relative objective gap (objective - best) / best at five lambda values; df = number of nonzero coefficients; seconds = time for the five fits.

| solver | lambda = 0.0707 | lambda = 0.011 | lambda = 0.00171 | lambda = 0.000266 | lambda = 4.14e-05 | seconds |
|---|---|---|---|---|---|---|
| quantreg::rq lasso (LP, exact) | 0 (df 6) | 0 (df 100) | 0 (df 183) | 0 (df 199) | 0 (df 200) | 0.52 |
| hdstats::hdqr (default) | 0.00013 (df 5) | 0.0047 (df 92) | 0.015 (df 184) | 0.018 (df 198) | 0.019 (df 200) | 0.31 |
| hdstats::hdqr (is_exact = TRUE) | 9.9e-05 (df 7) | 0.0017 (df 126) | 0.0047 (df 200) | 0.005 (df 200) | 0.0053 (df 200) | 0.33 |
| hqreg::hqreg_raw (default) | 0.00076 (df 5) | 0.016 (df 87) | 0.01 (df 186) | 0.0033 (df 199) | 0.003 (df 200) | 0.36 |
| hqreg::hqreg_raw (eps = 1e-10) | 0.00076 (df 5) | 0.016 (df 87) | 0.01 (df 181) | 0.0025 (df 198) | 0.00066 (df 200) | 4.12 |
| conquer::conquer.reg (Gaussian kernel, default) | 0.00038 (df 6) | 0.0093 (df 97) | 0.035 (df 195) | 0.057 (df 200) | 0.062 (df 200) | 0.10 |
| conquer::conquer.reg (uniform kernel, epsilon = 1e-6) | 0.00028 (df 6) | 0.0064 (df 91) | 0.02 (df 183) | 0.025 (df 198) | 0.026 (df 200) | 3.06 |

**SVM, n = 400, p = 200**: relative objective gap (objective - best) / best at five lambda values; df = number of nonzero coefficients; seconds = time for the five fits.

| solver | lambda = 0.111 | lambda = 0.0173 | lambda = 0.00269 | lambda = 0.000418 | lambda = 6.5e-05 | seconds |
|---|---|---|---|---|---|---|
| quantreg::rq lasso (LP, exact) | 0 (df 6) | 0 (df 119) | 0 (df 170) | 0 (df 182) | 0 (df 181) | 0.73 |
| hdstats::hdsvm (default, hval = 1) | 0.092 (df 5) | 0.51 (df 104) | 0.82 (df 172) | 1.8 (df 186) | 2.1 (df 195) | 0.64 |
| hdstats::hdsvm (hval = 0.1) | 0.18 (df 6) | 0.59 (df 121) | 1.2 (df 182) | 6.7 (df 197) | 42 (df 199) | 0.87 |
| hdstats::hdsvm (is_exact = TRUE) | 0.18 (df 6) | 0.59 (df 132) | 1.1 (df 197) | 6.2 (df 187) | 2.1 (df 195) | 1.87 |
| sparseSVM::sparseSVM (gamma = 0.1, default) | 0.0012 (df 5) | 0.43 (df 153) | 1.7 (df 196) | 16 (df 200) | 1.1e+02 (df 200) | 0.01 |
| sparseSVM::sparseSVM (gamma = 0.01, eps = 1e-8) | 0.0019 (df 6) | 0.37 (df 164) | 1.8 (df 199) | 16 (df 200) | 1.1e+02 (df 200) | 0.04 |

**Huber, n = 200, p = 500**: relative objective gap (objective - best) / best at five lambda values; df = number of nonzero coefficients; seconds = time for the five fits.

| solver | lambda = 0.31 | lambda = 0.122 | lambda = 0.0482 | lambda = 0.019 | lambda = 0.00749 | seconds |
|---|---|---|---|---|---|---|
| hdstats::hdhuber (default) | 7.2e-09 (df 2) | 5.7e-09 (df 23) | 3.2e-07 (df 111) | 4.6e-06 (df 174) | 6.7e-05 (df 190) | 0.12 |
| hdstats::hdhuber (eps = 1e-12) | 6.3e-13 (df 2) | 0 (df 23) | 0 (df 111) | 0 (df 171) | 0 (df 188) | 0.57 |
| hqreg::hqreg_raw (default) | 4e-12 (df 2) | 9.5e-09 (df 23) | 3.5e-06 (df 110) | 3.8e-05 (df 175) | 0.00066 (df 198) | 0.05 |
| hqreg::hqreg_raw (eps = 1e-10) | 0 (df 2) | 7.7e-11 (df 23) | 6.1e-09 (df 111) | 1.5e-07 (df 171) | 1.8e-06 (df 189) | 0.22 |

**Quantile (tau = 0.5), n = 200, p = 500**: relative objective gap (objective - best) / best at five lambda values; df = number of nonzero coefficients; seconds = time for the five fits.

| solver | lambda = 0.149 | lambda = 0.0586 | lambda = 0.0231 | lambda = 0.00911 | lambda = 0.00359 | seconds |
|---|---|---|---|---|---|---|
| quantreg::rq lasso (LP, exact) | 0 (df 2) | 0 (df 44) | 0 (df 148) | 0 (df 197) | 0 (df 208) | 2.88 |
| hdstats::hdqr (default) | 0.00015 (df 2) | 0.0023 (df 40) | 0.013 (df 152) | 0.037 (df 203) | 0.055 (df 219) | 0.34 |
| hdstats::hdqr (is_exact = TRUE) | 0.00015 (df 2) | 0.00063 (df 59) | 0.0021 (df 240) | 0.0047 (df 377) | 0.0064 (df 460) | 0.33 |
| hqreg::hqreg_raw (default) | 0.00046 (df 2) | 0.0047 (df 36) | 0.0089 (df 150) | 0.01 (df 241) | 0.1 (df 347) | 1.06 |
| hqreg::hqreg_raw (eps = 1e-10) | 0.00046 (df 2) | 0.0047 (df 36) | 0.0089 (df 148) | 0.0024 (df 215) | 0.041 (df 283) | 5.45 |
| conquer::conquer.reg (Gaussian kernel, default) | 0.0003 (df 2) | 0.0059 (df 40) | 0.038 (df 154) | 0.12 (df 295) | 0.23 (df 380) | 0.19 |
| conquer::conquer.reg (uniform kernel, epsilon = 1e-6) | 0.00019 (df 2) | 0.0049 (df 43) | 0.024 (df 144) | 0.061 (df 195) | 0.09 (df 204) | 5.99 |

**Quantile (tau = 0.3), n = 200, p = 500**: relative objective gap (objective - best) / best at five lambda values; df = number of nonzero coefficients; seconds = time for the five fits.

| solver | lambda = 0.133 | lambda = 0.0525 | lambda = 0.0207 | lambda = 0.00816 | lambda = 0.00322 | seconds |
|---|---|---|---|---|---|---|
| quantreg::rq lasso (LP, exact) | 0 (df 2) | 0 (df 47) | 0 (df 147) | 0 (df 193) | 0 (df 208) | 2.96 |
| hdstats::hdqr (default) | 6.3e-05 (df 2) | 0.0028 (df 43) | 0.014 (df 150) | 0.035 (df 206) | 0.048 (df 221) | 0.41 |
| hdstats::hdqr (is_exact = TRUE) | 6.3e-05 (df 2) | 0.00073 (df 63) | 0.0024 (df 228) | 0.0046 (df 384) | 0.0063 (df 462) | 0.45 |
| hqreg::hqreg_raw (default) | 0.00022 (df 2) | 0.0094 (df 39) | 0.019 (df 146) | 0.01 (df 230) | 0.078 (df 342) | 0.75 |
| hqreg::hqreg_raw (eps = 1e-10) | 0.00022 (df 2) | 0.0094 (df 39) | 0.018 (df 141) | 0.0077 (df 200) | 0.032 (df 275) | 4.86 |
| conquer::conquer.reg (Gaussian kernel, default) | 0.00011 (df 2) | 0.0066 (df 41) | 0.042 (df 158) | 0.18 (df 247) | 0.39 (df 385) | 0.19 |
| conquer::conquer.reg (uniform kernel, epsilon = 1e-6) | 8.7e-05 (df 2) | 0.0041 (df 42) | 0.029 (df 139) | 0.088 (df 190) | 0.26 (df 207) | 6.20 |

**SVM, n = 200, p = 500**: relative objective gap (objective - best) / best at five lambda values; df = number of nonzero coefficients; seconds = time for the five fits.

| solver | lambda = 0.158 | lambda = 0.0621 | lambda = 0.0245 | lambda = 0.00966 | lambda = 0.00381 | seconds |
|---|---|---|---|---|---|---|
| quantreg::rq lasso (LP, exact) | 0 (df 0) | 0 (df 89) | 0 (df 142) | 0 (df 149) | 0 (df 147) | 2.85 |
| hdstats::hdsvm (default, hval = 1) | 0.24 (df 7) | 0.67 (df 72) | 1.2 (df 114) | 1.6 (df 130) | 1.9 (df 138) | 0.05 |
| hdstats::hdsvm (hval = 0.1) | 0.24 (df 7) | 0.77 (df 88) | 1.8 (df 140) | 4.5 (df 156) | 11 (df 174) | 0.09 |
| hdstats::hdsvm (is_exact = TRUE) | 0.23 (df 0) | 0.77 (df 108) | 1.9 (df 212) | 4.8 (df 174) | 12 (df 184) | 0.26 |
| sparseSVM::sparseSVM (gamma = 0.1, default) | 0.48 (df 6) | 3.3 (df 135) | 10 (df 337) | 33 (df 444) | 1e+02 (df 448) | 0.04 |
| sparseSVM::sparseSVM (gamma = 0.01, eps = 1e-8) | 0.47 (df 6) | 2.9 (df 137) | 9.5 (df 362) | 32 (df 455) | 98 (df 458) | 1.07 |

**regression, normal** (n = 200, p = 500, s = 10, 1 replications; mean (s.e.); test error = mean absolute error on 1000 test observations).

| method | L2 error | test error | TPR | FPR | df | CV seconds |
|---|---|---|---|---|---|---|
| hdstats::hdhuber (delta = 1) | 0.629 (NA) | 0.916 (NA) | 1.000 (NA) | 0.333 (NA) | 40.0 | 0.3 |
| hqreg (huber, gamma = 1) | 0.631 (NA) | 0.916 (NA) | 1.000 (NA) | 0.333 (NA) | 40.0 | 0.1 |
| conquer (lasso, tau = 0.5) | 0.651 (NA) | 0.943 (NA) | 1.000 (NA) | 0.400 (NA) | 46.0 | 0.4 |
| glmnet (least squares lasso) | 0.676 (NA) | 0.909 (NA) | 1.000 (NA) | 0.300 (NA) | 37.0 | 0.1 |
| hqreg (quantile, tau = 0.5) | 0.693 (NA) | 0.976 (NA) | 1.000 (NA) | 0.667 (NA) | 70.0 | 5.7 |
| hdstats::hdqr (tau = 0.5) | 0.850 (NA) | 0.992 (NA) | 1.000 (NA) | 0.567 (NA) | 61.0 | 0.9 |
| hdstats::hdrr (rank) | 0.924 (NA) | 1.017 (NA) | 1.000 (NA) | 0.578 (NA) | 62.0 | 21.7 |

**regression, t3** (n = 200, p = 500, s = 10, 1 replications; mean (s.e.); test error = mean absolute error on 1000 test observations).

| method | L2 error | test error | TPR | FPR | df | CV seconds |
|---|---|---|---|---|---|---|
| conquer (lasso, tau = 0.5) | 1.186 (NA) | 1.323 (NA) | 1.000 (NA) | 0.356 (NA) | 42.0 | 0.4 |
| hdstats::hdhuber (delta = 1) | 1.196 (NA) | 1.312 (NA) | 1.000 (NA) | 0.311 (NA) | 38.0 | 0.5 |
| hqreg (huber, gamma = 1) | 1.197 (NA) | 1.312 (NA) | 1.000 (NA) | 0.311 (NA) | 38.0 | 0.1 |
| glmnet (least squares lasso) | 1.394 (NA) | 1.365 (NA) | 1.000 (NA) | 0.378 (NA) | 44.0 | 0.1 |
| hqreg (quantile, tau = 0.5) | 1.408 (NA) | 1.385 (NA) | 1.000 (NA) | 0.311 (NA) | 38.0 | 7.6 |
| hdstats::hdrr (rank) | 1.578 (NA) | 1.429 (NA) | 1.000 (NA) | 0.411 (NA) | 47.0 | 23.2 |
| hdstats::hdqr (tau = 0.5) | 1.835 (NA) | 1.528 (NA) | 1.000 (NA) | 0.367 (NA) | 43.0 | 1.2 |

**regression, contaminated** (n = 200, p = 500, s = 10, 1 replications; mean (s.e.); test error = mean absolute error on 1000 test observations).

| method | L2 error | test error | TPR | FPR | df | CV seconds |
|---|---|---|---|---|---|---|
| conquer (lasso, tau = 0.5) | 0.847 (NA) | 1.693 (NA) | 1.000 (NA) | 0.311 (NA) | 38.0 | 0.5 |
| hdstats::hdhuber (delta = 1) | 0.855 (NA) | 1.679 (NA) | 1.000 (NA) | 0.344 (NA) | 41.0 | 0.7 |
| hqreg (huber, gamma = 1) | 0.866 (NA) | 1.679 (NA) | 1.000 (NA) | 0.356 (NA) | 42.0 | 0.1 |
| hqreg (quantile, tau = 0.5) | 0.903 (NA) | 1.719 (NA) | 1.000 (NA) | 0.333 (NA) | 40.0 | 6.4 |
| hdstats::hdrr (rank) | 1.072 (NA) | 1.721 (NA) | 1.000 (NA) | 0.389 (NA) | 45.0 | 26.9 |
| hdstats::hdqr (tau = 0.5) | 1.091 (NA) | 1.750 (NA) | 1.000 (NA) | 0.411 (NA) | 47.0 | 1.1 |
| glmnet (least squares lasso) | 1.772 (NA) | 2.080 (NA) | 0.800 (NA) | 0.333 (NA) | 38.0 | 0.1 |

**classification, flip 10%** (n = 200, p = 500, s = 10, 1 replications; mean (s.e.); test error = misclassification rate on 1000 test observations).

| method | test error | TPR | FPR | df | CV seconds |
|---|---|---|---|---|---|
| glmnet (logistic lasso) | 0.320 (NA) | 1.000 (NA) | 0.278 (NA) | 35.0 | 0.2 |
| hdstats::hdsvm | 0.325 (NA) | 0.900 (NA) | 0.300 (NA) | 36.0 | 1.6 |
| hdstats::hdsvm (is_exact = TRUE) | 0.331 (NA) | 1.000 (NA) | 0.378 (NA) | 44.0 | 20.2 |
| gcdnet (Huberized squared hinge) | 0.347 (NA) | 1.000 (NA) | 0.878 (NA) | 89.0 | 1.0 |
| sparseSVM | 0.420 (NA) | 1.000 (NA) | 0.944 (NA) | 95.0 | 0.0 |

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

