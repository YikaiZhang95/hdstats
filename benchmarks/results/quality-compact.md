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

