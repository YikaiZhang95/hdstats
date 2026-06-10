## ---------------------------------------------------------------------------
## Reproduce the worked examples from Section 5 of the hdstats paper, verbatim,
## and check that the printed outputs match the manuscript.
## Expected: 5.1 misclassification error == 0 ; 5.2 selected vars qr=115 huber=80
## ---------------------------------------------------------------------------
library("hdstats")

## 5.1 -- Simulated high-dimensional classification --------------------------
set.seed(1)
n <- 200; p <- 1000
x <- matrix(rnorm(n * p), n, p)
beta_true <- c(rep(2, 10), rep(0, p - 10))
y <- ifelse(drop(x %*% beta_true) + rnorm(n) > 0, 1, -1)
cv_svm <- cv.hdsvm(x, y, nfolds = 5)
pred <- predict(cv_svm, newx = x, s = "lambda.min", type = "class")
err <- mean(pred != y)
stopifnot(err == 0)
cat("5.1 training misclassification error:", err, " (paper: 0)  OK\n")

## 5.2 -- Simulated high-dimensional regression ------------------------------
set.seed(2)
n <- 200; p <- 800
x <- matrix(rnorm(n * p), n, p)
beta_true <- c(runif(15, 1, 2), rep(0, p - 15))
y <- drop(x %*% beta_true + rt(n, df = 3))
fit_qr <- hdqr(x, y, tau = 0.5, lam2 = 0.01)
fit_huber <- hdhuber(x, y, delta = 1, lam2 = 0.01)
sel <- c(qr = sum(coef(fit_qr, s = fit_qr$lambda[30])[-1] != 0),
         huber = sum(coef(fit_huber, s = fit_huber$lambda[30])[-1] != 0))
stopifnot(sel["qr"] == 115, sel["huber"] == 80)
cat("5.2 selected variables:  qr =", sel["qr"], " huber =", sel["huber"],
    " (paper: 115, 80)  OK\n")

cat("\nAll paper examples reproduced exactly under",
    R.version.string, "with hdstats", as.character(packageVersion("hdstats")), "\n")
