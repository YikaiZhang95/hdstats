test_that("hdrr is median regression on pairwise differences", {
  d <- sim_reg(n = 30, p = 10)
  fit <- hdrr(d$x, d$y, lam2 = 0.01, nlambda = 10)
  expect_s3_class(fit, "hdrr")
  expect_s3_class(fit, "hdqr")
  expect_equal(dim(fit$beta), c(10, 10))
  ## same slopes as hdqr(tau = 0.5) on the differenced data
  pr <- which(upper.tri(matrix(0, 30, 30)), arr.ind = TRUE)
  xd <- d$x[pr[, 1], ] - d$x[pr[, 2], ]
  yd <- d$y[pr[, 1]] - d$y[pr[, 2]]
  ref <- hdqr(xd, yd, tau = 0.5, lam2 = 0.01, nlambda = 10)
  expect_equal(fit$lambda, ref$lambda)
  expect_equal(as.matrix(fit$beta), as.matrix(ref$beta))
  ## a different ordering of the pairs changes only rounding
  pr2 <- t(combn(30, 2))
  ref2 <- hdqr(d$x[pr2[, 1], ] - d$x[pr2[, 2], ], d$y[pr2[, 1]] - d$y[pr2[, 2]],
               tau = 0.5, lam2 = 0.01, nlambda = 10)
  expect_equal(as.matrix(fit$beta), as.matrix(ref2$beta), tolerance = 1e-4)
  ## intercept is the median residual
  res <- d$y - d$x %*% as.matrix(fit$beta)
  expect_equal(unname(fit$b0), unname(apply(res, 2, median)))
  ## coef/predict methods of hdqr apply
  cf <- coef(fit, s = fit$lambda[5])
  expect_equal(dim(cf), c(11, 1))
  expect_equal(predict(fit, d$x[1:3, ], s = fit$lambda[5]),
               as.matrix(cbind(1, d$x[1:3, ]) %*% cf), ignore_attr = TRUE)
  expect_error(hdrr(d$x[1, , drop = FALSE], d$y[1]), "at least 2")
})

test_that("unpenalized hdrr minimizes Jaeckel's Wilcoxon dispersion", {
  d <- sim_reg(n = 40, p = 3, seed = 5)
  fit <- hdrr(d$x, d$y, lambda = c(1, 0), lam2 = 0, standardize = FALSE,
              eps = 1e-12)
  b <- as.vector(fit$beta[, 2])
  pr <- t(combn(40, 2))
  disp <- function(beta) {
    r <- d$y - drop(d$x %*% beta)
    sum(abs(r[pr[, 1]] - r[pr[, 2]]))
  }
  ref <- optim(c(0, 0, 0), disp, method = "Nelder-Mead",
               control = list(reltol = 1e-12, maxit = 5000))
  expect_lt(disp(b) - ref$value, 1e-3 * ref$value)
  expect_equal(b, ref$par, tolerance = 0.05)
})

test_that("cv.hdrr works and inherits the cv.hdqr methods", {
  d <- sim_reg(n = 30, p = 8)
  set.seed(1)
  cv <- cv.hdrr(d$x, d$y, lam2 = 0.01, nfolds = 3, nlambda = 10)
  expect_s3_class(cv, "cv.hdrr")
  expect_s3_class(cv, "cv.hdqr")
  expect_true(cv$lambda.min %in% cv$lambda)
  expect_true(cv$lambda.1se >= cv$lambda.min)
  expect_equal(dim(coef(cv, s = "lambda.min")), c(9, 1))
  expect_length(predict(cv, d$x[1:4, ], s = "lambda.1se"), 4)
  expect_output(print(cv), "Call:")
  pdf(NULL); on.exit(dev.off())
  expect_silent(plot(cv))
})
