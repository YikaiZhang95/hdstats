test_that("hdhuber returns a well-formed solution path", {
  d <- sim_reg()
  fit <- hdhuber(d$x, d$y, nlambda = 20, lam2 = 0.01)
  expect_s3_class(fit, "hdhuber")
  expect_s4_class(fit$beta, "dgCMatrix")
  expect_equal(dim(fit$beta), c(ncol(d$x), 20))
  expect_length(fit$b0, 20)
  expect_length(fit$lambda, 20)
  expect_true(all(diff(fit$lambda) < 0))
  expect_equal(unname(fit$df), unname(colSums(as.matrix(fit$beta) != 0)))
  expect_equal(fit$jerr, 0)
  expect_equal(fit$delta, 1)
  expect_equal(fit$lam2, 0.01)
  ## the largest lambda gives the null model
  expect_equal(fit$df[1], 0)
  ## the first variables to enter the path are true signal variables
  expect_true(all(which(as.matrix(fit$beta)[, 3] != 0) %in% 1:3))
})

test_that("hdhuber satisfies the KKT conditions of the penalized problem", {
  d <- sim_reg(n = 80, p = 20)
  delta <- 0.8; lam2 <- 0.05
  pf2 <- runif(20, 0.5, 2)
  fit <- hdhuber(d$x, d$y, delta = delta, lam2 = lam2, pf2 = pf2,
                 standardize = FALSE, nlambda = 15, eps = 1e-12)
  for (l in c(5, 10, 15)) {
    cf <- path_coef(fit, l)
    lambda <- fit$lambda[l]
    r <- d$y - cf$b0 - drop(d$x %*% cf$beta)
    psi <- pmin(pmax(r, -delta), delta)
    grad <- -colMeans(d$x * psi) + lam2 * pf2 * cf$beta
    active <- cf$beta != 0
    expect_lt(abs(mean(psi)), 1e-6)
    expect_lt(max(abs(grad[active] + lambda * sign(cf$beta[active]))), 1e-6)
    if (any(!active)) expect_lt(max(abs(grad[!active])), lambda + 1e-6)
  }
})

test_that("unpenalized hdhuber agrees with a general-purpose optimizer", {
  d <- sim_reg(n = 120, p = 4)
  fit <- hdhuber(d$x, d$y, lambda = c(1, 0), lam2 = 0, delta = 0.5,
                 standardize = FALSE, eps = 1e-12)
  cf <- path_coef(fit, 2)
  obj <- function(par) huber_obj(par[1], par[-1], d$x, d$y, 0.5, 0, 0)
  ref <- optim(c(0, rep(0, 4)), obj, method = "BFGS",
               control = list(reltol = 1e-14, maxit = 1000))
  expect_equal(c(cf$b0, cf$beta), ref$par, tolerance = 1e-4)
  expect_lte(obj(c(cf$b0, cf$beta)), ref$value + 1e-8)
})

test_that("the strong rule does not change the solution", {
  d <- sim_reg(n = 50, p = 60)
  f1 <- hdhuber(d$x, d$y, nlambda = 20, lam2 = 0.01, is_strong = TRUE, eps = 1e-12)
  f2 <- hdhuber(d$x, d$y, nlambda = 20, lam2 = 0.01, is_strong = FALSE, eps = 1e-12)
  expect_equal(f1$lambda, f2$lambda)
  expect_equal(as.matrix(f1$beta), as.matrix(f2$beta), tolerance = 1e-4)
  expect_equal(f1$b0, f2$b0, tolerance = 1e-4)
})

test_that("standardization is undone in the returned coefficients", {
  d <- sim_reg(n = 80, p = 10)
  xs <- scale(d$x)
  ## with standardized inputs, standardize = TRUE and FALSE solve the same
  ## problem (up to the 1/n vs 1/(n-1) scaling of the sd), so paths match
  attr(xs, "scaled:center") <- NULL; attr(xs, "scaled:scale") <- NULL
  xs <- xs * sqrt(nrow(xs) / (nrow(xs) - 1))
  f1 <- hdhuber(xs, d$y, nlambda = 10, standardize = TRUE)
  f2 <- hdhuber(xs, d$y, nlambda = 10, standardize = FALSE)
  expect_equal(f1$lambda, f2$lambda, tolerance = 1e-8)
  expect_equal(as.matrix(f1$beta), as.matrix(f2$beta), tolerance = 1e-5)
})

test_that("penalty factors, exclusion, and user lambda are honoured", {
  d <- sim_reg(n = 60, p = 15)
  ## exclude
  fit <- hdhuber(d$x, d$y, nlambda = 10, exclude = c(2, 5))
  expect_true(all(as.matrix(fit$beta)[c(2, 5), ] == 0))
  ## pf = 0 means never penalized: variable is in every model
  pf <- rep(1, 15); pf[4] <- 0
  fit <- hdhuber(d$x, d$y, nlambda = 10, pf = pf)
  expect_true(all(as.matrix(fit$beta)[4, ] != 0))
  ## user-supplied lambda is sorted decreasingly and used as is
  lam <- c(0.05, 0.5, 0.2)
  fit <- hdhuber(d$x, d$y, lambda = lam)
  expect_equal(fit$lambda, sort(lam, decreasing = TRUE))
  ## a pf matrix with one column per lambda is accepted
  pfm <- matrix(1, 15, 3)
  expect_silent(hdhuber(d$x, d$y, lambda = lam, pf = pfm))
  expect_error(hdhuber(d$x, d$y, lambda = lam, pf = matrix(1, 15, 2)),
               "pf should be a matrix")
})

test_that("input vectors are not modified in place", {
  d <- sim_reg(n = 40, p = 8)
  pf2 <- c(-1, rep(1, 7)); pf2_copy <- pf2
  pf <- c(-0.5, rep(1, 7)); pf_copy <- pf
  fit <- hdhuber(d$x, d$y, nlambda = 5, pf = pf, pf2 = pf2)
  expect_identical(pf2, pf2_copy)
  expect_identical(pf, pf_copy)
  ## negative factors are treated as zero (no penalty on that variable)
  expect_true(all(as.matrix(fit$beta)[1, ] != 0))
})

test_that("hdhuber validates its inputs", {
  d <- sim_reg(n = 30, p = 5)
  expect_error(hdhuber(d$x, d$y[-1]), "different number of observations")
  expect_error(hdhuber(d$x, d$y, delta = -1), "delta")
  expect_error(hdhuber(d$x, d$y, lambda.factor = 1), "lambda.factor")
  expect_error(hdhuber(d$x, d$y, lambda = c(1, -1)), "non-negative")
  expect_error(hdhuber(d$x, d$y, pf = rep(1, 4)), "L1 penalty factor")
  expect_error(hdhuber(d$x, d$y, pf2 = rep(1, 4)), "L2 penalty factor")
  expect_error(hdhuber(d$x, d$y, exclude = 9), "out of range")
  expect_error(hdhuber(d$x, d$y, pf = rep(0, 5)), "penalty factors")
  expect_error(hdhuber(matrix(1, 30, 2), d$y), "zero variance")
  expect_error(hdhuber(d$x, factor(d$y > 0)), "numeric")
})

test_that("non-convergence is reported as a warning with partial output", {
  d <- sim_reg(n = 40, p = 10)
  expect_warning(fit <- hdhuber(d$x, d$y, nlambda = 20, maxit = 40),
                 "Convergence not reached")
  expect_true(length(fit$lambda) < 20)
  expect_equal(fit$jerr, -1)
})
