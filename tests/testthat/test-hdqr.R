test_that("hdqr returns a well-formed solution path", {
  d <- sim_reg()
  fit <- hdqr(d$x, d$y, nlambda = 20)
  expect_s3_class(fit, "hdqr")
  expect_equal(dim(fit$beta), c(ncol(d$x), 20))
  expect_length(fit$b0, 20)
  expect_true(all(diff(fit$lambda) < 0))
  expect_equal(unname(fit$df), unname(colSums(as.matrix(fit$beta) != 0)))
  expect_equal(fit$tau, 0.5)
  expect_equal(fit$jerr, 0)
  expect_equal(fit$df[1], 0)
  ## the default tau is 0.5
  fit2 <- hdqr(d$x, d$y, tau = 0.5, nlambda = 20)
  expect_equal(as.matrix(fit$beta), as.matrix(fit2$beta))
  ## the null model intercept minimizes the check loss (a tau-quantile of y)
  for (tau in c(0.5, 0.8)) {
    f <- hdqr(d$x, d$y, tau = tau, nlambda = 5)
    expect_equal(mean(check_loss_ref(d$y - f$b0[1], tau)),
                 mean(check_loss_ref(d$y - quantile(d$y, tau, type = 1), tau)),
                 tolerance = 1e-6)
  }
})

test_that("hdqr agrees with quantreg on unpenalized and lasso problems", {
  skip_if_not_installed("quantreg")
  d <- sim_reg(n = 200, p = 5, seed = 3)
  for (tau in c(0.3, 0.5, 0.75)) {
    fit <- hdqr(d$x, d$y, tau = tau, lambda = c(1, 0), lam2 = 0,
                standardize = FALSE, eps = 1e-12)
    cf <- path_coef(fit, 2)
    rq <- quantreg::rq(d$y ~ d$x, tau = tau)
    obj_fit <- qr_obj(cf$b0, cf$beta, d$x, d$y, tau, 0, 0)
    obj_rq <- qr_obj(coef(rq)[1], coef(rq)[-1], d$x, d$y, tau, 0, 0)
    expect_lt(obj_fit - obj_rq, 1e-3 * obj_rq)
    expect_equal(c(cf$b0, cf$beta), unname(coef(rq)), tolerance = 0.05)
  }
  ## lasso: quantreg penalizes sum(rho) + lambda * sum|beta|, hdstats
  ## penalizes mean(rho) + lambda * sum|beta|
  n <- nrow(d$x)
  lam_rq <- 20
  rq <- quantreg::rq(d$y ~ d$x, tau = 0.5, method = "lasso", lambda = lam_rq)
  fit <- hdqr(d$x, d$y, tau = 0.5, lambda = c(1, lam_rq / n), lam2 = 0,
              standardize = FALSE, eps = 1e-12)
  cf <- path_coef(fit, 2)
  obj_fit <- qr_obj(cf$b0, cf$beta, d$x, d$y, 0.5, lam_rq / n, 0)
  obj_rq <- qr_obj(coef(rq)[1], coef(rq)[-1], d$x, d$y, 0.5, lam_rq / n, 0)
  expect_lt(obj_fit - obj_rq, 1e-3 * obj_rq)
})

test_that("hdqr solutions are (near) minimizers of the penalized check loss", {
  d <- sim_reg(n = 80, p = 20)
  lam2 <- 0.05
  fit <- hdqr(d$x, d$y, tau = 0.4, lam2 = lam2, nlambda = 15,
              standardize = FALSE, eps = 1e-12)
  set.seed(11)
  for (l in c(8, 15)) {
    cf <- path_coef(fit, l)
    lambda <- fit$lambda[l]
    obj0 <- qr_obj(cf$b0, cf$beta, d$x, d$y, 0.4, lambda, lam2)
    ## the zero model is worse
    expect_lt(obj0, qr_obj(median(d$y), rep(0, 20), d$x, d$y, 0.4, lambda, lam2))
    ## no random perturbation of moderate size improves the objective
    for (k in 1:30) {
      db <- rnorm(20, sd = 0.02); db0 <- rnorm(1, sd = 0.02)
      expect_gte(qr_obj(cf$b0 + db0, cf$beta + db, d$x, d$y, 0.4, lambda, lam2),
                 obj0 - 1e-5)
    }
  }
})

test_that("is_exact = TRUE runs and stays close to the smoothed solution", {
  d <- sim_reg(n = 60, p = 10)
  f1 <- hdqr(d$x, d$y, nlambda = 8, is_exact = FALSE)
  f2 <- hdqr(d$x, d$y, nlambda = 8, is_exact = TRUE)
  expect_equal(f1$lambda, f2$lambda)
  expect_equal(as.matrix(f1$beta), as.matrix(f2$beta), tolerance = 0.05)
})

test_that("hdqr validates its inputs", {
  d <- sim_reg(n = 30, p = 5)
  expect_error(hdqr(d$x, d$y, tau = 0), "tau")
  expect_error(hdqr(d$x, d$y, tau = 1), "tau")
  expect_error(hdqr(d$x, d$y, hval = 0), "hval")
  expect_error(hdqr(d$x, d$y[-1]), "different number of observations")
  expect_error(hdqr(d$x, d$y, pf = rep(0, 5)), "penalty factors")
  expect_error(hdqr(d$x, d$y, pf2 = rep(-1, 5)), "penalty factors")
  expect_error(hdqr(matrix(1, 30, 2), d$y), "zero variance")
})

test_that("nc.hdqr honours pen, aval, and lambda handling", {
  d <- sim_reg(n = 80, p = 20)
  lambda <- 10^seq(0, -2, length.out = 8)
  f_scad <- nc.hdqr(d$x, d$y, lambda = lambda, lam2 = 0.01, pen = "scad")
  expect_s3_class(f_scad, "nc.hdqr")
  expect_equal(f_scad$lambda, lambda)
  expect_equal(f_scad$nc.lambda, lambda)
  expect_equal(f_scad$pen, "scad")
  expect_equal(f_scad$aval, 3.7)
  ## case-insensitive penalty names
  f_SCAD <- nc.hdqr(d$x, d$y, lambda = lambda, lam2 = 0.01, pen = "SCAD")
  expect_equal(as.matrix(f_SCAD$beta), as.matrix(f_scad$beta))
  expect_error(nc.hdqr(d$x, d$y, lambda = lambda, pen = "lasso"))
  ## aval changes the fit
  f_a <- nc.hdqr(d$x, d$y, lambda = lambda, lam2 = 0.01, pen = "scad", aval = 20)
  expect_equal(f_a$aval, 20)
  expect_false(isTRUE(all.equal(as.matrix(f_a$beta), as.matrix(f_scad$beta))))
  expect_error(nc.hdqr(d$x, d$y, lambda = lambda, pen = "scad", aval = 1), "aval")
  expect_error(nc.hdqr(d$x, d$y, lambda = lambda, pen = "mcp", aval = 0.5), "aval")
  ## MCP
  f_mcp <- nc.hdqr(d$x, d$y, lambda = lambda, lam2 = 0.01, pen = "mcp")
  expect_equal(f_mcp$aval, 2)
  expect_equal(dim(f_mcp$beta), c(20, 8))
  ## unsorted lambda is sorted, results identical
  f_u <- nc.hdqr(d$x, d$y, lambda = rev(lambda), lam2 = 0.01)
  expect_equal(f_u$lambda, lambda)
  expect_equal(as.matrix(f_u$beta), as.matrix(f_scad$beta))
  ## lambda = NULL uses the lasso sequence
  f_null <- nc.hdqr(d$x, d$y, lam2 = 0.01, nlambda = 10)
  expect_length(f_null$lambda, 10)
  expect_true(all(diff(f_null$lambda) < 0))
  ## ini_beta requires lambda and must have the right length
  expect_error(nc.hdqr(d$x, d$y, ini_beta = rep(0, 20)), "lambda")
  expect_warning(nc.hdqr(d$x, d$y, lambda = lambda, ini_beta = rep(0, 3)),
                 "ini_beta")
  ## a nonconvex penalty shrinks large coefficients less than the lasso
  f_lasso <- hdqr(d$x, d$y, lambda = lambda, lam2 = 0.01)
  l <- 4
  expect_gt(sum(abs(as.matrix(f_scad$beta)[1:3, l])),
            sum(abs(as.matrix(f_lasso$beta)[1:3, l])))
})
