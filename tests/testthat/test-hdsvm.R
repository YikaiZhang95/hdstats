test_that("hdsvm returns a well-formed solution path", {
  d <- sim_class()
  fit <- hdsvm(d$x, d$y, nlambda = 20, lam2 = 0.01)
  expect_s3_class(fit, "hdsvm")
  expect_equal(dim(fit$beta), c(ncol(d$x), 20))
  expect_length(fit$b0, 20)
  expect_true(all(diff(fit$lambda) < 0))
  expect_equal(unname(fit$df), unname(colSums(as.matrix(fit$beta) != 0)))
  expect_equal(fit$jerr, 0)
  expect_equal(fit$df[1], 0)
  expect_equal(fit$hval, 1)
  pr <- predict(fit, d$x, s = fit$lambda[10])
  expect_true(all(pr %in% c(-1, 1)))
  expect_gt(mean(pr == d$y), 0.8)
})

test_that("factor and 0/1 responses give the same fit as -1/1", {
  d <- sim_class()
  f1 <- hdsvm(d$x, d$y, nlambda = 10)
  f2 <- hdsvm(d$x, factor(d$y, labels = c("no", "yes")), nlambda = 10)
  f3 <- hdsvm(d$x, (d$y + 1) / 2, nlambda = 10)
  expect_equal(as.matrix(f1$beta), as.matrix(f2$beta))
  expect_equal(as.matrix(f1$beta), as.matrix(f3$beta))
  expect_error(hdsvm(d$x, sample(1:3, nrow(d$x), TRUE)), "two levels")
  expect_error(hdsvm(d$x, rep(1, nrow(d$x))), "two levels")
})

test_that("hdsvm with is_exact = TRUE matches the exact SVM solution of e1071", {
  skip_if_not_installed("e1071")
  d <- sim_class(n = 100, p = 20)
  n <- nrow(d$x)
  for (lam2 in c(0.05, 0.5)) {
    ## ridge-only problem: (1/n) sum hinge + lam2/2 ||beta||^2 is the standard
    ## linear SVM with cost = 1 / (n * lam2)
    fit <- hdsvm(d$x, d$y, lambda = c(1, 0), lam2 = lam2, standardize = FALSE,
                 eps = 1e-12, is_exact = TRUE)
    cf <- path_coef(fit, 2)
    m <- e1071::svm(d$x, factor(d$y), kernel = "linear", cost = 1 / (n * lam2),
                    scale = FALSE, tolerance = 1e-6)
    w <- drop(t(m$coefs) %*% m$SV)
    b <- -m$rho
    if (cor(drop(d$x %*% w + b), d$y) < 0) { w <- -w; b <- -b }
    obj_fit <- svm_obj(cf$b0, cf$beta, d$x, d$y, 0, lam2)
    obj_ref <- svm_obj(b, w, d$x, d$y, 0, lam2)
    expect_lt(abs(obj_fit - obj_ref), 1e-4 * obj_ref)
    expect_equal(cf$beta, unname(w), tolerance = 0.05)
  }
})

test_that("the smoothed hdsvm solution is close to the exact hinge-loss minimizer", {
  d <- sim_class(n = 100, p = 20)
  lam2 <- 0.05
  f_smooth <- hdsvm(d$x, d$y, lam2 = lam2, nlambda = 15, standardize = FALSE)
  f_exact <- hdsvm(d$x, d$y, lam2 = lam2, nlambda = 15, standardize = FALSE,
                   is_exact = TRUE)
  expect_equal(f_smooth$lambda, f_exact$lambda)
  for (l in c(8, 15)) {
    lambda <- f_smooth$lambda[l]
    cs <- path_coef(f_smooth, l)
    ce <- path_coef(f_exact, l)
    o_s <- svm_obj(cs$b0, cs$beta, d$x, d$y, lambda, lam2)
    o_e <- svm_obj(ce$b0, ce$beta, d$x, d$y, lambda, lam2)
    o_0 <- svm_obj(0, rep(0, 20), d$x, d$y, lambda, lam2)
    ## the exact refinement improves on the smoothed solution, which is
    ## itself within a few percent of the optimum and far from the null model
    expect_lt(o_e, o_s + 1e-6)
    expect_lt(o_s, 1.1 * o_e)
    expect_lt(o_s, 0.75 * o_0)
    big <- abs(ce$beta) > 0.2
    expect_equal(sign(cs$beta[big]), sign(ce$beta[big]))
  }
})

test_that("hdsvm validates its inputs", {
  d <- sim_class(n = 30, p = 5)
  expect_error(hdsvm(d$x, d$y[-1]), "different number of observations")
  expect_error(hdsvm(d$x, d$y, hval = -1), "hval")
  expect_error(hdsvm(d$x, d$y, pf = rep(0, 5)), "penalty factors")
  expect_error(hdsvm(d$x, d$y, lambda = c(0.1, -0.1)), "non-negative")
})

test_that("nc.hdsvm honours pen, aval, and lambda handling", {
  d <- sim_class(n = 100, p = 20)
  lambda <- 10^seq(0, -2, length.out = 6)
  f_scad <- nc.hdsvm(d$x, d$y, lambda = lambda, lam2 = 0.01, pen = "scad")
  expect_s3_class(f_scad, "nc.hdsvm")
  expect_equal(f_scad$lambda, lambda)
  expect_equal(f_scad$aval, 3.7)
  f_MCP <- nc.hdsvm(d$x, d$y, lambda = lambda, lam2 = 0.01, pen = "MCP")
  expect_equal(f_MCP$pen, "mcp")
  expect_equal(f_MCP$aval, 2)
  f_a <- nc.hdsvm(d$x, d$y, lambda = lambda, lam2 = 0.01, pen = "mcp", aval = 5)
  expect_false(isTRUE(all.equal(as.matrix(f_a$beta), as.matrix(f_MCP$beta))))
  expect_error(nc.hdsvm(d$x, d$y, lambda = lambda, pen = "foo"))
  f_null <- nc.hdsvm(d$x, d$y, lam2 = 0.01, nlambda = 8)
  expect_length(f_null$lambda, 8)
  pr <- predict(f_scad, d$x, s = lambda[4])
  expect_true(all(pr %in% c(-1, 1)))
})
