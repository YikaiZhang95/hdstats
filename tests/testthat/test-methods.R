test_that("coef and predict are consistent and interpolate lambda", {
  d <- sim_reg(n = 60, p = 15)
  fits <- list(hdhuber(d$x, d$y, nlambda = 10),
               hdqr(d$x, d$y, nlambda = 10),
               hdsvm(d$x, sign(d$y), nlambda = 10))
  for (fit in fits) {
    cf <- coef(fit)
    expect_equal(dim(cf), c(16, 10))
    expect_equal(rownames(cf)[1], "(Intercept)")
    ## predictions are cbind(1, x) %*% coef
    pr <- if (inherits(fit, "hdsvm")) predict(fit, d$x, type = "loss") else predict(fit, d$x)
    expect_equal(pr, as.matrix(cbind(1, d$x) %*% cf), ignore_attr = TRUE)
    ## exact lambda values reproduce the path
    cf5 <- coef(fit, s = fit$lambda[5])
    expect_equal(as.vector(cf5), as.vector(cf[, 5]))
    ## halfway between two lambdas gives the average of the coefficients
    s <- mean(fit$lambda[5:6])
    cfm <- coef(fit, s = s)
    expect_equal(as.vector(cfm), as.vector((cf[, 5] + cf[, 6]) / 2), tolerance = 1e-8)
    ## values outside the range are clamped
    expect_equal(as.vector(coef(fit, s = 1e6)), as.vector(cf[, 1]))
    expect_equal(as.vector(coef(fit, s = 0)), as.vector(cf[, 10]))
    ## nonzero indices
    nz <- coef(fit, type = "nonzero")
    expect_length(nz, 10)
    expect_equal(unname(lengths(nz)), unname(fit$df))
    expect_output(print(fit), "Call:")
  }
})

test_that("cross-validation objects have the documented structure", {
  d <- sim_reg(n = 60, p = 15)
  set.seed(2)
  cvs <- list(cv.hdhuber(d$x, d$y, nlambda = 10, delta = 0.5),
              cv.hdqr(d$x, d$y, nlambda = 10, tau = 0.6),
              cv.hdsvm(d$x, sign(d$y), nlambda = 10))
  fitname <- c("hdhuber.fit", "hdqr.fit", "hdsvm.fit")
  for (i in seq_along(cvs)) {
    cv <- cvs[[i]]
    expect_true(all(c("lambda", "cvm", "cvsd", "cvupper", "cvlower", "nzero",
                      "name", "lambda.min", "lambda.1se", "cvm.min",
                      "cvm.1se", "call", fitname[i]) %in% names(cv)))
    expect_length(cv$cvm, 10)
    expect_equal(cv$cvupper, cv$cvm + cv$cvsd)
    expect_true(all(cv$cvsd >= 0))
    expect_true(cv$lambda.min %in% cv$lambda)
    expect_true(cv$lambda.1se %in% cv$lambda)
    expect_gte(cv$lambda.1se, cv$lambda.min)
    expect_equal(cv$cvm.min, min(cv$cvm))
    expect_equal(cv$cvm.1se, cv$cvm[match(cv$lambda.1se, cv$lambda)])
    ## the CV error at lambda.min beats the null model's error
    expect_lt(cv$cvm.min, cv$cvm[1])
    ## methods
    expect_equal(as.vector(coef(cv, s = "lambda.min")),
                 as.vector(coef(cv[[fitname[i]]], s = cv$lambda.min)))
    expect_equal(as.vector(coef(cv)), as.vector(coef(cv, s = "lambda.1se")))
    expect_equal(dim(coef(cv, s = c(0.1, 0.2))), c(16, 2))
    expect_length(predict(cv, d$x[1:5, ], s = "lambda.min"), 5)
    expect_error(coef(cv, s = TRUE), "Invalid form")
    expect_error(predict(cv, d$x, s = list()), "Invalid form")
    expect_output(print(cv), "Measure:")
    pdf(NULL)
    expect_silent(plot(cv))
    expect_silent(plot(cv, sign.lambda = -1))
    dev.off()
  }
  ## a user-supplied foldid is respected
  foldid <- rep(1:4, length.out = 60)
  cv <- cv.hdhuber(d$x, d$y, nlambda = 10, foldid = foldid)
  expect_length(cv$cvm, 10)
  ## parallel folds give the same result as sequential folds
  cv2 <- cv.hdhuber(d$x, d$y, nlambda = 10, foldid = foldid, ncores = 2)
  expect_equal(cv2$cvm, cv$cvm)
  expect_equal(cv2$cvsd, cv$cvsd)
  cvq <- cv.hdqr(d$x, d$y, nlambda = 10, foldid = foldid, lam2 = 0.1)
  cvq2 <- cv.hdqr(d$x, d$y, nlambda = 10, foldid = foldid, lam2 = 0.1, ncores = 2)
  expect_equal(cvq2$cvm, cvq$cvm)
  expect_error(cv.hdhuber(d$x, d$y, nlambda = 10, ncores = 0), "positive integer")
  expect_error(cv.hdhuber(d$x, d$y, nfolds = 2), "at least 3")
  ## cv.hdsvm with a factor response
  cvf <- cv.hdsvm(d$x, factor(sign(d$y)), nlambda = 10)
  expect_true(all(predict(cvf, d$x, s = "lambda.min") %in% c(-1, 1)))
  expect_true(all(predict(cvf, d$x, s = "lambda.min", type = "loss") != 0))
})

test_that("nonconvex cross-validation objects work end to end", {
  d <- sim_reg(n = 60, p = 15)
  lambda <- 10^seq(0, -1.5, length.out = 6)
  set.seed(3)
  cvq <- cv.nc.hdqr(d$x, d$y, lambda = lambda, lam2 = 0.01, pen = "scad")
  expect_s3_class(cvq, "cv.nc.hdqr")
  expect_equal(cvq$lambda, lambda)
  expect_true(cvq$lambda.min %in% lambda)
  expect_equal(dim(coef(cvq, s = "lambda.min")), c(16, 1))
  expect_length(predict(cvq, d$x[1:3, ], s = c(0.1, 0.2)), 6)
  expect_output(print(cvq), "Measure:")
  expect_output(print(cvq$nchdqr.fit), "Call:")
  pdf(NULL); expect_silent(plot(cvq)); dev.off()

  cvs <- cv.nc.hdsvm(d$x, sign(d$y), lambda = lambda, lam2 = 0.01, pen = "mcp")
  expect_s3_class(cvs, "cv.nc.hdsvm")
  expect_true(all(predict(cvs, d$x, s = "lambda.1se") %in% c(-1, 1)))
  expect_equal(dim(coef(cvs, s = "lambda.min")), c(16, 1))
  expect_output(print(cvs), "Measure:")
  expect_output(print(cvs$nchdsvm.fit), "Call:")
  pdf(NULL); expect_silent(plot(cvs)); dev.off()
})
