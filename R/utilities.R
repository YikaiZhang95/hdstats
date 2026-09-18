######################################################################
## Internal helper functions.
##
## Several of these functions are minor modifications of, or directly
## copied from, the glmnet package:
##   Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
##   Regularization Paths for Generalized Linear Models via Coordinate
##   Descent. Journal of Statistical Software, 33(1), 1-22.
##   doi:10.18637/jss.v033.i01
######################################################################

error.bars <- function(x, upper, lower, width = 0.02, ...) {
  ## adapted from glmnet
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

getmin <- function(lambda, cvm, cvsd) {
  ## adapted from glmnet
  cvmin <- min(cvm, na.rm = TRUE)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin], na.rm = TRUE)
  cvmin2 <- min(cvm[!is.na(cvsd)])
  lambda.min2 <- max(lambda[cvm[!is.na(cvsd)] <= cvmin2], na.rm = TRUE)
  idmin <- match(lambda.min2, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm[!is.na(cvsd)] <= semin
  lambda.1se <- max(lambda[idmin])
  id1se <- match(lambda.1se, lambda)
  cv.1se <- cvm[id1se]
  list(lambda.min = lambda.min, lambda.1se = lambda.1se,
       cvm.min = cvmin, cvm.1se = cv.1se)
}

lambda.interp <- function(lambda, s) {
  ## adapted from glmnet
  ## lambda is the (decreasing) sequence produced by the model;
  ## s is the new vector at which evaluations are required.
  ## Returns left and right indices and the fraction such that the
  ## interpolated value is frac * left + (1 - frac) * right.
  if (length(lambda) == 1) {
    nums <- length(s)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    s[s > max(lambda)] <- max(lambda)
    s[s < min(lambda)] <- min(lambda)
    k <- length(lambda)
    sfrac <- (lambda[1] - s) / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}

## Translate the integer error flag returned by the C++ kernels into a
## message. n = 1 signals an error, n = -1 a warning, n = 0 nothing.
err_pmax <- function(n, nalam, maxit, pmax) {
  msg <- ""
  if (n > 0) {
    if (n < 7777)
      msg <- "Memory allocation error"
    if (n == 7777)
      msg <- "All used predictors have zero variance"
    if (n == 10000)
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in hdstats C++ code:", msg)
  }
  if (n < 0) {
    if (n > -10000)
      msg <- paste0("Convergence not reached after maxit = ", maxit,
                    " coordinate-descent passes; solutions for the ", nalam,
                    " largest lambda values were returned")
    if (n < -10000)
      msg <- paste0("Number of nonzero coefficients along the path exceeds pmax = ",
                    pmax, " at the ", -n - 10000,
                    "th lambda value; solutions for larger lambdas were returned")
    n <- -1
    msg <- paste("from hdstats C++ code:", msg)
  }
  list(n = n, msg = msg)
}

## Convert the raw list returned by the C++ kernels into the coefficient
## representation used by the fitted objects (sparse beta, b0, df, lambda).
getoutput <- function(fit, maxit, pmax, nvars, vnames) {
  nalam <- fit$nalam
  errmsg <- err_pmax(fit$jerr, nalam, maxit, pmax)
  if (errmsg$n == 1)
    stop(errmsg$msg, call. = FALSE)
  if (nalam < 1L)
    stop("no lambda value was completed; increase maxit or eps", call. = FALSE)
  if (errmsg$n == -1)
    warning(errmsg$msg, call. = FALSE)
  nbeta <- fit$nbeta[seq(nalam)]
  nbetamax <- max(nbeta)
  lam <- fit$alam[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  dd <- c(nvars, nalam)
  if (nbetamax > 0) {
    beta <- matrix(fit$beta[seq(pmax * nalam)],
                   pmax, nalam)[seq(nbetamax), , drop = FALSE]
    df <- apply(abs(beta) > 0, 2, sum)
    ja <- fit$ibeta[seq(nbetamax)]
    oja <- order(ja)
    ja <- rep(ja[oja], nalam)
    ibeta <- cumsum(c(1, rep(nbetamax, nalam)))
    beta <- new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames),
                x = as.vector(beta[oja, ]),
                p = as.integer(ibeta - 1), i = as.integer(ja - 1))
  } else {
    beta <- zeromat(nvars, nalam, vnames, stepnames)
    df <- rep(0L, nalam)
  }
  b0 <- fit$b0
  if (!is.null(b0)) {
    b0 <- b0[seq(nalam)]
    names(b0) <- stepnames
  }
  list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
}

## The first lambda of an automatically generated sequence is +Inf in
## effect (all coefficients zero); replace it by a finite extrapolation of
## the geometric sequence so that log(lambda) can be plotted.
lamfix <- function(lam) {
  llam <- log(lam)
  lam[1] <- exp(2 * llam[2] - llam[3])
  lam
}

nonzero <- function(beta, bystep = FALSE) {
  ## adapted from glmnet; beta should be in 'dgCMatrix' format
  ns <- ncol(beta)
  if (nrow(beta) == 1) {
    if (bystep) {
      apply(beta, 2, function(x) if (abs(x) > 0) 1 else NULL)
    } else {
      if (any(abs(beta) > 0)) 1 else NULL
    }
  } else {
    beta <- t(beta)
    which <- diff(beta@p)
    which <- seq(which)[which > 0]
    if (bystep) {
      nzel <- function(x, which) if (any(x)) which[x] else NULL
      beta <- abs(as.matrix(beta[, which])) > 0
      if (ns == 1) {
        apply(beta, 2, nzel, which)
      } else apply(beta, 1, nzel, which)
    } else which
  }
}

zeromat <- function(nvars, nalam, vnames, stepnames) {
  ca <- rep(0, nalam)
  ia <- seq(nalam + 1)
  ja <- rep(1, nalam)
  dd <- c(nvars, nalam)
  new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames),
      x = as.vector(ca), p = as.integer(ia - 1), i = as.integer(ja - 1))
}

## Loss functions used for cross-validation -----------------------------

check_loss <- function(r, tau) {
  0.5 * (abs(r) + (2.0 * tau - 1.0) * r)
}

huber_loss <- function(u, delta) {
  ## Huber loss (quadratic near 0, linear beyond delta)
  ifelse(abs(u) <= delta,
         0.5 * u^2,
         delta * (abs(u) - 0.5 * delta))
}

## Derivatives of the nonconvex penalties (for the LLA weights) ----------

deriv_scad <- function(u, lambda, a = 3.7) {
  u <- abs(u) # u must be nonnegative
  lambda * (u <= lambda) + (a * lambda - u) / (a - 1) *
    (u > lambda) * (u <= a * lambda)
}

deriv_mcp <- function(u, lambda, a = 2) {
  u <- abs(u) # u must be nonnegative
  (lambda - u / a) * (u <= a * lambda)
}
