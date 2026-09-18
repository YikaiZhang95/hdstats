#' Cross-validation for Selecting the Tuning Parameter in Penalized Quantile Regression
#'
#' Performs k-fold cross-validation for \code{\link{hdqr}}. The held-out
#' prediction error is measured with the check loss at the same quantile
#' level \code{tau} used for fitting.
#'
#' @param x A numerical matrix with \eqn{n} rows (observations) and \eqn{p}
#'   columns (variables).
#' @param y Numeric response vector of length \eqn{n}.
#' @param lambda Optional; a user-supplied sequence of \code{lambda} values.
#'   If \code{NULL}, \code{\link{hdqr}} selects its own sequence.
#' @param tau The quantile level used in the loss function. Default is 0.5.
#' @param nfolds Number of folds for cross-validation. Default is 5.
#'   The smallest allowed value is 3.
#' @param foldid Optional vector of values between 1 and \code{nfolds}
#'   identifying the fold of each observation. If provided, it overrides
#'   \code{nfolds}.
#' @param ncores Number of processes used to fit the folds in parallel
#'   (with the \pkg{parallel} package: forked processes on Unix-alikes, a
#'   socket cluster on Windows). Default is 1, i.e. the folds are fitted
#'   sequentially.
#' @param ... Additional arguments passed to \code{\link{hdqr}}.
#'
#' @details
#' The function first fits \code{\link{hdqr}} on the full data to obtain
#' the \code{lambda} sequence, then refits the model \code{nfolds} times,
#' each time leaving out one fold. The cross-validation error curve is the
#' average held-out check loss over all observations, and its standard error
#' is computed from the held-out losses.
#'
#' @return
#' An object with S3 class \code{"cv.hdqr"} consisting of
#'   \item{lambda}{the \code{lambda} values used in the fits.}
#'   \item{cvm}{mean cross-validation error, a vector of length
#'     \code{length(lambda)}.}
#'   \item{cvsd}{estimated standard error of \code{cvm}.}
#'   \item{cvupper}{upper curve: \code{cvm + cvsd}.}
#'   \item{cvlower}{lower curve: \code{cvm - cvsd}.}
#'   \item{nzero}{number of nonzero coefficients at each \code{lambda}.}
#'   \item{name}{a text string describing the error measure (for plotting).}
#'   \item{call}{the call that produced this object.}
#'   \item{hdqr.fit}{a fitted \code{\link{hdqr}} object for the full data.}
#'   \item{lambda.min}{the \code{lambda} achieving the minimum
#'     cross-validation error.}
#'   \item{lambda.1se}{the largest \code{lambda} whose cross-validation error
#'     is within one standard error of the minimum.}
#'   \item{cvm.min}{cross-validation error at \code{lambda.min}.}
#'   \item{cvm.1se}{cross-validation error at \code{lambda.1se}.}
#' @seealso \code{\link{hdqr}}, \code{\link{coef.cv.hdqr}},
#'   \code{\link{predict.cv.hdqr}}, \code{\link{plot.cv.hdqr}}
#' @keywords models regression
#' @export
#' @examples
#' set.seed(315)
#' n <- 100
#' p <- 400
#' x <- matrix(data = rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
#' beta_star <- c(c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.3), rep(0, (p - 7)))
#' eps <- rnorm(n, mean = 0, sd = 1)
#' y <- x %*% beta_star + eps
#' cv.fit <- cv.hdqr(x = x, y = y, tau = 0.5)
#' cv.fit
#' plot(cv.fit)
cv.hdqr <- function(x, y, lambda = NULL, tau = 0.5, nfolds = 5L, foldid,
                    ncores = 1L, ...) {
  ####################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  ####################################################################
  hdqr.object <- hdqr(x, y, lambda = lambda, tau = tau, ...)
  lambda <- hdqr.object$lambda
  nz <- sapply(coef(hdqr.object, type = "nonzero"), length)
  if (missing(foldid))
    foldid <- sample(rep(seq(nfolds), length = x.row)) else nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be at least 3; nfolds = 5 recommended.")
  ## fit the model nfold times and save them
  dots <- list(...)
  outlist <- cv_folds(nfolds, function(i) {
    which <- foldid == i
    do.call(hdqr, c(list(x = x[!which, , drop = FALSE], y = y[!which],
                         tau = tau, lambda = lambda), dots))
  }, ncores)
  ## select the lambda according to predmat
  cvstuff <- cvpath.hdqr(outlist, x, y, tau, lambda, foldid, x.row, ...)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd,
              cvupper = cvm + cvsd, cvlower = cvm - cvsd,
              nzero = nz, name = cvname,
              hdqr.fit = hdqr.object)
  obj <- c(out, as.list(getmin(lambda, cvm, cvsd)))
  obj$call <- match.call()
  class(obj) <- "cv.hdqr"
  obj
}

cvpath.hdqr <- function(outlist, x, y, tau, lambda, foldid, x.row, ...) {
  nfolds <- max(foldid)
  predmat <- matrix(NA, x.row, length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    whichfold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[whichfold, , drop = FALSE])
    nlami <- length(fitobj$lambda)
    predmat[whichfold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- check_loss(y - predmat, tau)
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- colMeans(cvraw, na.rm = TRUE)
  scaled <- scale(cvraw, cvm, FALSE)^2
  cvsd <- sqrt(colMeans(scaled, na.rm = TRUE) / (N - 1))
  list(cvm = cvm, cvsd = cvsd, cvraw = cvraw,
       name = paste0("Check loss (tau = ", format(tau), ")"))
}
