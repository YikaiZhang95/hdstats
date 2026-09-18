#' Cross-validation for Selecting the Tuning Parameter of Nonconvex Penalized SVM
#'
#' Performs k-fold cross-validation for \code{\link{nc.hdsvm}}. The held-out
#' prediction error is measured with the hinge loss.
#'
#' @param x A numerical matrix with \eqn{n} rows (observations) and \eqn{p}
#'   columns (variables).
#' @param y Binary response of length \eqn{n}; see \code{\link{hdsvm}}.
#' @param lambda Optional user-supplied sequence of \code{lambda} values
#'   for the nonconvex penalty; see \code{\link{nc.hdsvm}}.
#' @param nfolds Number of folds for cross-validation. Default is 5.
#'   The smallest allowed value is 3.
#' @param foldid Optional vector of values between 1 and \code{nfolds}
#'   identifying the fold of each observation. If provided, it overrides
#'   \code{nfolds}.
#' @param ncores Number of processes used to fit the folds in parallel
#'   (with the \pkg{parallel} package: forked processes on Unix-alikes, a
#'   socket cluster on Windows). Default is 1, i.e. the folds are fitted
#'   sequentially.
#' @param ... Additional arguments passed to \code{\link{nc.hdsvm}}
#'   (e.g. \code{pen}, \code{aval}, \code{lam2}, \code{lla_step}).
#'
#' @details
#' The function first fits \code{\link{nc.hdsvm}} on the full data, then
#' refits the model \code{nfolds} times, each time leaving out one fold.
#' The cross-validation error curve is the average held-out hinge loss over
#' all observations, and its standard error is computed from the held-out
#' losses.
#'
#' @return
#' An object with S3 class \code{"cv.nc.hdsvm"} consisting of
#'   \item{lambda}{the \code{lambda} values used in the fits.}
#'   \item{cvm}{mean cross-validation error, a vector of length
#'     \code{length(lambda)}.}
#'   \item{cvsd}{estimated standard error of \code{cvm}.}
#'   \item{cvupper}{upper curve: \code{cvm + cvsd}.}
#'   \item{cvlower}{lower curve: \code{cvm - cvsd}.}
#'   \item{nzero}{number of nonzero coefficients at each \code{lambda}.}
#'   \item{name}{a text string describing the error measure (for plotting).}
#'   \item{call}{the call that produced this object.}
#'   \item{nchdsvm.fit}{a fitted \code{\link{nc.hdsvm}} object for the full
#'     data.}
#'   \item{lambda.min}{the \code{lambda} achieving the minimum
#'     cross-validation error.}
#'   \item{lambda.1se}{the largest \code{lambda} whose cross-validation error
#'     is within one standard error of the minimum.}
#'   \item{cvm.min}{cross-validation error at \code{lambda.min}.}
#'   \item{cvm.1se}{cross-validation error at \code{lambda.1se}.}
#' @seealso \code{\link{nc.hdsvm}}, \code{\link{coef.cv.nc.hdsvm}},
#'   \code{\link{predict.cv.nc.hdsvm}}, \code{\link{plot.cv.nc.hdsvm}}
#' @keywords models classification
#' @export
#' @examples
#' set.seed(315)
#' n <- 100
#' p <- 400
#' x1 <- matrix(rnorm(n / 2 * p, -0.25, 0.1), n / 2)
#' x2 <- matrix(rnorm(n / 2 * p, 0.25, 0.1), n / 2)
#' x <- rbind(x1, x2)
#' beta <- 0.1 * rnorm(p)
#' prob <- plogis(c(x %*% beta))
#' y <- 2 * rbinom(n, 1, prob) - 1
#' lambda <- 10^(seq(0, -2, length.out = 10))
#' \donttest{
#' cv.nc.fit <- cv.nc.hdsvm(x = x, y = y, lambda = lambda, lam2 = 0.01,
#'                          pen = "scad")
#' cv.nc.fit
#' }
cv.nc.hdsvm <- function(x, y, lambda = NULL, nfolds = 5L, foldid, ncores = 1L, ...) {
  ####################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  ####################################################################
  nc.hdsvm.object <- nc.hdsvm(x, y, lambda = lambda, ...)
  lambda <- nc.hdsvm.object$nc.lambda
  nz <- sapply(coef(nc.hdsvm.object, type = "nonzero"), length)
  if (missing(foldid))
    foldid <- sample(rep(seq(nfolds), length = x.row)) else nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be at least 3; nfolds = 5 recommended.")
  ## fit the model nfold times and save them
  dots <- list(...)
  outlist <- cv_folds(nfolds, function(i) {
    which <- foldid == i
    do.call(nc.hdsvm, c(list(x = x[!which, , drop = FALSE], y = y[!which],
                             lambda = lambda), dots))
  }, ncores)
  ## select the lambda according to predmat
  cvstuff <- cvpath.hdsvm(outlist, x, y, lambda, foldid, x.row, ...)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd,
              cvupper = cvm + cvsd, cvlower = cvm - cvsd, nzero = nz,
              name = cvname, nchdsvm.fit = nc.hdsvm.object)
  obj <- c(out, as.list(getmin(lambda, cvm, cvsd)))
  obj$call <- match.call()
  class(obj) <- "cv.nc.hdsvm"
  obj
}
