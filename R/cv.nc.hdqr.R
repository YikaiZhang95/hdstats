#' Cross-validation for Selecting the Tuning Parameter of Nonconvex Penalized Quantile Regression
#'
#' Performs k-fold cross-validation for \code{\link{nc.hdqr}}. The held-out
#' prediction error is measured with the check loss at the quantile level
#' \code{tau}.
#'
#' @param x A numerical matrix with \eqn{n} rows (observations) and \eqn{p}
#'   columns (variables).
#' @param y Numeric response vector of length \eqn{n}.
#' @param lambda Optional user-supplied sequence of \code{lambda} values
#'   for the nonconvex penalty; see \code{\link{nc.hdqr}}.
#' @param tau The quantile level used in the loss function. Default is 0.5.
#' @param nfolds Number of folds for cross-validation. Default is 5.
#'   The smallest allowed value is 3.
#' @param foldid Optional vector of values between 1 and \code{nfolds}
#'   identifying the fold of each observation. If provided, it overrides
#'   \code{nfolds}.
#' @param ... Additional arguments passed to \code{\link{nc.hdqr}}
#'   (e.g. \code{pen}, \code{aval}, \code{lam2}, \code{lla_step}).
#'
#' @details
#' The function first fits \code{\link{nc.hdqr}} on the full data, then
#' refits the model \code{nfolds} times, each time leaving out one fold.
#' The cross-validation error curve is the average held-out check loss over
#' all observations, and its standard error is computed from the held-out
#' losses.
#'
#' @return
#' An object with S3 class \code{"cv.nc.hdqr"} consisting of
#'   \item{lambda}{the \code{lambda} values used in the fits.}
#'   \item{cvm}{mean cross-validation error, a vector of length
#'     \code{length(lambda)}.}
#'   \item{cvsd}{estimated standard error of \code{cvm}.}
#'   \item{cvupper}{upper curve: \code{cvm + cvsd}.}
#'   \item{cvlower}{lower curve: \code{cvm - cvsd}.}
#'   \item{nzero}{number of nonzero coefficients at each \code{lambda}.}
#'   \item{name}{a text string describing the error measure (for plotting).}
#'   \item{call}{the call that produced this object.}
#'   \item{nchdqr.fit}{a fitted \code{\link{nc.hdqr}} object for the full
#'     data.}
#'   \item{lambda.min}{the \code{lambda} achieving the minimum
#'     cross-validation error.}
#'   \item{lambda.1se}{the largest \code{lambda} whose cross-validation error
#'     is within one standard error of the minimum.}
#'   \item{cvm.min}{cross-validation error at \code{lambda.min}.}
#'   \item{cvm.1se}{cross-validation error at \code{lambda.1se}.}
#' @seealso \code{\link{nc.hdqr}}, \code{\link{coef.cv.nc.hdqr}},
#'   \code{\link{predict.cv.nc.hdqr}}, \code{\link{plot.cv.nc.hdqr}}
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
#' lambda <- 10^(seq(0, -2, length.out = 10))
#' \donttest{
#' cv.nc.fit <- cv.nc.hdqr(x = x, y = y, tau = 0.5, lambda = lambda,
#'                         lam2 = 0.01, pen = "scad")
#' cv.nc.fit
#' }
cv.nc.hdqr <- function(x, y, lambda = NULL, tau = 0.5, nfolds = 5L, foldid, ...) {
  ####################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  ####################################################################
  nc.hdqr.object <- nc.hdqr(x, y, lambda = lambda, tau = tau, ...)
  lambda <- nc.hdqr.object$nc.lambda
  nz <- sapply(coef(nc.hdqr.object, type = "nonzero"), length)
  if (missing(foldid))
    foldid <- sample(rep(seq(nfolds), length = x.row)) else nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be at least 3; nfolds = 5 recommended.")
  outlist <- as.list(seq(nfolds))
  ## fit the model nfold times and save them
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- nc.hdqr(x = x[!which, , drop = FALSE],
                            y = y[!which], tau = tau, lambda = lambda, ...)
  }
  ## select the lambda according to predmat
  cvstuff <- cvpath.hdqr(outlist, x, y, tau, lambda, foldid, x.row, ...)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd,
              cvupper = cvm + cvsd, cvlower = cvm - cvsd, nzero = nz,
              name = cvname, nchdqr.fit = nc.hdqr.object)
  obj <- c(out, as.list(getmin(lambda, cvm, cvsd)))
  obj$call <- match.call()
  class(obj) <- "cv.nc.hdqr"
  obj
}
