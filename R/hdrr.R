#' Canonical (Wilcoxon) Rank Regression via Quantile Regression
#'
#' Fits high-dimensional penalized canonical (Wilcoxon) rank regression by
#' exploiting its equivalence with median quantile regression on the pairwise
#' differences of the data.
#'
#' Jaeckel's Wilcoxon dispersion (equivalently, the Gehan rank statistic) is,
#' up to a positive constant,
#' \deqn{D(\beta) = \sum_{i<j} |(y_i - y_j) - (x_i - x_j)^\top \beta|.}
#' This is exactly the median (\eqn{\tau = 0.5}) quantile-regression check loss
#' applied to the differenced response \eqn{\tilde y_{ij} = y_i - y_j} and
#' differenced design \eqn{\tilde x_{ij} = x_i - x_j} over all pairs
#' \eqn{i < j}. The intercept cancels in the differences, so \code{hdrr} fits
#' the slopes with \code{\link{hdqr}} at \code{tau = 0.5} on the differenced
#' data and then recovers the model intercept as the median of the residuals
#' \eqn{\mathrm{median}(y - X\hat\beta)} at each \code{lambda} (the standard
#' Hettmansperger--McKean location estimate).
#'
#' Note that the differenced design has \eqn{n(n-1)/2} rows, so memory and run
#' time scale as \eqn{O(n^2 p)}.
#'
#' @param x An \eqn{n \times p} matrix of predictors, one observation per row.
#' @param y Response vector of length \eqn{n}.
#' @param lam2 The L2 (ridge) regularization parameter passed to \code{hdqr}.
#'   Default is \code{0.01}.
#' @param ... Further arguments passed to \code{\link{hdqr}}, e.g.
#'   \code{nlambda}, \code{lambda}, \code{lambda.factor}, \code{pf}, \code{pf2},
#'   \code{standardize}, \code{eps}, or \code{maxit}.
#'
#' @return An object of class \code{c("hdrr", "hdqr")}. Because it carries the
#'   standard \code{b0}/\code{beta}/\code{lambda} fields, the \code{\link{coef.hdqr}}
#'   and \code{\link{predict.hdqr}} methods apply directly. Components:
#'   \item{b0}{Intercept sequence (recovered as \code{median(y - x \%*\% beta)}).}
#'   \item{beta}{A \eqn{p \times} \code{length(lambda)} sparse matrix of slopes.}
#'   \item{lambda}{The \code{lambda} sequence used.}
#'   \item{df}{Number of nonzero coefficients at each \code{lambda}.}
#'   \item{qrfit}{The underlying \code{hdqr} fit on the differenced data.}
#'
#' @seealso \code{\link{hdqr}}, \code{\link{coef.hdqr}}, \code{\link{predict.hdqr}}
#'
#' @importFrom stats median
#' @export
#' @examples
#' set.seed(1)
#' n <- 50; p <- 100
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] * 2 - x[, 2] * 1.5 + rnorm(n)
#' fit <- hdrr(x, y, lam2 = 0.01)
#' coef(fit, s = fit$lambda[5])
hdrr <- function(x, y, lam2 = 0.01, ...) {
  x <- as.matrix(x)
  y <- drop(y)
  n <- nrow(x)
  if (length(y) != n)
    stop("x and y have different number of observations")
  if (n < 2L)
    stop("need at least 2 observations for rank regression")

  # all i < j pairwise differences (Gehan/Wilcoxon dispersion)
  pr <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  i1 <- pr[, 1]; i2 <- pr[, 2]
  xd <- x[i1, , drop = FALSE] - x[i2, , drop = FALSE]
  yd <- y[i1] - y[i2]

  # Wilcoxon rank regression == median (tau = 0.5) QR on the differences
  qrfit <- hdqr(xd, yd, tau = 0.5, lam2 = lam2, ...)

  # recover the location intercept: alpha(lambda) = median(y - x beta)
  bmat <- as.matrix(qrfit$beta)
  resid <- as.vector(y) - x %*% bmat
  b0 <- apply(resid, 2, median)
  names(b0) <- colnames(bmat)

  out <- list(b0 = b0, beta = qrfit$beta, lambda = qrfit$lambda,
              df = qrfit$df, npasses = qrfit$npasses, jerr = qrfit$jerr,
              qrfit = qrfit, call = match.call())
  class(out) <- c("hdrr", "hdqr")
  out
}
