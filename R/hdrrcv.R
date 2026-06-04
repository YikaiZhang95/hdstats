#' Cross-validation for Penalized Wilcoxon Rank Regression
#'
#' Performs k-fold cross-validation for \code{\link{hdrr}}. The held-out
#' prediction error is measured with the median (\eqn{\tau = 0.5}) check loss,
#' i.e. mean absolute error, consistent with the median-regression formulation
#' underlying canonical Wilcoxon rank regression.
#'
#' @param x A numerical matrix with \eqn{n} rows (observations) and \eqn{p}
#'   columns (variables).
#' @param y Response variable of length \eqn{n}.
#' @param lambda Optional; a user-supplied sequence of \code{lambda} values. If
#'   \code{NULL}, \code{\link{hdrr}} selects its own sequence.
#' @param nfolds Number of folds for cross-validation. Defaults to 5.
#' @param foldid Optional vector identifying the fold of each observation. If
#'   provided, it overrides \code{nfolds}.
#' @param ... Additional arguments passed to \code{\link{hdrr}} (e.g.
#'   \code{lam2}, \code{pf}, \code{standardize}).
#'
#' @return An object with S3 class \code{c("cv.hdrr", "cv.hdqr")}. Because it
#'   stores the full-data fit and carries the standard cross-validation fields,
#'   the \code{\link{coef.cv.hdqr}} and \code{\link{predict.cv.hdqr}} methods
#'   apply directly. Components:
#'   \item{lambda}{Candidate \code{lambda} values.}
#'   \item{cvm}{Mean cross-validation error.}
#'   \item{cvsd}{Standard error of the mean cross-validation error.}
#'   \item{cvupper, cvlower}{\code{cvm} plus/minus \code{cvsd}.}
#'   \item{lambda.min}{\code{lambda} achieving the minimum CV error.}
#'   \item{lambda.1se}{Largest \code{lambda} within one standard error of the minimum.}
#'   \item{nzero}{Number of non-zero coefficients at each \code{lambda}.}
#'   \item{hdqr.fit}{The full-data \code{\link{hdrr}} fit (an object of class
#'     \code{c("hdrr", "hdqr")}).}
#'
#' @seealso \code{\link{hdrr}}, \code{\link{cv.hdqr}}
#' @keywords models regression
#' @export
#' @examples
#' set.seed(1)
#' n <- 50; p <- 100
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] * 2 - x[, 2] * 1.5 + rnorm(n)
#' cv.fit <- cv.hdrr(x, y, lam2 = 0.01)
#' coef(cv.fit, s = "lambda.min")
cv.hdrr <- function(x, y, lambda = NULL, nfolds = 5L, foldid, ...) {
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row)
    stop("x and y have different number of observations.")

  hdrr.object <- hdrr(x, y, lambda = lambda, ...)
  lambda <- hdrr.object$lambda
  nz <- sapply(coef(hdrr.object, type = "nonzero"), length)

  if (missing(foldid))
    foldid <- sample(rep(seq(nfolds), length = x.row)) else nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=5 recommended.")

  outlist <- as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    which <- foldid == i
    outlist[[i]] <- hdrr(x = x[!which, , drop = FALSE],
                         y = y[!which], lambda = lambda, ...)
  }

  # rank regression == median QR, so score held-out error with tau = 0.5
  cvstuff <- cvpath.hdqr(outlist, x, y, 0.5, lambda, foldid, x.row, ...)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd

  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd,
              cvupper = cvm + cvsd, cvlower = cvm - cvsd,
              nzero = nz, name = "Wilcoxon rank (MAE)",
              hdqr.fit = hdrr.object)
  obj <- c(out, as.list(getmin(lambda, cvm, cvsd)))
  class(obj) <- c("cv.hdrr", "cv.hdqr")
  obj
}
