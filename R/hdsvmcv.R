#' Cross-validation for Selecting the Tuning Parameter in the Penalized SVM
#'
#' Performs k-fold cross-validation for \code{\link{hdsvm}}. The held-out
#' prediction error is measured with the hinge loss.
#'
#' @param x A numerical matrix with \eqn{n} rows (observations) and \eqn{p}
#'   columns (variables).
#' @param y Binary response of length \eqn{n}; see \code{\link{hdsvm}}.
#' @param lambda Optional; a user-supplied sequence of \code{lambda} values.
#'   If \code{NULL}, \code{\link{hdsvm}} selects its own sequence.
#' @param nfolds Number of folds for cross-validation. Default is 5.
#'   The smallest allowed value is 3.
#' @param foldid Optional vector of values between 1 and \code{nfolds}
#'   identifying the fold of each observation. If provided, it overrides
#'   \code{nfolds}.
#' @param ncores Number of processes used to fit the folds in parallel
#'   (with the \pkg{parallel} package: forked processes on Unix-alikes, a
#'   socket cluster on Windows). Default is 1, i.e. the folds are fitted
#'   sequentially.
#' @param ... Additional arguments passed to \code{\link{hdsvm}}.
#'
#' @details
#' The function first fits \code{\link{hdsvm}} on the full data to obtain
#' the \code{lambda} sequence, then refits the model \code{nfolds} times,
#' each time leaving out one fold. The cross-validation error curve is the
#' average held-out hinge loss over all observations, and its standard error
#' is computed from the held-out losses.
#'
#' @return
#' An object with S3 class \code{"cv.hdsvm"} consisting of
#'   \item{lambda}{the \code{lambda} values used in the fits.}
#'   \item{cvm}{mean cross-validation error, a vector of length
#'     \code{length(lambda)}.}
#'   \item{cvsd}{estimated standard error of \code{cvm}.}
#'   \item{cvupper}{upper curve: \code{cvm + cvsd}.}
#'   \item{cvlower}{lower curve: \code{cvm - cvsd}.}
#'   \item{nzero}{number of nonzero coefficients at each \code{lambda}.}
#'   \item{name}{a text string describing the error measure (for plotting).}
#'   \item{call}{the call that produced this object.}
#'   \item{hdsvm.fit}{a fitted \code{\link{hdsvm}} object for the full data.}
#'   \item{lambda.min}{the \code{lambda} achieving the minimum
#'     cross-validation error.}
#'   \item{lambda.1se}{the largest \code{lambda} whose cross-validation error
#'     is within one standard error of the minimum.}
#'   \item{cvm.min}{cross-validation error at \code{lambda.min}.}
#'   \item{cvm.1se}{cross-validation error at \code{lambda.1se}.}
#' @seealso \code{\link{hdsvm}}, \code{\link{coef.cv.hdsvm}},
#'   \code{\link{predict.cv.hdsvm}}, \code{\link{plot.cv.hdsvm}}
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
#' cv.fit <- cv.hdsvm(x, y, lam2 = 0.01)
#' cv.fit
#' plot(cv.fit)
cv.hdsvm <- function(x, y, lambda = NULL, nfolds = 5L, foldid, ncores = 1L, ...) {
  ####################################################################
  ## data setup
  y <- drop(y)
  x <- as.matrix(x)
  x.row <- as.integer(NROW(x))
  if (length(y) != x.row)
    stop("x and y have different number of observations.")
  ####################################################################
  hdsvm.object <- hdsvm(x, y, lambda = lambda, ...)
  lambda <- hdsvm.object$lambda
  nz <- sapply(coef(hdsvm.object, type = "nonzero"), length)
  if (missing(foldid))
    foldid <- sample(rep(seq(nfolds), length = x.row)) else nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be at least 3; nfolds = 5 recommended.")
  ## fit the model nfold times and save them
  dots <- list(...)
  outlist <- cv_folds(nfolds, function(i) {
    which <- foldid == i
    do.call(hdsvm, c(list(x = x[!which, , drop = FALSE], y = y[!which],
                          lambda = lambda), dots))
  }, ncores)
  ## select the lambda according to predmat
  cvstuff <- cvpath.hdsvm(outlist, x, y, lambda, foldid, x.row, ...)
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd,
              cvupper = cvm + cvsd, cvlower = cvm - cvsd, nzero = nz,
              name = cvname, hdsvm.fit = hdsvm.object)
  obj <- c(out, as.list(getmin(lambda, cvm, cvsd)))
  obj$call <- match.call()
  class(obj) <- "cv.hdsvm"
  obj
}

cvpath.hdsvm <- function(outlist, x, y, lambda, foldid, x.row, ...) {
  nfolds <- max(foldid)
  ## the fold models were fitted on the +/-1 coding; use the same coding here
  y <- c(-1, 1)[as.factor(y)]
  predmat <- matrix(NA, x.row, length(lambda))
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    whichfold <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[whichfold, , drop = FALSE], type = "loss")
    nlami <- length(fitobj$lambda)
    predmat[whichfold, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- svm_loss(y * predmat)
  N <- length(y) - apply(is.na(predmat), 2, sum)
  cvm <- colMeans(cvraw, na.rm = TRUE)
  scaled <- scale(cvraw, cvm, FALSE)^2
  cvsd <- sqrt(colMeans(scaled, na.rm = TRUE) / (N - 1))
  list(cvm = cvm, cvsd = cvsd, cvraw = cvraw, name = "Hinge loss")
}

svm_loss <- function(tval) pmax(1 - tval, 0)
