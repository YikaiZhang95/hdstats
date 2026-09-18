#' Print a Fitted Model or Cross-validation Object
#'
#' For a fitted solution path (\code{hdhuber}, \code{hdsvm}, \code{hdqr},
#' \code{hdrr}, \code{nc.hdsvm}, or \code{nc.hdqr} object), prints the call
#' and a table with the number of nonzero coefficients (\code{Df}) at each
#' \code{lambda}. For a cross-validation object, prints the call, the error
#' measure, and a summary of the fit at \code{lambda.min} and
#' \code{lambda.1se}.
#'
#' @param x A fitted model or cross-validation object.
#' @param digits Number of significant digits in the printout.
#' @param ... Not used.
#' @return The object \code{x}, invisibly.
#' @name print.hdstats
#' @examples
#' set.seed(1)
#' n <- 50; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] * 2 - x[, 2] * 1.5 + rnorm(n)
#' fit <- hdhuber(x, y, nlambda = 10)
#' print(fit)
#' cv.fit <- cv.hdhuber(x, y, nlambda = 10)
#' print(cv.fit)
NULL

print_path <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  lambda <- if (!is.null(x$nc.lambda)) x$nc.lambda else x$lambda
  print(cbind(Df = x$df, Lambda = signif(lambda, digits)))
  invisible(x)
}

print_cv <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  cat("Measure:", x$name, "\n\n")
  optlams <- c(x$lambda.min, x$lambda.1se)
  which <- match(optlams, x$lambda)
  mat <- data.frame(Lambda = signif(optlams, digits), Index = which,
                    Measure = signif(x$cvm[which], digits),
                    SE = signif(x$cvsd[which], digits),
                    Nonzero = x$nzero[which])
  rownames(mat) <- c("min", "1se")
  print(mat)
  invisible(x)
}

#' @rdname print.hdstats
#' @method print hdhuber
#' @export
print.hdhuber <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_path(x, digits, ...)

#' @rdname print.hdstats
#' @method print hdsvm
#' @export
print.hdsvm <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_path(x, digits, ...)

#' @rdname print.hdstats
#' @method print hdqr
#' @export
print.hdqr <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_path(x, digits, ...)

#' @rdname print.hdstats
#' @method print nc.hdsvm
#' @export
print.nc.hdsvm <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_path(x, digits, ...)

#' @rdname print.hdstats
#' @method print nc.hdqr
#' @export
print.nc.hdqr <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_path(x, digits, ...)

#' @rdname print.hdstats
#' @method print cv.hdhuber
#' @export
print.cv.hdhuber <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_cv(x, digits, ...)

#' @rdname print.hdstats
#' @method print cv.hdsvm
#' @export
print.cv.hdsvm <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_cv(x, digits, ...)

#' @rdname print.hdstats
#' @method print cv.hdqr
#' @export
print.cv.hdqr <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_cv(x, digits, ...)

#' @rdname print.hdstats
#' @method print cv.nc.hdsvm
#' @export
print.cv.nc.hdsvm <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_cv(x, digits, ...)

#' @rdname print.hdstats
#' @method print cv.nc.hdqr
#' @export
print.cv.nc.hdqr <- function(x, digits = max(3, getOption("digits") - 3), ...)
  print_cv(x, digits, ...)
