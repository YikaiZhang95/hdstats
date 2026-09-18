#' Plot the Cross-validation Curve
#'
#' Plots the cross-validation curve, with upper and lower standard-error
#' bars, as a function of \code{log(lambda)} for a cross-validated fit.
#' Vertical dotted lines mark \code{lambda.min} and \code{lambda.1se}, and
#' the number of nonzero coefficients is shown along the top axis.
#'
#' @param x A fitted \code{cv.hdhuber}, \code{cv.hdsvm}, \code{cv.hdqr},
#'   \code{cv.hdrr}, \code{cv.nc.hdsvm}, or \code{cv.nc.hdqr} object.
#' @param sign.lambda Either plot against \code{log(lambda)} (default) or
#'   its negative if \code{sign.lambda = -1}.
#' @param ... Other graphical parameters passed to \code{plot}.
#' @return No return value; called for its side effect of producing a plot.
#' @seealso \code{\link{cv.hdhuber}}, \code{\link{cv.hdsvm}},
#'   \code{\link{cv.hdqr}}, \code{\link{cv.hdrr}}
#' @name plot.cv
#' @examples
#' set.seed(315)
#' n <- 100
#' p <- 400
#' x <- matrix(data = rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
#' beta_star <- c(c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.3), rep(0, (p - 7)))
#' eps <- rnorm(n, mean = 0, sd = 1)
#' y <- x %*% beta_star + eps
#' cv.fit <- cv.hdqr(x = x, y = y, tau = 0.5)
#' plot(cv.fit)
NULL

plot_cv <- function(x, sign.lambda = 1, ...) {
  cvobj <- x
  xlab <- expression(log(lambda))
  if (sign.lambda < 0) xlab <- expression(-log(lambda))
  ylab <- if (is.null(cvobj$name)) "Cross-validation error" else cvobj$name
  plot.args <- list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm,
                    ylim = range(cvobj$cvupper, cvobj$cvlower),
                    xlab = xlab, ylab = ylab, type = "n")
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvupper, cvobj$cvlower,
             width = 0.01, col = "darkgrey")
  points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20, col = "red")
  axis(side = 3, at = sign.lambda * log(cvobj$lambda),
       labels = paste(cvobj$nzero), tick = FALSE, line = 0)
  abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
  abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
  invisible()
}

#' @rdname plot.cv
#' @method plot cv.hdhuber
#' @export
plot.cv.hdhuber <- function(x, sign.lambda = 1, ...) plot_cv(x, sign.lambda, ...)

#' @rdname plot.cv
#' @method plot cv.hdsvm
#' @export
plot.cv.hdsvm <- function(x, sign.lambda = 1, ...) plot_cv(x, sign.lambda, ...)

#' @rdname plot.cv
#' @method plot cv.hdqr
#' @export
plot.cv.hdqr <- function(x, sign.lambda = 1, ...) plot_cv(x, sign.lambda, ...)

#' @rdname plot.cv
#' @method plot cv.nc.hdsvm
#' @export
plot.cv.nc.hdsvm <- function(x, sign.lambda = 1, ...) plot_cv(x, sign.lambda, ...)

#' @rdname plot.cv
#' @method plot cv.nc.hdqr
#' @export
plot.cv.nc.hdqr <- function(x, sign.lambda = 1, ...) plot_cv(x, sign.lambda, ...)
