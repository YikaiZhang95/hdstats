#' Extract Coefficients from a \code{cv.nc.hdqr} Object
#'
#' Retrieves coefficients from a cross-validated \code{\link{nc.hdqr}}
#' model, using the stored full-data fit and the value of \code{lambda}
#' selected by cross-validation.
#'
#' @param object A fitted \code{\link{cv.nc.hdqr}} object.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#'   coefficients are desired. The default is \code{s = "lambda.1se"}, the
#'   largest value of \code{lambda} such that the cross-validation error is
#'   within one standard error of the minimum. Alternatively,
#'   \code{s = "lambda.min"} gives the value minimizing the cross-validation
#'   error. If \code{s} is numeric, these are taken as the actual values of
#'   \code{lambda} to use.
#' @param ... Not used.
#' @return A matrix of coefficients at the specified \code{lambda} values.
#' @seealso \code{\link{cv.nc.hdqr}}, \code{\link{predict.cv.nc.hdqr}}
#' @method coef cv.nc.hdqr
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
#'                         lam2 = 0.01)
#' coef(cv.nc.fit, s = c(0.02, 0.03))
#' }
coef.cv.nc.hdqr <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  coef(object$nchdqr.fit, s = cv_lambda(object, s), ...)
}

#' Make Predictions from a \code{cv.nc.hdqr} Object
#'
#' Generates predictions from a cross-validated \code{\link{nc.hdqr}}
#' model, using the stored full-data fit and the value of \code{lambda}
#' selected by cross-validation.
#'
#' @param object A fitted \code{\link{cv.nc.hdqr}} object.
#' @param newx Matrix of new predictor values at which predictions are to be
#'   made. This is a required argument.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#'   predictions are desired. The default is \code{s = "lambda.1se"}, the
#'   largest value of \code{lambda} such that the cross-validation error is
#'   within one standard error of the minimum. Alternatively,
#'   \code{s = "lambda.min"} gives the value minimizing the cross-validation
#'   error. If \code{s} is numeric, these are taken as the actual values of
#'   \code{lambda} to use.
#' @param ... Not used.
#' @return A matrix of predicted values, one row per row of \code{newx} and
#'   one column per value of \code{s}.
#' @seealso \code{\link{cv.nc.hdqr}}, \code{\link{coef.cv.nc.hdqr}}
#' @method predict cv.nc.hdqr
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
#'                         lam2 = 0.01)
#' predict(cv.nc.fit, newx = x[50:60, ], s = "lambda.min")
#' }
predict.cv.nc.hdqr <- function(object, newx, s = c("lambda.1se", "lambda.min"), ...) {
  predict(object$nchdqr.fit, newx, s = cv_lambda(object, s), ...)
}
