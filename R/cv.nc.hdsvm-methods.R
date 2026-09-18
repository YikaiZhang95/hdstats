#' Extract Coefficients from a \code{cv.nc.hdsvm} Object
#'
#' Retrieves coefficients from a cross-validated \code{\link{nc.hdsvm}}
#' model, using the stored full-data fit and the value of \code{lambda}
#' selected by cross-validation.
#'
#' @param object A fitted \code{\link{cv.nc.hdsvm}} object.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#'   coefficients are desired. The default is \code{s = "lambda.1se"}, the
#'   largest value of \code{lambda} such that the cross-validation error is
#'   within one standard error of the minimum. Alternatively,
#'   \code{s = "lambda.min"} gives the value minimizing the cross-validation
#'   error. If \code{s} is numeric, these are taken as the actual values of
#'   \code{lambda} to use.
#' @param ... Not used.
#' @return A matrix of coefficients at the specified \code{lambda} values.
#' @seealso \code{\link{cv.nc.hdsvm}}, \code{\link{predict.cv.nc.hdsvm}}
#' @method coef cv.nc.hdsvm
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
#' lambda <- 10^(seq(1, -4, length.out = 10))
#' \donttest{
#' cv.nc.fit <- cv.nc.hdsvm(x = x, y = y, lambda = lambda, lam2 = 0.01,
#'                          pen = "scad")
#' coef(cv.nc.fit, s = c(0.02, 0.03))
#' }
coef.cv.nc.hdsvm <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  coef(object$nchdsvm.fit, s = cv_lambda(object, s), ...)
}

#' Make Predictions from a \code{cv.nc.hdsvm} Object
#'
#' Generates predictions from a cross-validated \code{\link{nc.hdsvm}}
#' model, using the stored full-data fit and the value of \code{lambda}
#' selected by cross-validation.
#'
#' @param object A fitted \code{\link{cv.nc.hdsvm}} object.
#' @param newx Matrix of new predictor values at which predictions are to be
#'   made. This is a required argument.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#'   predictions are desired. The default is \code{s = "lambda.1se"}, the
#'   largest value of \code{lambda} such that the cross-validation error is
#'   within one standard error of the minimum. Alternatively,
#'   \code{s = "lambda.min"} gives the value minimizing the cross-validation
#'   error. If \code{s} is numeric, these are taken as the actual values of
#'   \code{lambda} to use.
#' @param type Type \code{"class"} (the default) returns the predicted class
#'   labels (\eqn{-1}/\eqn{+1}); type \code{"loss"} returns the linear
#'   predictor.
#' @param ... Not used.
#' @return A matrix of predicted values, one row per row of \code{newx} and
#'   one column per value of \code{s}.
#' @seealso \code{\link{cv.nc.hdsvm}}, \code{\link{coef.cv.nc.hdsvm}}
#' @method predict cv.nc.hdsvm
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
#' lambda <- 10^(seq(1, -4, length.out = 10))
#' \donttest{
#' cv.nc.fit <- cv.nc.hdsvm(x = x, y = y, lambda = lambda, lam2 = 0.01,
#'                          pen = "scad")
#' predict(cv.nc.fit, newx = x[50:60, ], s = "lambda.min")
#' }
predict.cv.nc.hdsvm <- function(object, newx, s = c("lambda.1se", "lambda.min"),
                                type = c("class", "loss"), ...) {
  type <- match.arg(type)
  predict(object$nchdsvm.fit, newx, s = cv_lambda(object, s), type = type, ...)
}
