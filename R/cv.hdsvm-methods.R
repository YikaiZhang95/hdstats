#' Extract Coefficients from a \code{cv.hdsvm} Object
#'
#' Retrieves coefficients from a cross-validated \code{\link{hdsvm}} model,
#' using the stored full-data fit and the value of \code{lambda} selected by
#' cross-validation.
#'
#' @param object A fitted \code{\link{cv.hdsvm}} object.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#'   coefficients are desired. The default is \code{s = "lambda.1se"}, the
#'   largest value of \code{lambda} such that the cross-validation error is
#'   within one standard error of the minimum. Alternatively,
#'   \code{s = "lambda.min"} gives the value minimizing the cross-validation
#'   error. If \code{s} is numeric, these are taken as the actual values of
#'   \code{lambda} to use.
#' @param ... Not used.
#' @return A matrix of coefficients at the specified \code{lambda} values.
#' @seealso \code{\link{cv.hdsvm}}, \code{\link{predict.cv.hdsvm}}
#' @method coef cv.hdsvm
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
#' coef(cv.fit, s = c(0.02, 0.03))
coef.cv.hdsvm <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  coef(object$hdsvm.fit, s = cv_lambda(object, s), ...)
}

#' Make Predictions from a \code{cv.hdsvm} Object
#'
#' Generates predictions from a cross-validated \code{\link{hdsvm}} model,
#' using the stored full-data fit and the value of \code{lambda} selected by
#' cross-validation.
#'
#' @param object A fitted \code{\link{cv.hdsvm}} object.
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
#' @seealso \code{\link{cv.hdsvm}}, \code{\link{coef.cv.hdsvm}}
#' @method predict cv.hdsvm
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
#' predict(cv.fit, newx = x[50:60, ], s = "lambda.min")
predict.cv.hdsvm <- function(object, newx, s = c("lambda.1se", "lambda.min"),
                             type = c("class", "loss"), ...) {
  type <- match.arg(type)
  predict(object$hdsvm.fit, newx, s = cv_lambda(object, s), type = type, ...)
}
