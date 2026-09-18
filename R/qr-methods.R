#' Extract Model Coefficients from a \code{hdqr} Object
#'
#' Retrieves the coefficients at specified values of \code{lambda} from a
#' fitted \code{\link{hdqr}} (or \code{\link{hdrr}}) object. If \code{s},
#' the vector of \code{lambda} values, contains values not used in the model
#' fitting, linear interpolation between the closest fitted \code{lambda}
#' values is used.
#'
#' @param object Fitted \code{\link{hdqr}} or \code{\link{hdrr}} object.
#' @param s Values of the penalty parameter \code{lambda} at which
#'   coefficients are requested. Default is the entire sequence used during
#'   the model fit.
#' @param type Type \code{"coefficients"} (the default) returns the
#'   coefficients at the requested values of \code{s}; type \code{"nonzero"}
#'   returns a list of the indices of the nonzero coefficients for each value
#'   of \code{s}.
#' @param ... Not used.
#' @return A matrix of coefficients (intercept in the first row, one column
#'   per value of \code{s}), or a list of indices when
#'   \code{type = "nonzero"}.
#' @seealso \code{\link{hdqr}}, \code{\link{predict.hdqr}}
#' @method coef hdqr
#' @export
#' @examples
#' set.seed(315)
#' n <- 100
#' p <- 400
#' x <- matrix(data = rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
#' beta_star <- c(c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.3), rep(0, (p - 7)))
#' eps <- rnorm(n, mean = 0, sd = 1)
#' y <- x %*% beta_star + eps
#' fit <- hdqr(x = x, y = y, tau = 0.5, lam2 = 0.01)
#' coefs <- coef(fit, s = fit$lambda[3:5])
coef.hdqr <- function(object, s = NULL, type = c("coefficients", "nonzero"), ...) {
  type <- match.arg(type)
  coef_path_type(object, s, type)
}

#' Make Predictions from a \code{hdqr} Object
#'
#' Produces fitted values for new predictor data using a fitted
#' \code{\link{hdqr}} (or \code{\link{hdrr}}) object at specified
#' \code{lambda} values.
#'
#' @param object Fitted \code{\link{hdqr}} or \code{\link{hdrr}} object.
#' @param newx Matrix of new predictor values at which predictions are to be
#'   made. This is a required argument.
#' @param s Values of the penalty parameter \code{lambda} at which
#'   predictions are requested. Default is the entire sequence used during
#'   the model fit.
#' @param ... Not used.
#' @return A matrix of predicted values, one row per row of \code{newx} and
#'   one column per value of \code{s}.
#' @seealso \code{\link{hdqr}}, \code{\link{coef.hdqr}}
#' @method predict hdqr
#' @export
#' @examples
#' set.seed(315)
#' n <- 100
#' p <- 400
#' x <- matrix(data = rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
#' beta_star <- c(c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.3), rep(0, (p - 7)))
#' eps <- rnorm(n, mean = 0, sd = 1)
#' y <- x %*% beta_star + eps
#' fit <- hdqr(x = x, y = y, tau = 0.5, lam2 = 0.01)
#' preds <- predict(fit, newx = tail(x), s = fit$lambda[3:5])
predict.hdqr <- function(object, newx, s = NULL, ...) {
  predict_path(object, newx, s)
}
