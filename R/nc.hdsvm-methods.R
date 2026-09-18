#' Extract Model Coefficients from a \code{nc.hdsvm} Object
#'
#' Retrieves the coefficients at specified values of \code{lambda} from a
#' fitted \code{\link{nc.hdsvm}} object. If \code{s}, the vector of
#' \code{lambda} values, contains values not used in the model fitting,
#' linear interpolation between the closest fitted \code{lambda} values is
#' used.
#'
#' @param object Fitted \code{\link{nc.hdsvm}} object.
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
#' @seealso \code{\link{nc.hdsvm}}, \code{\link{predict.nc.hdsvm}}
#' @method coef nc.hdsvm
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
#' lambda <- 10^(seq(1, -4, length.out = 30))
#' \donttest{
#' nc.fit <- nc.hdsvm(x = x, y = y, lambda = lambda, lam2 = 0.01, pen = "scad")
#' nc.coefs <- coef(nc.fit, s = nc.fit$lambda[3:5])
#' }
coef.nc.hdsvm <- function(object, s = NULL, type = c("coefficients", "nonzero"), ...) {
  type <- match.arg(type)
  coef_path_type(object, s, type, lambda = object$nc.lambda)
}

#' Make Predictions from a \code{nc.hdsvm} Object
#'
#' Produces class labels or the linear predictor for new predictor data
#' using a fitted \code{\link{nc.hdsvm}} object at specified \code{lambda}
#' values.
#'
#' @param object Fitted \code{\link{nc.hdsvm}} object.
#' @param newx Matrix of new predictor values at which predictions are to be
#'   made. This is a required argument.
#' @param s Values of the penalty parameter \code{lambda} at which
#'   predictions are requested. Default is the entire sequence used during
#'   the model fit.
#' @param type Type \code{"class"} (the default) returns the predicted class
#'   labels (\eqn{-1}/\eqn{+1}); type \code{"loss"} returns the linear
#'   predictor.
#' @param ... Not used.
#' @return A matrix of predicted values, one row per row of \code{newx} and
#'   one column per value of \code{s}.
#' @seealso \code{\link{nc.hdsvm}}, \code{\link{coef.nc.hdsvm}}
#' @method predict nc.hdsvm
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
#' lambda <- 10^(seq(1, -4, length.out = 30))
#' \donttest{
#' nc.fit <- nc.hdsvm(x = x, y = y, lambda = lambda, lam2 = 0.01, pen = "scad")
#' nc.preds <- predict(nc.fit, newx = tail(x), s = nc.fit$lambda[3:5])
#' }
predict.nc.hdsvm <- function(object, newx, s = NULL, type = c("class", "loss"), ...) {
  type <- match.arg(type)
  nfit <- predict_path(object, newx, s, lambda = object$nc.lambda)
  switch(type, loss = nfit, class = ifelse(nfit > 0, 1, -1))
}
