#' Extract Model Coefficients from a \code{hdsvm} Object
#'
#' Retrieves the coefficients at specified values of \code{lambda} from a
#' fitted \code{\link{hdsvm}} object. If \code{s}, the vector of
#' \code{lambda} values, contains values not used in the model fitting,
#' linear interpolation between the closest fitted \code{lambda} values is
#' used.
#'
#' @param object Fitted \code{\link{hdsvm}} object.
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
#' @seealso \code{\link{hdsvm}}, \code{\link{predict.hdsvm}}
#' @method coef hdsvm
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
#' fit <- hdsvm(x, y, lam2 = 0.01)
#' coefs <- coef(fit, s = fit$lambda[3:5])
coef.hdsvm <- function(object, s = NULL, type = c("coefficients", "nonzero"), ...) {
  type <- match.arg(type)
  coef_path_type(object, s, type)
}

#' Make Predictions from a \code{hdsvm} Object
#'
#' Produces class labels or the linear predictor for new predictor data
#' using a fitted \code{\link{hdsvm}} object at specified \code{lambda}
#' values.
#'
#' @param object Fitted \code{\link{hdsvm}} object.
#' @param newx Matrix of new predictor values at which predictions are to be
#'   made. This is a required argument.
#' @param s Values of the penalty parameter \code{lambda} at which
#'   predictions are requested. Default is the entire sequence used during
#'   the model fit.
#' @param type Type \code{"class"} (the default) returns the predicted class
#'   labels, coded as \eqn{-1} (first level of \code{y}) and \eqn{+1}
#'   (second level); type \code{"loss"} returns the linear predictor
#'   \eqn{b_0 + x^\top\beta}.
#' @param ... Not used.
#' @return A matrix of predicted values, one row per row of \code{newx} and
#'   one column per value of \code{s}.
#' @seealso \code{\link{hdsvm}}, \code{\link{coef.hdsvm}}
#' @method predict hdsvm
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
#' fit <- hdsvm(x, y, lam2 = 0.01)
#' preds <- predict(fit, newx = tail(x), s = fit$lambda[3:5])
predict.hdsvm <- function(object, newx, s = NULL, type = c("class", "loss"), ...) {
  type <- match.arg(type)
  nfit <- predict_path(object, newx, s)
  switch(type, loss = nfit, class = ifelse(nfit > 0, 1, -1))
}
