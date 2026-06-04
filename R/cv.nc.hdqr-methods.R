#' Extract Coefficients from a \code{cv.nc.hdqr} Object
#'
#' Retrieves coefficients at specified values of \code{lambda} from a fitted \code{cv.nc.hdqr} model.
#' Utilizes the stored \code{nchdqr.fit} object and the optimal \code{lambda} values determined during
#' the cross-validation process.
#'
#' @param object A fitted \code{cv.nc.hdqr} object from which coefficients are to be extracted.
#' @param s Specifies the \code{lambda} values at which coefficients are requested.
#'   The default is \code{s = "lambda.1se"}, representing the largest \code{lambda} such that the cross-validation
#'   error estimate is within one standard error of the minimum. Alternatively, \code{s = "lambda.min"}
#'   corresponds to the \code{lambda} yielding the minimum cross-validation error. If \code{s} is numeric, these
#'   values are directly used as the \code{lambda} values for coefficient extraction.
#' @param ... Not used.
#'
#' @return Returns a vector or matrix of coefficients corresponding to the specified `lambda` values.
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
#' tau <- 0.5
#' lam2 <- 0.01
#' lambda <- 10^(seq(1,-4, length.out=30))
#' \donttest{cv.nc.fit <- cv.nc.hdqr(x = x, y = y, tau = tau, lambda = lambda, lam2 = lam2)}
#' \donttest{coef(cv.nc.fit, s = c(0.02, 0.03))}


coef.cv.nc.hdqr <- function(object,
                            s = c("lambda.1se", "lambda.min"),
                            ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  } else stop("Invalid form for s")
  coef(object$nchdqr.fit, s = lambda, ...)
}

#' Make Predictions from a \code{cv.nc.hdqr} Object
#'
#' Generates predictions using a fitted \code{cv.nc.hdqr} object. This function utilizes the
#' stored \code{nchdqr.fit} object and an optimal value of \code{lambda} determined during the
#' cross-validation process.
#'
#' @param newx Matrix of new predictor values for which predictions are desired.
#'   This must be a matrix and is a required argument.
#' @param object A fitted \code{cv.nc.hdqr} object from which predictions are to be made.
#' @param s Specifies the value(s) of the penalty parameter \code{lambda} at which predictions
#'   are desired. The default is \code{s = "lambda.1se"}, representing the largest value of \code{lambda}
#'   such that the cross-validation error estimate is within one standard error of the minimum.
#'   Alternatively, \code{s = "lambda.min"} can be used, corresponding to the minimum of the
#'   cross-validation error estimate. If \code{s} is numeric, these are taken as the actual values
#'   of \code{lambda} to use for predictions.
#' @param ... Not used.
#' @return Returns a matrix or vector of predicted values corresponding to the specified
#'   `lambda` values.
#' @seealso \code{\link{cv.nc.hdqr}}, \code{\link{predict.cv.nc.hdqr}}
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
#' tau <- 0.5
#' lam2 <- 0.01
#' lambda <- 10^(seq(1,-4, length.out=10))
#' \donttest{cv.nc.fit <- cv.nc.hdqr(x = x, y = y, tau = tau, lambda = lambda, lam2 = lam2)}
#' \donttest{predict(cv.nc.fit, newx = x[50:60, ], s = "lambda.min")}

predict.cv.nc.hdqr <- function(object, newx,
                               s = c("lambda.1se", "lambda.min"),
                               ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  } else stop("Invalid form for s")
    predict(object$nchdqr.fit, newx, s = lambda, ...)
}