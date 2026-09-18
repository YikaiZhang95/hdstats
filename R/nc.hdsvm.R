#' Fit Penalized SVM with Nonconvex Penalties
#'
#' Fits the penalized support vector machine with the nonconvex SCAD or MCP
#' penalty by the local linear approximation (LLA) algorithm: each LLA step
#' solves a weighted lasso-penalized SVM with \code{\link{hdsvm}}, using as
#' weights the derivative of the nonconvex penalty evaluated at the current
#' coefficients.
#'
#' @param x Matrix of predictors, of dimension \eqn{n \times p}{n x p};
#'   each row is an observation.
#' @param y Binary response of length \eqn{n}; see \code{\link{hdsvm}}.
#' @param lambda Sequence of \code{lambda} values for the nonconvex
#'   penalty. If \code{NULL} (the default), the \code{lambda} sequence
#'   computed by \code{\link{cv.hdsvm}} for the initial lasso fit is used.
#'   The sequence is sorted in decreasing order.
#' @param pen Type of nonconvex penalty: \code{"scad"} (default) or
#'   \code{"mcp"} (case-insensitive).
#' @param aval The concavity parameter of the SCAD or MCP penalty. Default
#'   is 3.7 for SCAD and 2 for MCP.
#' @param lam2 Regularization parameter \eqn{\lambda_2} for the quadratic
#'   penalty on the coefficients. Only one value of \code{lam2} is used per
#'   fit. Default is 1.
#' @param ini_beta Optional vector of initial coefficients (of length
#'   \eqn{p}, without intercept) to start the LLA iterations. If
#'   \code{NULL} (the default), the coefficients of a lasso-penalized SVM
#'   at \code{lambda.1se} selected by \code{\link{cv.hdsvm}} are used.
#' @param lla_step Number of LLA steps. Default is 3.
#' @param ... Additional arguments passed to \code{\link{hdsvm}}.
#'
#' @details
#' The SCAD penalty (Fan and Li, 2001) has derivative
#' \deqn{p'_\lambda(t) = \lambda I(t \le \lambda) +
#'   \frac{(a\lambda - t)_+}{a - 1} I(t > \lambda), \quad a > 2,}
#' and the MCP penalty (Zhang, 2010) has derivative
#' \deqn{p'_\lambda(t) = (\lambda - t/a)_+, \quad a > 1,}
#' for \eqn{t \ge 0}. At each LLA step (Zou and Li, 2008) the penalty is
#' linearized around the current estimate \eqn{\tilde\beta}, giving the
#' weighted lasso penalty \eqn{\sum_j p'_\lambda(|\tilde\beta_j|)|\beta_j|}
#' which is solved by \code{\link{hdsvm}} with penalty factors
#' \code{pf} equal to the weights.
#'
#' @return
#' An object with S3 class \code{"nc.hdsvm"} consisting of
#'   \item{call}{the call that produced this object.}
#'   \item{b0}{intercept sequence of length \code{length(lambda)}.}
#'   \item{beta}{a \eqn{p \times} \code{length(lambda)} matrix of
#'     coefficients, stored as a sparse matrix (\code{dgCMatrix} class).
#'     To convert it into an ordinary matrix, use \code{as.matrix()}.}
#'   \item{lambda}{the sequence of \code{lambda} values used for the
#'     nonconvex penalty (also stored as \code{nc.lambda}).}
#'   \item{df}{the number of nonzero coefficients for each value of
#'     \code{lambda}.}
#'   \item{npasses}{the number of coordinate-descent passes for every
#'     \code{lambda} value in the last LLA step.}
#'   \item{jerr}{error flag, for warnings and errors; 0 if no error.}
#'   \item{pen}{the penalty type used.}
#'   \item{aval}{the concavity parameter used.}
#' @seealso \code{\link{hdsvm}}, \code{\link{cv.nc.hdsvm}},
#'   \code{\link{coef.nc.hdsvm}}, \code{\link{predict.nc.hdsvm}}
#' @references
#' Fan, J. and Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. \emph{Journal of the American
#' Statistical Association}, 96(456), 1348--1360.
#' \doi{10.1198/016214501753382273}
#'
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. \emph{The Annals of Statistics}, 38(2), 894--942.
#' \doi{10.1214/09-AOS729}
#'
#' Zou, H. and Li, R. (2008). One-step sparse estimates in nonconcave
#' penalized likelihood models. \emph{The Annals of Statistics}, 36(4),
#' 1509--1533. \doi{10.1214/009053607000000802}
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
#' lambda <- 10^(seq(1, -4, length.out = 30))
#' \donttest{
#' nc.fit <- nc.hdsvm(x = x, y = y, lambda = lambda, lam2 = 0.01, pen = "scad")
#' nc.fit
#' }
nc.hdsvm <- function(x, y, lambda = NULL, pen = c("scad", "mcp"), aval = NULL,
                     lam2 = 1, ini_beta = NULL, lla_step = 3, ...) {
  this.call <- match.call()
  pen <- match.arg(tolower(pen), c("scad", "mcp"))
  if (pen == "scad") {
    if (is.null(aval)) aval <- 3.7
    if (!is.numeric(aval) || length(aval) != 1L || aval <= 2)
      stop("aval must be a single number greater than 2 for the SCAD penalty")
    pen_deriv <- deriv_scad
  } else {
    if (is.null(aval)) aval <- 2
    if (!is.numeric(aval) || length(aval) != 1L || aval <= 1)
      stop("aval must be a single number greater than 1 for the MCP penalty")
    pen_deriv <- deriv_mcp
  }
  if (!is.numeric(lla_step) || length(lla_step) != 1L || lla_step < 1)
    stop("lla_step must be a positive integer")
  if (!is.null(ini_beta) && length(ini_beta) != NCOL(x)) {
    warning("wrong length of ini_beta; it is ignored")
    ini_beta <- NULL
  }
  if (is.null(lambda)) {
    if (!is.null(ini_beta))
      stop("lambda must be provided with ini_beta.")
    fit <- cv.hdsvm(x, y, lam2 = lam2, ...)
    lambda <- fit$lambda
    ini_beta <- coef(fit, s = "lambda.1se")[-1]
  } else {
    if (any(lambda < 0))
      stop("The values of lambda should be non-negative.")
    lambda <- rev(sort(as.double(lambda)))
    if (is.null(ini_beta)) {
      fit <- cv.hdsvm(x, y, lambda = lambda, lam2 = lam2, ...)
      ini_beta <- coef(fit, s = "lambda.1se")[-1]
    }
  }
  nlam <- length(lambda)
  pfmat <- matrix(NA_real_, NCOL(x), nlam)
  beta <- matrix(ini_beta, length(ini_beta), nlam)
  for (j in seq(lla_step)) {
    for (l in seq(nlam)) {
      pfmat[, l] <- as.vector(pen_deriv(abs(beta[, l, drop = FALSE]), lambda[l], aval))
    }
    fit <- hdsvm(x, y, lambda = rep(1, nlam), lam2 = lam2, pf = pfmat, ...)
    beta <- coef(fit)[-1, , drop = FALSE]
  }
  fit$lambda <- lambda
  fit$nc.lambda <- lambda
  fit$pen <- pen
  fit$aval <- aval
  fit$call <- this.call
  class(fit) <- c("nc.hdsvm")
  fit
}
