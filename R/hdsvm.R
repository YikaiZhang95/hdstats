#' Fit Penalized Support Vector Machines
#'
#' Fits the solution path of the elastic-net penalized linear support vector
#' machine (SVM) over a grid of values of the tuning parameter \code{lambda},
#' using the finite smoothing algorithm with coordinate descent, warm
#' starts, and active-set updates.
#'
#' @param x Matrix of predictors, of dimension \eqn{n \times p}{n x p};
#'   each row is an observation.
#' @param y Binary response of length \eqn{n}: either a factor with two
#'   levels or a vector with exactly two distinct values (e.g. \code{-1}/\code{1}
#'   or \code{0}/\code{1}). Internally the first level is coded as \eqn{-1}
#'   and the second as \eqn{+1}; class predictions are returned in this
#'   coding.
#' @param nlambda The number of \code{lambda} values (default is 100).
#' @param lambda.factor The factor for getting the minimal value
#'   in the \code{lambda} sequence, where
#'   \code{min(lambda) = lambda.factor * max(lambda)}
#'   and \code{max(lambda)} is the smallest value of \code{lambda}
#'   for which all coefficients (except the intercept)
#'   are penalized to zero. The default depends on the relationship
#'   between \eqn{n} (the number of rows in the design matrix) and
#'   \eqn{p} (the number of predictors): if \eqn{n < p}, it is
#'   \code{0.01}; otherwise it is \code{1e-04}. A very small value of
#'   \code{lambda.factor} will lead to a saturated fit. The argument has no
#'   effect if a \code{lambda} sequence is supplied by the user.
#' @param lambda A user-supplied \code{lambda} sequence. Typically,
#'   by leaving this option unspecified, users can have the program
#'   compute its own \code{lambda} sequence based on \code{nlambda}
#'   and \code{lambda.factor}. It is better to supply, if necessary,
#'   a decreasing sequence of \code{lambda} values than a single
#'   (small) value. The program will ensure that the user-supplied
#'   \code{lambda} sequence is sorted in decreasing order before
#'   fitting the model to take advantage of the warm-start technique.
#' @param lam2 Regularization parameter \eqn{\lambda_2} for the
#'   quadratic (ridge) penalty on the coefficients. Unlike \code{lambda},
#'   only one value of \code{lam2} is used for each fitting process.
#'   Default is 0 (lasso penalty).
#' @param hval Initial bandwidth of the uniform-kernel convolution smoothing
#'   of the hinge loss. Default is 1. Inside the algorithm the bandwidth is
#'   divided by 8 (at most three times) until the solution satisfies the
#'   KKT conditions of the original (non-smooth) problem to a fixed
#'   tolerance; see \code{is_exact}.
#' @param pf L1 penalty factor of length \eqn{p} used for the adaptive
#'   lasso or adaptive elastic net. Separate L1 penalty weights can be
#'   applied to each coefficient to allow different L1 shrinkage.
#'   Can be 0 for some variables (but not all), which imposes no
#'   shrinkage, and results in that variable always being included
#'   in the model. Default is 1 for all variables (and implicitly
#'   infinity for variables in the \code{exclude} list). A matrix with
#'   \code{nlambda} columns (one penalty-factor vector per \code{lambda}
#'   value) is also accepted.
#' @param pf2 L2 penalty factor of length \eqn{p} used for adaptive
#'   elastic net. Separate L2 penalty weights can be applied to
#'   each coefficient to allow different L2 shrinkage.
#'   Can be 0 for some variables, which imposes no shrinkage.
#'   Default is 1 for all variables.
#' @param exclude Indices of variables to be excluded from the model.
#'   Default is none. Equivalent to an infinite penalty factor.
#' @param dfmax The maximum number of variables allowed in the model.
#'   Useful for very large \eqn{p} when a partial path is desired.
#'   Default is \eqn{p+1}.
#' @param pmax The maximum number of coefficients allowed ever
#'   to be nonzero along the solution path. For example, once
#'   \eqn{\beta_j} enters the model, no matter how many times it
#'   exits or re-enters the model through the path, it will be
#'   counted only once. Default is \code{min(dfmax * 1.2, p)}.
#' @param standardize Logical flag for variable standardization,
#'   prior to fitting the model sequence. The coefficients are
#'   always returned on the original scale. Default is \code{TRUE}.
#' @param eps Convergence threshold for coordinate descent: iterations
#'   stop when the largest squared change of any coefficient within one
#'   pass falls below \code{eps}. Default is \code{1e-08}.
#' @param maxit Maximum number of coordinate-descent passes over the data,
#'   summed over the whole \code{lambda} path. Default is \code{1e+06}.
#' @param sigma Penalty parameter of the quadratic term in the augmented
#'   Lagrangian used by the exact projection step (only used when
#'   \code{is_exact = TRUE}). Must be positive. Default is 0.9.
#' @param is_exact Logical. If \code{FALSE} (the default), the minimizer of
#'   the smoothed problem at the last bandwidth is returned; it satisfies
#'   the KKT conditions of the hinge-loss problem approximately and its
#'   hinge-loss objective is typically within a few percent of the optimum.
#'   If \code{TRUE}, an additional ADMM-type projection step is applied at
#'   small bandwidths, which yields the exact solution of the hinge-loss
#'   problem (for \code{lambda = 0} it agrees with the solution of the
#'   standard linear SVM quadratic program) at a higher computational cost.
#'
#' @details
#' The penalized SVM problem solved is
#' \deqn{\min_{b_0, \beta} \frac{1}{n}\sum_{i=1}^n
#'   \max\{1 - y_i (b_0 + x_i^\top \beta), 0\}
#'   + \lambda_1 |pf_1 \circ \beta|_1
#'   + \frac{\lambda_2}{2} \|\sqrt{pf_2} \circ \beta\|_2^2,}
#' where \eqn{\circ} denotes the Hadamard product. The hinge loss is
#' replaced by its convolution with a uniform kernel of bandwidth
#' \code{hval}, which yields a smooth surrogate that is minimized by
#' coordinate descent; the bandwidth is then shrunk and the fit is
#' refined (the finite smoothing algorithm of Tang, Zhang, and Wang, 2024).
#'
#' For faster computation, if the algorithm is not converging or
#' running slowly, consider increasing \code{eps}, increasing
#' \code{sigma}, decreasing \code{nlambda}, or increasing
#' \code{lambda.factor} before increasing \code{maxit}.
#'
#' @return
#' An object with S3 class \code{"hdsvm"} consisting of
#'   \item{call}{the call that produced this object.}
#'   \item{b0}{intercept sequence of length \code{length(lambda)}.}
#'   \item{beta}{a \eqn{p \times} \code{length(lambda)} matrix of
#'     coefficients, stored as a sparse matrix (\code{dgCMatrix} class,
#'     the standard class for sparse numeric matrices in the \code{Matrix}
#'     package). To convert it into an ordinary matrix, use
#'     \code{as.matrix()}.}
#'   \item{lambda}{the actual sequence of \code{lambda} values used.}
#'   \item{df}{the number of nonzero coefficients for each value of
#'     \code{lambda}.}
#'   \item{dim}{dimension of the coefficient matrix.}
#'   \item{npasses}{the number of coordinate-descent passes for every
#'     \code{lambda} value.}
#'   \item{jerr}{error flag, for warnings and errors; 0 if no error.}
#'   \item{hval}{the initial smoothing bandwidth used.}
#'   \item{lam2}{the \eqn{\lambda_2} value used.}
#' @seealso \code{\link{cv.hdsvm}}, \code{\link{nc.hdsvm}},
#'   \code{\link{coef.hdsvm}}, \code{\link{predict.hdsvm}}
#' @references
#' Tang, Q., Zhang, Y., and Wang, B. (2024). Finite smoothing algorithm for
#' high-dimensional support vector machines and quantile regression.
#' \emph{Proceedings of the 41st International Conference on Machine
#' Learning (ICML)}. \url{https://openreview.net/forum?id=RvwMTDYTOb}
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
#' fit <- hdsvm(x, y, lam2 = 0.01)
#' fit
#' predict(fit, newx = x[1:5, ], s = fit$lambda[10])
hdsvm <- function(x, y, nlambda = 100,
                  lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
                  lambda = NULL, lam2 = 0, hval = 1,
                  pf = rep(1, nvars), pf2 = rep(1, nvars),
                  exclude, dfmax = nvars + 1,
                  pmax = min(dfmax * 1.2, nvars), standardize = TRUE,
                  eps = 1e-08, maxit = 1e+06, sigma = 0.9, is_exact = FALSE) {
  ####################################################################
  this.call <- match.call()
  y <- drop(y)
  if (nlevels(as.factor(y)) != 2L)
    stop("y should be a factor with two levels (or a vector with exactly two distinct values).")
  y <- c(-1, 1)[as.factor(y)]
  x <- as.matrix(x)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)
  if (is.null(vnames))
    vnames <- paste0("V", seq(nvars))
  if (length(y) != nobs)
    stop("x and y have different number of observations")
  if (!is.numeric(hval) || length(hval) != 1L || hval <= 0)
    stop("hval must be a single positive number")
  ####################################################################
  ## parameter setup
  alpha <- NULL # user can change this to tune elastic-net based on alpha.
  if (!is.null(alpha)) {
    alpha <- as.double(alpha)
    if (alpha <= 0 || alpha > 1)
      stop("alpha: 0 < alpha <= 1")
    if (!is.null(lam2))
      warning("lam2 has no effect.")
    lam2 <- -1.0
  } else {
    if (!is.null(lam2)) {
      if (lam2 < 0) stop("lam2 is non-negative.")
      alpha <- -1.0
    } else {
      alpha <- 1.0 # default lasso
    }
  }
  maxit <- as.integer(maxit)
  isd <- as.integer(standardize)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  if (!missing(exclude)) {
    jd <- match(exclude, seq(nvars), 0)
    if (!all(jd > 0))
      stop("Some excluded variables are out of range.")
    jd <- as.integer(c(length(jd), jd))
  } else jd <- as.integer(0)
  ####################################################################
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1)
      stop("lambda.factor should be less than 1.")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    ## flmin = 1 if user defines lambda
    flmin <- as.double(1)
    if (any(lambda < 0))
      stop("The values of lambda should be non-negative.")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  pfncol <- NCOL(pf)
  pf <- matrix(as.double(pf), ncol = pfncol)
  if (NROW(pf) != nvars)
    stop("The size of L1 penalty factor must be the same with the number of input variables.")
  if (pfncol != 1 && pfncol != nlam)
    stop("pf should be a matrix with 1 or length(lambda) columns.")
  if (length(pf2) != nvars)
    stop("The size of L2 penalty factor must be the same with the number of input variables.")
  pf2 <- as.double(pf2)
  ####################################################################
  fit <- hdsvm_cd_cpp(alpha, lam2, hval, nobs, nvars,
                      x, as.double(y), jd, pfncol, pf, pf2, dfmax,
                      pmax, nlam, flmin, ulam, eps, isd, maxit,
                      as.double(sigma), as.integer(is_exact))
  outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  fit <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  if (is.null(lambda))
    fit$lambda <- lamfix(fit$lambda)
  fit$hval <- hval
  fit$lam2 <- lam2
  fit$call <- this.call
  ####################################################################
  class(fit) <- c("hdsvm")
  fit
}
