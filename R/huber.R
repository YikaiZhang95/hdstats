#' Fit Penalized Huber Regression
#'
#' Fits the solution path of elastic-net penalized Huber regression over a
#' grid of values of the tuning parameter \code{lambda}, using the finite
#' smoothing algorithm with coordinate descent, warm starts, and (optionally)
#' the strong rule for screening predictors.
#'
#' @param x Matrix of predictors, of dimension \eqn{n \times p}{n x p};
#'   each row is an observation.
#' @param y Numeric response vector of length \eqn{n}.
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
#' @param delta The Huber transition parameter \eqn{\delta > 0}: the loss is
#'   quadratic for residuals with \eqn{|r| \le \delta} and linear beyond.
#'   Default is 1.
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
#' @param is_strong Logical; if \code{TRUE} (the default) the strong rule
#'   is used to screen predictors at each \code{lambda} value, with KKT
#'   checks guaranteeing the same solution as without screening.
#'
#' @details
#' The penalized Huber regression problem solved is
#' \deqn{\min_{b_0, \beta} \frac{1}{n}\sum_{i=1}^n
#'   \ell_\delta(y_i - b_0 - x_i^\top\beta)
#'   + \lambda_1 |pf_1 \circ \beta|_1
#'   + \frac{\lambda_2}{2} \|\sqrt{pf_2} \circ \beta\|_2^2,}
#' where \eqn{\ell_\delta(r) = r^2/2} if \eqn{|r| \le \delta} and
#' \eqn{\ell_\delta(r) = \delta |r| - \delta^2/2} otherwise, and
#' \eqn{\circ} denotes the Hadamard product. The path is computed for a
#' decreasing sequence of \eqn{\lambda_1} values (argument \code{lambda})
#' with a single value of \eqn{\lambda_2} (argument \code{lam2}).
#'
#' The Huber loss is smooth, so it is minimized directly by coordinate
#' descent with a majorization step size. The compute kernel is written
#' in C++.
#'
#' @return
#' An object with S3 class \code{"hdhuber"} consisting of
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
#'   \item{npasses}{total number of coordinate-descent passes over the data
#'     for the whole path.}
#'   \item{jerr}{error flag, for warnings and errors; 0 if no error.}
#'   \item{delta}{the Huber parameter used.}
#'   \item{lam2}{the \eqn{\lambda_2} value used.}
#'
#' @seealso \code{\link{cv.hdhuber}}, \code{\link{coef.hdhuber}},
#'   \code{\link{predict.hdhuber}}, \code{\link{hdqr}}, \code{\link{hdsvm}}
#' @references
#' Tang, Q., Zhang, Y., and Wang, B. (2024). Finite smoothing algorithm for
#' high-dimensional support vector machines and quantile regression.
#' \emph{Proceedings of the 41st International Conference on Machine
#' Learning (ICML)}. \url{https://openreview.net/forum?id=RvwMTDYTOb}
#'
#' Huber, P. J. (1964). Robust estimation of a location parameter.
#' \emph{The Annals of Mathematical Statistics}, 35(1), 73--101.
#' \doi{10.1214/aoms/1177703732}
#' @keywords models regression robust
#' @export
#' @examples
#' set.seed(315)
#' n <- 100
#' p <- 400
#' x <- matrix(data = rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
#' beta_star <- c(c(2, 1.5, 0.8, 1, 1.75, 0.75, 0.3), rep(0, (p - 7)))
#' eps <- rnorm(n, mean = 0, sd = 1)
#' y <- x %*% beta_star + eps
#' fit <- hdhuber(x, y, delta = 0.5, lam2 = 0.01)
#' fit
#' coef(fit, s = fit$lambda[10])[1:8, ]
hdhuber <- function(x, y, nlambda = 100,
                    lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
                    lambda = NULL, lam2 = 0, delta = 1,
                    pf = rep(1, nvars), pf2 = rep(1, nvars),
                    exclude, dfmax = nvars + 1,
                    pmax = min(dfmax * 1.2, nvars), standardize = TRUE,
                    eps = 1e-08, maxit = 1e+06, is_strong = TRUE) {
  ####################################################################
  this.call <- match.call()
  y <- drop(y)
  x <- as.matrix(x)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)
  if (is.null(vnames))
    vnames <- paste0("V", seq(nvars))
  if (length(y) != nobs)
    stop("x and y have different number of observations")
  if (!is.numeric(y))
    stop("y must be numeric")
  if (!is.numeric(delta) || length(delta) != 1L || delta <= 0)
    stop("delta must be a single positive number")
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
  fit <- huber_cd_cpp(alpha, lam2, delta, nobs, nvars,
                      x, as.double(y), jd, pfncol, pf, pf2, dfmax,
                      pmax, nlam, flmin, ulam, eps, isd, maxit,
                      as.integer(is_strong))
  outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  fit <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  if (is.null(lambda))
    fit$lambda <- lamfix(fit$lambda)
  fit$delta <- delta
  fit$lam2 <- lam2
  fit$call <- this.call
  ####################################################################
  class(fit) <- c("hdhuber")
  fit
}
