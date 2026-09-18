#' hdstats: High-Dimensional Penalized Huber, Quantile, Rank, and Support Vector Regression
#'
#' Fits the entire regularization path of elastic-net penalized Huber
#' regression, quantile regression, Wilcoxon rank regression, and support
#' vector machines for high-dimensional data, with optional nonconvex SCAD
#' and MCP penalties and k-fold cross-validation for tuning. All solvers use
#' coordinate descent with warm starts and active-set or strong-rule
#' screening; the non-smooth hinge and check losses are handled by the
#' finite smoothing algorithm of Tang, Zhang, and Wang (2024). The compute
#' kernels are written in C++ via \pkg{Rcpp}.
#'
#' @section Model fitting functions:
#' \describe{
#'   \item{\code{\link{hdhuber}}}{penalized Huber regression.}
#'   \item{\code{\link{hdqr}}}{penalized quantile regression.}
#'   \item{\code{\link{hdrr}}}{penalized Wilcoxon (canonical) rank
#'     regression, solved as median regression on pairwise differences.}
#'   \item{\code{\link{hdsvm}}}{penalized linear support vector machine.}
#'   \item{\code{\link{nc.hdqr}}, \code{\link{nc.hdsvm}}}{nonconvex (SCAD,
#'     MCP) penalized versions via local linear approximation.}
#' }
#' Each fitting function has a cross-validation counterpart
#' (\code{\link{cv.hdhuber}}, \code{\link{cv.hdqr}}, \code{\link{cv.hdrr}},
#' \code{\link{cv.hdsvm}}, \code{\link{cv.nc.hdqr}},
#' \code{\link{cv.nc.hdsvm}}), and all fitted objects have \code{coef},
#' \code{predict}, and \code{print} methods; cross-validation objects also
#' have a \code{plot} method.
#'
#' @references
#' Tang, Q., Zhang, Y., and Wang, B. (2024). Finite smoothing algorithm for
#' high-dimensional support vector machines and quantile regression.
#' \emph{Proceedings of the 41st International Conference on Machine
#' Learning (ICML)}. \url{https://openreview.net/forum?id=RvwMTDYTOb}
#'
#' Friedman, J., Hastie, T., and Tibshirani, R. (2010). Regularization
#' paths for generalized linear models via coordinate descent.
#' \emph{Journal of Statistical Software}, 33(1), 1--22.
#' \doi{10.18637/jss.v033.i01}
#'
#' @useDynLib hdstats, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import Matrix
#' @importFrom methods new rbind2
#' @importFrom stats approx coef predict median
#' @importFrom graphics segments points axis abline
#' @importFrom parallel mclapply makeCluster stopCluster parLapply
#' @keywords internal
"_PACKAGE"
