## Speed benchmark: time to compute a 100-value solution path with each
## package, on the same data, single-threaded. Run from the package root:
##   Rscript benchmarks/01-speed.R
## Environment variables: HDSTATS_BENCH_REPS (default 3), HDSTATS_BENCH_QUICK
## (any value: small sizes for a smoke test), HDSTATS_BENCH_LIB.
source("benchmarks/setup.R")
suppressPackageStartupMessages({ library(hqreg); library(sparseSVM); library(gcdnet); library(LiblineaR); library(Rfit) })

REPS <- as.integer(Sys.getenv("HDSTATS_BENCH_REPS", "3"))
QUICK <- nzchar(Sys.getenv("HDSTATS_BENCH_QUICK"))
SIZES <- if (QUICK) list(c(200, 500), c(400, 1000)) else
  list(c(200, 1000), c(500, 2000), c(1000, 5000), c(2000, 2000))
NLAMBDA <- 100
res <- list()
add <- function(model, package, method, n, p, seconds, note = "") {
  res[[length(res) + 1]] <<- data.frame(model = model, package = package, method = method,
                                        n = n, p = p, seconds = seconds, note = note,
                                        stringsAsFactors = FALSE)
  cat(sprintf("%-9s %-10s %-45s n=%5d p=%5d  %8.3f s\n", model, package, method, n, p, seconds))
}

for (sz in SIZES) {
  n <- sz[1]; p <- sz[2]
  lam.min <- if (n < p) 0.01 else 1e-4     # hdstats' default lambda.min ratio
  d <- sim_reg(n, p, s = 10, error = "t3", seed = 1)

  ## ---------------- Huber regression ----------------
  delta <- 1
  t <- time_it(fh <- hdhuber(d$x, d$y, delta = delta, nlambda = NLAMBDA), REPS)
  add("Huber", "hdstats", "hdhuber()", n, p, t)
  t <- time_it(hqreg(d$x, d$y, method = "huber", gamma = delta, nlambda = NLAMBDA, lambda.min = lam.min), REPS)
  add("Huber", "hqreg", "hqreg(method = 'huber')", n, p, t)

  ## ---------------- Quantile regression ----------------
  t <- time_it(fq <- hdqr(d$x, d$y, tau = 0.5, nlambda = NLAMBDA), REPS)
  add("Quantile", "hdstats", "hdqr()", n, p, t)
  if (n * p <= 1e6) {
    t <- time_it(hdqr(d$x, d$y, tau = 0.5, nlambda = NLAMBDA, is_exact = TRUE), REPS)
    add("Quantile", "hdstats", "hdqr(is_exact = TRUE)", n, p, t)
  }
  t <- time_it(hqreg(d$x, d$y, method = "quantile", tau = 0.5, nlambda = NLAMBDA, lambda.min = lam.min), REPS)
  add("Quantile", "hqreg", "hqreg(method = 'quantile')", n, p, t)
  if (has("conquer")) {
    t <- time_it(conquer::conquer.reg(d$x, d$y, lambda = fq$lambda, tau = 0.5), REPS)
    add("Quantile", "conquer", "conquer.reg() on hdqr's lambda sequence", n, p, t)
  }
  if (FALSE && has("rqPen")) {   # rqPen: alg = "huber" is hqreg, alg = "br" is quantreg; not built here
    t <- time_it(rqPen::rq.pen(d$x, d$y, tau = 0.5, nlambda = NLAMBDA, eps = lam.min, alg = "huber"), REPS)
    add("Quantile", "rqPen", "rq.pen(alg = 'huber')", n, p, t)
  }
  if (n * p <= 5e5) {   # LP solver: 10 lambda values, no warm starts
    lam10 <- fq$lambda[round(seq(1, NLAMBDA, length.out = 10))]
    t <- time_it(for (l in lam10) rq(d$y ~ d$x, tau = 0.5, method = "lasso", lambda = c(0, rep(n * l, p))), 1)
    add("Quantile", "quantreg", "rq(method = 'lasso'), 10 lambda values", n, p, t, "LP, 10 lambdas")
  }

  ## ---------------- SVM ----------------
  dc <- sim_class(n, p, s = 10, seed = 1)
  t <- time_it(fs <- hdsvm(dc$x, dc$y, nlambda = NLAMBDA), REPS)
  add("SVM", "hdstats", "hdsvm()", n, p, t)
  if (n * p <= 1e6) {
    t <- time_it(hdsvm(dc$x, dc$y, nlambda = NLAMBDA, is_exact = TRUE), REPS)
    add("SVM", "hdstats", "hdsvm(is_exact = TRUE)", n, p, t)
  }
  t <- time_it(sparseSVM(dc$x, dc$y, nlambda = NLAMBDA, lambda.min = lam.min), REPS)
  add("SVM", "sparseSVM", "sparseSVM()", n, p, t)
  t <- time_it(gcdnet(dc$x, dc$y, method = "hhsvm", nlambda = NLAMBDA), REPS)
  add("SVM", "gcdnet", "gcdnet(method = 'hhsvm'), Huberized squared hinge", n, p, t, "different loss")
  t <- time_it(LiblineaR(dc$x, dc$y, type = 5, cost = 1, bias = 1), REPS)
  add("SVM", "LiblineaR", "LiblineaR(type = 5), one cost value", n, p, t, "single fit, squared hinge")
}

## ---------------- Rank regression (n > p) ----------------
RSIZES <- if (QUICK) list(c(100, 10)) else list(c(100, 10), c(200, 20), c(400, 20))
for (sz in RSIZES) {
  n <- sz[1]; p <- sz[2]
  d <- sim_reg(n, p, s = min(10, p), error = "t3", seed = 2, rho = 0)
  t <- time_it(hdrr(d$x, d$y, nlambda = NLAMBDA), REPS)
  add("Rank", "hdstats", "hdrr(), 100-value path", n, p, t)
  t <- time_it(rfit(d$y ~ d$x), REPS)
  add("Rank", "Rfit", "rfit(), one unpenalized fit", n, p, t, "single fit")
  pr <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  xd <- d$x[pr[, 1], ] - d$x[pr[, 2], ]; yd <- d$y[pr[, 1]] - d$y[pr[, 2]]
  t <- time_it(rq(yd ~ xd - 1, tau = 0.5), REPS)
  add("Rank", "quantreg", "rq() on all pairwise differences, one fit", n, p, t, "single fit, LP")
}

res <- do.call(rbind, res)
res$platform <- paste(platform_info()[c("cpu")], collapse = "; ")
write.csv(res, "benchmarks/results/speed.csv", row.names = FALSE)
cat("\nwritten benchmarks/results/speed.csv\n")
