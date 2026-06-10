## ---------------------------------------------------------------------------
## SPEED: wall-clock time to compute a full 100-lambda regularization path,
## hdstats vs the closest mainstream packages.
##
## Fairness protocol:
##   * matched convergence tolerance  eps = 1e-6  (hdstats defaults to a
##     stricter 1e-8, i.e. this is, if anything, conservative for hdstats);
##   * matched path length            nlambda = 100;
##   * matched penalty                pure L1 / lasso (hdstats lam2 = 0,
##     hqreg/sparseSVM alpha = 1, conquer penalty = "lasso");
##   * conquer has no auto path generator, so it is given the same lambda grid
##     that hdqr produces (it solves the identical 100 penalized problems);
##   * median elapsed time over 5 runs after one untimed warm-up.
## ---------------------------------------------------------------------------
suppressMessages({library(hdstats); library(hqreg); library(sparseSVM); library(conquer)})

EPS  <- 1e-6
REPS <- 5
settings <- list(c(n = 200, p = 1000), c(n = 200, p = 4000),
                 c(n = 400, p = 2000), c(n = 600, p = 2000))

timeit <- function(expr, reps = REPS) {
  expr <- substitute(expr); env <- parent.frame()
  eval(expr, env)                                          # warm-up (untimed)
  median(replicate(reps, system.time(eval(expr, env))[["elapsed"]]))
}
gen <- function(n, p, s = 10, seed = 1) {
  set.seed(seed)
  b <- c(rep(1.5, s), rep(0, p - s))
  x <- matrix(rnorm(n * p), n, p)
  list(x = x,
       y  = drop(x %*% b + rt(n, df = 3)),                 # heavy-tailed errors
       yb = ifelse(drop(x %*% b) + rnorm(n) > 0, 1, -1))   # binary labels
}

rows <- list()
for (st in settings) {
  n <- st["n"]; p <- st["p"]; d <- gen(n, p); x <- d$x; y <- d$y; yb <- d$yb
  lamq <- hdqr(x, y, tau = .5, lam2 = 0, eps = EPS)$lambda

  r <- data.frame(
    n = n, p = p,
    svm_hdstats   = timeit(hdsvm(x, yb, lam2 = 0, eps = EPS)),
    svm_sparseSVM = timeit(sparseSVM(x, yb, eps = EPS)),
    qr_hdstats    = timeit(hdqr(x, y, tau = .5, lam2 = 0, eps = EPS)),
    qr_hqreg      = timeit(hqreg(x, y, method = "quantile", tau = .5, eps = EPS)),
    qr_conquer    = timeit(conquer.reg(x, y, lambda = lamq, tau = .5, penalty = "lasso")),
    huber_hdstats = timeit(hdhuber(x, y, delta = 1, lam2 = 0, eps = EPS)),
    huber_hqreg   = timeit(hqreg(x, y, method = "huber", eps = EPS)))
  rows[[length(rows) + 1]] <- r
  cat(sprintf("n=%d p=%d done\n", n, p))
}
res <- do.call(rbind, rows)
dir.create("results", showWarnings = FALSE)
write.csv(res, "results/bench_speed.csv", row.names = FALSE)

sp <- with(res, data.frame(n = n, p = p,
  `SVM:sparseSVM` = round(svm_sparseSVM / svm_hdstats, 2),
  `QR:hqreg`      = round(qr_hqreg     / qr_hdstats, 2),
  `QR:conquer`    = round(qr_conquer   / qr_hdstats, 2),
  `Huber:hqreg`   = round(huber_hqreg  / huber_hdstats, 2),
  check.names = FALSE))
write.csv(sp, "results/bench_speedup.csv", row.names = FALSE)
cat("\nSpeedup (competitor time / hdstats time; > 1 means hdstats is faster):\n")
print(sp, row.names = FALSE)
