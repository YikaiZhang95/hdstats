## Optimization quality: objective value attained by each solver for the
## same penalized problem (same lambda values, no standardization), compared
## with the exact optimum. Run from the package root:
##   Rscript benchmarks/02-quality.R
##
## Exact references: for the Huber loss all solvers converge to the same
## optimum; for the quantile and hinge losses the problem is a linear program
## solved with quantreg::rq(method = "lasso"). Note that quantreg's lasso
## penalty is lambda/2 * sum|beta| on the summed loss, so lambda_rq = 2 n lambda
## reproduces (1/n) sum loss + lambda sum|beta|. The hinge loss max(1 - m, 0)
## is the check loss rho_tau(1 - m) at tau = 1 (we use tau = 1 - 1e-6).
source("benchmarks/setup.R")
suppressPackageStartupMessages({ library(hqreg); library(sparseSVM) })
coef <- stats::coef

idx <- c(10, 30, 50, 70, 90)      # positions along a 100-value path
res <- list()
add <- function(regime, model, solver, lambda, obj, df, secs)
  res[[length(res) + 1]] <<- data.frame(regime = regime, model = model, solver = solver, lambda = lambda,
                                        objective = obj, df = df, seconds = secs, stringsAsFactors = FALSE)
run_all <- function(regime, model, fits, lam, objfun, extract) {
  for (nm in names(fits)) {
    f <- tryCatch(fits[[nm]](), error = function(e) { message(nm, ": ", conditionMessage(e)); NULL })
    if (is.null(f)) next
    secs <- system.time(fits[[nm]]())[["elapsed"]]
    cf <- extract(f)
    for (k in seq_along(lam))
      add(regime, model, nm, lam[k], objfun(cf[1, k], cf[-1, k], lam[k]), sum(abs(cf[-1, k]) > 1e-8), secs)
  }
}
ext_hd <- function(f) as.matrix(coef(f))
ext_beta <- function(f) as.matrix(f$beta)

for (sz in list(c(400, 200), c(200, 500))) {
  n <- sz[1]; p <- sz[2]
  regime <- sprintf("n = %d, p = %d", n, p)
  d <- sim_reg(n, p, s = 10, error = "t3", seed = 3)
  dc <- sim_class(n, p, s = 10, seed = 3)

  ## ---------------- Huber (delta = 1): all solvers target the same optimum ----------------
  delta <- 1
  ref <- hdhuber(d$x, d$y, delta = delta, standardize = FALSE, nlambda = 100)
  lam <- ref$lambda[idx]
  run_all(regime, "Huber", list(
    "hdstats::hdhuber (default)" = function() hdhuber(d$x, d$y, delta = delta, lambda = lam, standardize = FALSE),
    "hdstats::hdhuber (eps = 1e-12)" = function() hdhuber(d$x, d$y, delta = delta, lambda = lam, standardize = FALSE, eps = 1e-12),
    "hqreg::hqreg_raw (default)" = function() { f <- hqreg_raw(d$x, d$y, method = "huber", gamma = delta, lambda = lam / delta); list(beta = f$beta) },
    "hqreg::hqreg_raw (eps = 1e-10)" = function() { f <- hqreg_raw(d$x, d$y, method = "huber", gamma = delta, lambda = lam / delta, eps = 1e-10); list(beta = f$beta) }),
    lam, function(b0, b, l) obj_huber(b0, b, d$x, d$y, delta, l), function(f) if (inherits(f, "hdhuber")) ext_hd(f) else ext_beta(f))

  ## ---------------- Quantile regression: LP reference ----------------
  for (tau in c(0.5, 0.3)) {
    ref <- hdqr(d$x, d$y, tau = tau, lam2 = 0, standardize = FALSE, nlambda = 100)
    lam <- ref$lambda[idx]
    fits <- list(
      "quantreg::rq lasso (LP, exact)" = function() list(beta = sapply(lam, function(l)
        coef(rq(d$y ~ d$x, tau = tau, method = "lasso", lambda = c(0, rep(2 * n * l, p)))))),
      "hdstats::hdqr (default)" = function() hdqr(d$x, d$y, tau = tau, lambda = lam, lam2 = 0, standardize = FALSE),
      "hdstats::hdqr (is_exact = TRUE)" = function() hdqr(d$x, d$y, tau = tau, lambda = lam, lam2 = 0, standardize = FALSE, is_exact = TRUE),
      "hqreg::hqreg_raw (default)" = function() { f <- hqreg_raw(d$x, d$y, method = "quantile", tau = tau, lambda = lam); list(beta = f$beta) },
      "hqreg::hqreg_raw (eps = 1e-10)" = function() { f <- hqreg_raw(d$x, d$y, method = "quantile", tau = tau, lambda = lam, eps = 1e-10); list(beta = f$beta) })
    if (has("conquer")) {
      cq <- function(...) { f <- conquer::conquer.reg(d$x, d$y, lambda = lam, tau = tau, ...); list(beta = f$coeff[, match(lam, sort(lam)), drop = FALSE]) }
      fits[["conquer::conquer.reg (Gaussian kernel, default)"]] <- function() cq()
      fits[["conquer::conquer.reg (uniform kernel, epsilon = 1e-6)"]] <- function() cq(kernel = "uniform", epsilon = 1e-6, iteMax = 5000)
    }
    run_all(regime, paste0("Quantile (tau = ", tau, ")"), fits, lam,
            function(b0, b, l) obj_qr(b0, b, d$x, d$y, tau, l), function(f) if (inherits(f, "hdqr")) ext_hd(f) else ext_beta(f))
  }

  ## ---------------- SVM: LP reference ----------------
  ref <- hdsvm(dc$x, dc$y, lam2 = 0, standardize = FALSE, nlambda = 100)
  lam <- ref$lambda[idx]
  z <- cbind(dc$y, dc$y * dc$x)
  fits <- list(
    "quantreg::rq lasso (LP, exact)" = function() list(beta = sapply(lam, function(l)
      coef(rq(rep(1, n) ~ z - 1, tau = 1 - 1e-6, method = "lasso", lambda = c(0, rep(2 * n * l, p)))))),
    "hdstats::hdsvm (default, hval = 1)" = function() hdsvm(dc$x, dc$y, lambda = lam, lam2 = 0, standardize = FALSE),
    "hdstats::hdsvm (hval = 0.1)" = function() hdsvm(dc$x, dc$y, lambda = lam, lam2 = 0, standardize = FALSE, hval = 0.1),
    "hdstats::hdsvm (is_exact = TRUE)" = function() hdsvm(dc$x, dc$y, lambda = lam, lam2 = 0, standardize = FALSE, is_exact = TRUE),
    ## sparseSVM codes the first level of y as +1, hence the sign flip of its weights
    "sparseSVM::sparseSVM (gamma = 0.1, default)" = function() { f <- sparseSVM(dc$x, dc$y, lambda = lam, preprocess = "none"); list(beta = -f$weights) },
    "sparseSVM::sparseSVM (gamma = 0.01, eps = 1e-8)" = function() { f <- sparseSVM(dc$x, dc$y, lambda = lam, preprocess = "none", gamma = 0.01, eps = 1e-8, max.iter = 1e5); list(beta = -f$weights) })
  run_all(regime, "SVM", fits, lam, function(b0, b, l) obj_svm(b0, b, dc$x, dc$y, l), function(f) if (inherits(f, "hdsvm")) ext_hd(f) else ext_beta(f))
}

res <- do.call(rbind, res)
best <- ave(res$objective, res$regime, res$model, res$lambda, FUN = min)
res$rel_gap <- (res$objective - best) / best
write.csv(res, "benchmarks/results/quality.csv", row.names = FALSE)
summ <- do.call(rbind, lapply(split(res, list(res$regime, res$model, res$solver), drop = TRUE), function(g)
  data.frame(regime = g$regime[1], model = g$model[1], solver = g$solver[1],
             max_rel_gap = max(g$rel_gap), median_rel_gap = median(g$rel_gap), seconds = g$seconds[1])))
summ <- summ[order(summ$regime, summ$model, summ$max_rel_gap), ]
rownames(summ) <- NULL
print(summ, digits = 3)
write.csv(summ, "benchmarks/results/quality-summary.csv", row.names = FALSE)
