## Statistical performance: estimation error, variable selection and test
## error at the cross-validation-selected lambda of each package, over
## replicated simulations.
##   Rscript benchmarks/03-accuracy.R
## Environment variables: HDSTATS_BENCH_REPS (default 10), HDSTATS_BENCH_QUICK.
source("benchmarks/setup.R")
suppressPackageStartupMessages({ library(hqreg); library(sparseSVM); library(gcdnet) })

REPS <- as.integer(Sys.getenv("HDSTATS_BENCH_REPS", "10"))
QUICK <- nzchar(Sys.getenv("HDSTATS_BENCH_QUICK"))
## HDSTATS_BENCH_PART: "all", "regression", "extra" (lam2 = 0 variants),
## "extra2" (rank regression with lam2 = 0) or "classification"; the parts
## write separate files that 04-figures.R and 05-tables.R merge.
PART <- Sys.getenv("HDSTATS_BENCH_PART", "all")
do_part <- function(x) PART == "all" || PART == x
n <- 200; p <- if (QUICK) 100 else 500; s <- 10; ntest <- 1000
NF <- 5
errors <- c("normal", "t3", "contaminated")

metrics <- function(b0, beta, d) {
  beta <- as.numeric(beta); sel <- which(abs(beta) > 1e-8); S <- seq_len(d$s)
  c(l2 = sqrt(sum((beta - d$beta)^2)),
    tpr = length(intersect(sel, S)) / length(S),
    fpr = length(setdiff(sel, S)) / (length(beta) - length(S)),
    df = length(sel),
    test_mae = mean(abs(d$ytest - b0 - drop(d$xtest %*% beta))))
}

res <- list()
add <- function(setting, error, method, rep, m, secs)
  res[[length(res) + 1]] <<- data.frame(setting = setting, error = error, method = method, rep = rep,
                                        t(m), seconds = secs, stringsAsFactors = FALSE)

## ---------------- regression ----------------
for (err in errors) for (r in seq_len(REPS)) {
  if (!(do_part("regression") || do_part("extra") || do_part("extra2"))) break
  d <- sim_reg(n, p, s = s, error = err, seed = 100 * r, ntest = ntest)
  foldid <- sample(rep(seq_len(NF), length.out = n))
  run <- function(method, expr) {
    secs <- system.time(cf <- expr)[["elapsed"]]
    add("regression", err, method, r, metrics(cf[1], cf[-1], d), secs)
  }
  if (do_part("regression")) {
  run("glmnet (least squares lasso)", {
    cv <- cv.glmnet(d$x, d$y, foldid = foldid); as.numeric(stats::coef(cv, s = "lambda.min")) })
  run("hdstats::hdhuber (delta = 1)", {
    cv <- cv.hdhuber(d$x, d$y, delta = 1, foldid = foldid); as.numeric(stats::coef(cv, s = "lambda.min")) })
  run("hqreg (huber, gamma = 1)", {
    cv <- cv.hqreg(d$x, d$y, method = "huber", gamma = 1, nfolds = NF, fold.id = foldid)
    as.numeric(stats::coef(cv, lambda = cv$lambda.min)) })
  run("hdstats::hdqr (tau = 0.5)", {
    cv <- cv.hdqr(d$x, d$y, tau = 0.5, foldid = foldid); as.numeric(stats::coef(cv, s = "lambda.min")) })
  run("hqreg (quantile, tau = 0.5)", {
    cv <- cv.hqreg(d$x, d$y, method = "quantile", tau = 0.5, nfolds = NF, fold.id = foldid)
    as.numeric(stats::coef(cv, lambda = cv$lambda.min)) })
  if (has("conquer")) run("conquer (lasso, tau = 0.5)", {
    cv <- conquer::conquer.cv.reg(d$x, d$y, tau = 0.5, kfolds = NF); as.numeric(cv$coeff.min) })
  run("hdstats::hdrr (rank)", {
    cv <- cv.hdrr(d$x, d$y, foldid = foldid, nlambda = 30, lambda.factor = 0.05)
    as.numeric(stats::coef(cv, s = "lambda.min")) })
  }
  ## lasso-only variants: hdqr() and hdrr() default to lam2 = 0.01
  if (do_part("extra")) run("hdstats::hdqr (tau = 0.5, lam2 = 0)", {
    cv <- cv.hdqr(d$x, d$y, tau = 0.5, lam2 = 0, foldid = foldid); as.numeric(stats::coef(cv, s = "lambda.min")) })
  if (do_part("extra2")) run("hdstats::hdrr (rank, lam2 = 0)", {
    cv <- cv.hdrr(d$x, d$y, lam2 = 0, foldid = foldid, nlambda = 30, lambda.factor = 0.05)
    as.numeric(stats::coef(cv, s = "lambda.min")) })
  cat(sprintf("regression  error=%-12s rep %d done\n", err, r))
}

## ---------------- classification ----------------
## Class labels are obtained with each package's own predict() method so that
## label-coding conventions (sparseSVM codes the first level as +1) do not
## matter; the selection metrics only depend on which coefficients are nonzero.
for (r in seq_len(REPS)) {
  if (!do_part("classification")) break
  dc <- sim_class(n, p, s = s, seed = 200 * r, ntest = ntest)
  foldid <- sample(rep(seq_len(NF), length.out = n))
  cmetrics <- function(beta, pred) {
    beta <- as.numeric(beta); sel <- which(abs(beta) > 1e-8); S <- seq_len(dc$s)
    c(l2 = NA, tpr = length(intersect(sel, S)) / length(S), fpr = length(setdiff(sel, S)) / (p - length(S)),
      df = length(sel), test_mae = mean(as.numeric(as.character(pred)) != dc$ytest))   # misclassification rate
  }
  run <- function(method, expr) {
    secs <- system.time(out <- expr)[["elapsed"]]
    add("classification", "flip 5%", method, r, cmetrics(out$beta, out$pred), secs)
  }
  run("glmnet (logistic lasso)", {
    cv <- cv.glmnet(dc$x, dc$y, family = "binomial", foldid = foldid)
    list(beta = stats::coef(cv, s = "lambda.min")[-1], pred = stats::predict(cv, dc$xtest, s = "lambda.min", type = "class")) })
  run("hdstats::hdsvm", {
    cv <- cv.hdsvm(dc$x, dc$y, foldid = foldid)
    list(beta = stats::coef(cv, s = "lambda.min")[-1], pred = stats::predict(cv, dc$xtest, s = "lambda.min")) })
  run("hdstats::hdsvm (is_exact = TRUE)", {
    cv <- cv.hdsvm(dc$x, dc$y, foldid = foldid, is_exact = TRUE)
    list(beta = stats::coef(cv, s = "lambda.min")[-1], pred = stats::predict(cv, dc$xtest, s = "lambda.min")) })
  run("sparseSVM", {
    cv <- cv.sparseSVM(dc$x, dc$y, nfolds = NF, fold.id = foldid)
    list(beta = stats::coef(cv$fit, lambda = cv$lambda.min)[-1],
         pred = stats::predict(cv$fit, dc$xtest, lambda = cv$lambda.min, type = "class")) })
  run("gcdnet (Huberized squared hinge)", {
    cv <- cv.gcdnet(dc$x, dc$y, method = "hhsvm", foldid = foldid)
    list(beta = gcdnet::coef(cv, s = "lambda.min")[-1], pred = gcdnet::predict(cv, dc$xtest, s = "lambda.min", type = "class")) })
  cat(sprintf("classification rep %d done\n", r))
}

res <- do.call(rbind, res)
write.csv(res, sprintf("benchmarks/results/accuracy-raw-%s.csv", PART), row.names = FALSE)
agg <- aggregate(cbind(l2, tpr, fpr, df, test_mae, seconds) ~ setting + error + method, data = res, FUN = mean, na.action = na.pass)
print(agg[order(agg$setting, agg$error, agg$method), ], digits = 3)
