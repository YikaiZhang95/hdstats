## ---------------------------------------------------------------------------
## ACCURACY: does the speed cost anything statistically?
##
## Identical protocol for every package: fit the full path on a training set,
## pick lambda by minimizing error on an independent validation set, then score
## the selected model on an independent test set. Averaged over 10 replications.
##
## Reported: test error (misclassification for SVM; MAE for QR/Huber),
##           number of selected variables, and slope estimation error
##           ||beta_hat - beta*||_2 (not identified for SVM, shown as NA).
## ---------------------------------------------------------------------------
suppressMessages({library(hdstats); library(hqreg); library(sparseSVM); library(conquer)})

n <- 200; p <- 1000; s <- 10; REPS <- 10
b_star <- c(rep(1.5, s), rep(0, p - s))

## full coefficient path (intercept in row 1) from each package
path_hdsvm <- function(x, y) { f <- hdsvm(x, y, lam2 = 0, eps = 1e-6); rbind(f$b0, as.matrix(f$beta)) }
path_ssvm  <- function(x, y) as.matrix(sparseSVM(x, y, eps = 1e-6)$weights)
path_hdqr  <- function(x, y) { f <- hdqr(x, y, tau = .5, lam2 = 0, eps = 1e-6); rbind(f$b0, as.matrix(f$beta)) }
path_hqrq  <- function(x, y) as.matrix(hqreg(x, y, method = "quantile", tau = .5, eps = 1e-6)$beta)
path_conq  <- function(x, y) conquer.reg(x, y, tau = .5, penalty = "lasso",
                                          lambda = hdqr(x, y, tau = .5, lam2 = 0, eps = 1e-6)$lambda)$coeff
path_hdhub <- function(x, y) { f <- hdhuber(x, y, delta = 1, lam2 = 0, eps = 1e-6); rbind(f$b0, as.matrix(f$beta)) }
path_hqrh  <- function(x, y) as.matrix(hqreg(x, y, method = "huber", eps = 1e-6)$beta)

mae    <- function(pred, y) mean(abs(pred - y))
mcr    <- function(pred, y) mean(sign(pred) != y)
esterr <- function(b, bt) sqrt(sum((b - bt)^2))

evaluate <- function(B, xva, xte, yva, yte, errfun, btrue = NULL) {
  k  <- which.min(apply(cbind(1, xva) %*% B, 2, errfun, y = yva))
  bk <- B[, k]
  c(test = errfun(cbind(1, xte) %*% bk, yte),
    nsel = sum(bk[-1] != 0),
    est  = if (is.null(btrue)) NA else esterr(bk[-1], btrue))
}

acc <- function(model) {
  M <- matrix(NA, REPS, 3, dimnames = list(NULL, c("test", "nsel", "est")))
  for (r in 1:REPS) {
    set.seed(100 + r)
    mk <- function() matrix(rnorm(n * p), n, p)
    xtr <- mk(); xva <- mk(); xte <- mk()
    if (model %in% c("hdsvm", "sparseSVM")) {
      g <- function(X) ifelse(drop(X %*% b_star) + rnorm(n) > 0, 1, -1)
      ef <- mcr; bt <- NULL
    } else {
      g <- function(X) drop(X %*% b_star + rt(n, df = 3))
      ef <- mae; bt <- b_star
    }
    ytr <- g(xtr); yva <- g(xva); yte <- g(xte)
    B <- switch(model,
      hdsvm = path_hdsvm(xtr, ytr), sparseSVM = path_ssvm(xtr, ytr),
      hdqr = path_hdqr(xtr, ytr),  hqreg_qr  = path_hqrq(xtr, ytr),
      conquer = path_conq(xtr, ytr),
      hdhuber = path_hdhub(xtr, ytr), hqreg_hub = path_hqrh(xtr, ytr))
    M[r, ] <- evaluate(B, xva, xte, yva, yte, ef, bt)
  }
  c(test_mean = mean(M[, "test"]), test_sd = sd(M[, "test"]),
    nsel_mean = mean(M[, "nsel"]), est_mean = mean(M[, "est"]))
}

models <- c("hdsvm", "sparseSVM", "hdqr", "hqreg_qr", "conquer", "hdhuber", "hqreg_hub")
out <- as.data.frame(round(t(sapply(models, acc)), 4))
out$model <- rownames(out)
out <- out[, c("model", "test_mean", "test_sd", "nsel_mean", "est_mean")]
dir.create("results", showWarnings = FALSE)
write.csv(out, "results/bench_accuracy.csv", row.names = FALSE)
print(out, row.names = FALSE)
cat("\ntest = misclassification (SVM) or test MAE (QR/Huber);",
    "est = ||b_hat - b*||_2 on slopes.  True model: s =", s, "of p =", p, "\n")
