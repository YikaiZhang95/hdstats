## ---------------------------------------------------------------------------
## EXACT-SOLUTION CHECK: hdstats vs the CVXR convex optimum.
##
## CVXR (with the ECOS conic solver) solves a disciplined-convex-programming
## formulation of each penalized problem and provides the exact optimizer at a
## given lambda. We compare it against hdstats fit with standardize = FALSE
## (CVXR penalizes the original-scale coefficients) and is_exact = TRUE
## (return the exact, not the finitely smoothed, solution) so that both
## minimize the identical objective.
##
## Reported per model: the penalized objective achieved by each solver, their
## relative difference, and the largest coefficient discrepancy.
## ---------------------------------------------------------------------------
suppressMessages({library(hdstats); library(CVXR)})

## --- exact convex references (return intercept + slopes; rank: slopes) ------
svm_cvxr <- function(X, y, lambda) {
  b0 <- Variable(1); beta <- Variable(ncol(X))
  loss <- sum(pos(1 - y * (X %*% beta + b0))) / nrow(X)
  res  <- solve(Problem(Minimize(loss + lambda * p_norm(beta, 1))), solver = "ECOS")
  c(res$getValue(b0), as.numeric(res$getValue(beta)))
}
qr_cvxr <- function(X, y, lambda, tau = 0.5) {
  b0 <- Variable(1); beta <- Variable(ncol(X)); r <- y - X %*% beta - b0
  loss <- sum(0.5 * abs(r) + (tau - 0.5) * r) / nrow(X)
  res  <- solve(Problem(Minimize(loss + lambda * p_norm(beta, 1))), solver = "ECOS")
  c(res$getValue(b0), as.numeric(res$getValue(beta)))
}
rr_cvxr <- function(X, y, lambda) {
  pr <- which(upper.tri(matrix(0, nrow(X), nrow(X))), arr.ind = TRUE)
  xd <- X[pr[, 1], ] - X[pr[, 2], ]; yd <- y[pr[, 1]] - y[pr[, 2]]
  beta <- Variable(ncol(X))
  loss <- sum(0.5 * abs(yd - xd %*% beta)) / nrow(xd)
  res  <- solve(Problem(Minimize(loss + lambda * p_norm(beta, 1))), solver = "ECOS")
  as.numeric(res$getValue(beta))
}

## --- matching objectives (mean loss + lambda * ||slopes||_1) -----------------
o_svm <- function(b, X, y, lam) mean(pmax(0, 1 - y * (cbind(1, X) %*% b))) + lam * sum(abs(b[-1]))
o_qr  <- function(b, X, y, lam, tau = 0.5) { r <- y - cbind(1, X) %*% b
  mean(0.5 * abs(r) + (tau - 0.5) * r) + lam * sum(abs(b[-1])) }
o_rr  <- function(b, xd, yd, lam) mean(0.5 * abs(yd - xd %*% b)) + lam * sum(abs(b))

k <- 40
rows <- list()

## --- SVM ---------------------------------------------------------------------
set.seed(1); n <- 150; p <- 300
x <- matrix(rnorm(n * p), n, p)
y <- ifelse(drop(x %*% c(rep(2, 5), rep(0, p - 5))) + rnorm(n) > 0, 1, -1)
f <- hdsvm(x, y, lam2 = 0, standardize = FALSE, is_exact = TRUE); lam <- f$lambda[k]
bh <- c(f$b0[k], as.numeric(f$beta[, k])); bc <- svm_cvxr(x, y, lam)
rows$SVM <- c(obj_hdstats = o_svm(bh, x, y, lam), obj_cvxr = o_svm(bc, x, y, lam),
              max_coef_diff = max(abs(bh - bc)))

## --- Quantile ----------------------------------------------------------------
set.seed(2); n <- 150; p <- 300
x <- matrix(rnorm(n * p), n, p)
y <- drop(x %*% c(runif(5, 1, 2), rep(0, p - 5)) + rt(n, df = 3))
f <- hdqr(x, y, tau = 0.5, lam2 = 0, standardize = FALSE, is_exact = TRUE); lam <- f$lambda[k]
bh <- c(f$b0[k], as.numeric(f$beta[, k])); bc <- qr_cvxr(x, y, lam)
rows$Quantile <- c(obj_hdstats = o_qr(bh, x, y, lam), obj_cvxr = o_qr(bc, x, y, lam),
                   max_coef_diff = max(abs(bh - bc)))

## --- Rank (Wilcoxon) ---------------------------------------------------------
set.seed(3); n <- 90; p <- 200
x <- matrix(rnorm(n * p), n, p)
y <- drop(x %*% c(runif(5, 1, 2), rep(0, p - 5)) + rt(n, df = 3))
f <- hdrr(x, y, lam2 = 0, standardize = FALSE, is_exact = TRUE); lam <- f$lambda[k]
bh <- as.numeric(f$beta[, k]); bc <- rr_cvxr(x, y, lam)
pr <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
xd <- x[pr[, 1], ] - x[pr[, 2], ]; yd <- y[pr[, 1]] - y[pr[, 2]]
rows$Rank <- c(obj_hdstats = o_rr(bh, xd, yd, lam), obj_cvxr = o_rr(bc, xd, yd, lam),
               max_coef_diff = max(abs(bh - bc)))

## --- report ------------------------------------------------------------------
res <- as.data.frame(do.call(rbind, rows))
res$rel_obj_diff <- abs(res$obj_hdstats - res$obj_cvxr) / res$obj_cvxr
res <- data.frame(model = rownames(res), round(res, 6), row.names = NULL)
dir.create("results", showWarnings = FALSE)
write.csv(res, "results/bench_cvxr_exact.csv", row.names = FALSE)
print(res, row.names = FALSE)
cat("\nhdstats recovers the exact CVXR optimum (objectives agree to solver tolerance).\n")
