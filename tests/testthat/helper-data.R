## Small simulated data sets and reference objective functions used by the
## tests. Nothing here is exported.

sim_reg <- function(n = 60, p = 30, seed = 1) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(2, -1.5, 1, rep(0, p - 3))
  y <- drop(x %*% beta) + rnorm(n)
  list(x = x, y = y, beta = beta)
}

sim_class <- function(n = 80, p = 30, seed = 2) {
  set.seed(seed)
  x <- matrix(rnorm(n * p), n, p)
  beta <- c(1.5, -1.5, 1, rep(0, p - 3))
  prob <- plogis(drop(x %*% beta))
  y <- 2 * rbinom(n, 1, prob) - 1
  list(x = x, y = y)
}

## Objective functions on the original scale; they match the models fitted
## with standardize = FALSE (the intercept is never penalized).
huber_loss_ref <- function(r, delta)
  ifelse(abs(r) <= delta, 0.5 * r^2, delta * (abs(r) - 0.5 * delta))

check_loss_ref <- function(r, tau) r * (tau - (r < 0))

hinge_loss_ref <- function(m) pmax(1 - m, 0)

penalty_ref <- function(beta, lambda, lam2, pf = 1, pf2 = 1)
  lambda * sum(pf * abs(beta)) + 0.5 * lam2 * sum(pf2 * beta^2)

huber_obj <- function(b0, beta, x, y, delta, lambda, lam2, pf = 1, pf2 = 1) {
  r <- y - b0 - drop(x %*% beta)
  mean(huber_loss_ref(r, delta)) + penalty_ref(beta, lambda, lam2, pf, pf2)
}

qr_obj <- function(b0, beta, x, y, tau, lambda, lam2, pf = 1, pf2 = 1) {
  r <- y - b0 - drop(x %*% beta)
  mean(check_loss_ref(r, tau)) + penalty_ref(beta, lambda, lam2, pf, pf2)
}

svm_obj <- function(b0, beta, x, y, lambda, lam2, pf = 1, pf2 = 1) {
  m <- y * (b0 + drop(x %*% beta))
  mean(hinge_loss_ref(m)) + penalty_ref(beta, lambda, lam2, pf, pf2)
}

## coefficient vector (without intercept) and intercept at path index l
path_coef <- function(fit, l) {
  list(b0 = unname(fit$b0[l]), beta = as.vector(fit$beta[, l]))
}
