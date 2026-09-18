## Common helpers for the hdstats benchmarks.
##
## The competitor packages are expected to be installed; the environment
## variable HDSTATS_BENCH_LIB may point to an additional library path.

lib <- Sys.getenv("HDSTATS_BENCH_LIB", unset = NA)
if (!is.na(lib) && nzchar(lib)) .libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(hdstats)
  library(glmnet)
  library(quantreg)
})
has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

## ---- data generation ---------------------------------------------------

## Correlated (AR(1), rho = 0.5) Gaussian design, sparse coefficient vector
## with s nonzero entries, and errors from the requested family.
sim_reg <- function(n, p, s = 10, error = c("normal", "t3", "t2", "contaminated"),
                    rho = 0.5, seed = 1, ntest = 0) {
  error <- match.arg(error)
  set.seed(seed)
  gen_x <- function(m) {
    z <- matrix(rnorm(m * p), m, p)
    if (rho > 0) {           # AR(1) correlation via a one-pass filter
      for (j in 2:p) z[, j] <- rho * z[, j - 1] + sqrt(1 - rho^2) * z[, j]
    }
    z
  }
  beta <- numeric(p)
  beta[seq_len(s)] <- rep(c(2, -1.5, 1, -0.8, 0.6), length.out = s)
  gen_e <- function(m) switch(error,
    normal = rnorm(m),
    t3 = rt(m, df = 3),
    t2 = rt(m, df = 2),
    contaminated = ifelse(runif(m) < 0.1, rnorm(m, sd = 10), rnorm(m)))
  x <- gen_x(n); y <- drop(x %*% beta) + gen_e(n)
  out <- list(x = x, y = y, beta = beta, s = s, error = error)
  if (ntest > 0) {
    out$xtest <- gen_x(ntest)
    out$ytest <- drop(out$xtest %*% beta) + gen_e(ntest)
  }
  out
}

## Binary classification with a sparse linear signal (two blocks of equal
## coefficients with opposite signs); Gaussian noise before taking the sign
## and labels flipped with probability `flip` make the problem
## non-separable.
sim_class <- function(n, p, s = 10, rho = 0.5, flip = 0.05, sd_noise = 1, seed = 1, ntest = 0) {
  set.seed(seed)
  gen_x <- function(m) {
    z <- matrix(rnorm(m * p), m, p)
    if (rho > 0) for (j in 2:p) z[, j] <- rho * z[, j - 1] + sqrt(1 - rho^2) * z[, j]
    z
  }
  beta <- numeric(p)
  beta[seq_len(s)] <- rep(c(1, -1), each = ceiling(s / 2))[seq_len(s)]
  gen_y <- function(x) {
    y <- sign(drop(x %*% beta) + rnorm(nrow(x), sd = sd_noise))
    y[y == 0] <- 1
    fl <- runif(nrow(x)) < flip
    y[fl] <- -y[fl]
    y
  }
  x <- gen_x(n); y <- gen_y(x)
  out <- list(x = x, y = y, beta = beta, s = s)
  if (ntest > 0) { out$xtest <- gen_x(ntest); out$ytest <- gen_y(out$xtest) }
  out
}

## ---- objective functions (raw scale, intercept unpenalized) ------------

huber_loss <- function(r, delta) ifelse(abs(r) <= delta, r^2 / 2, delta * (abs(r) - delta / 2))
check_loss <- function(r, tau) r * (tau - (r < 0))
hinge_loss <- function(m) pmax(1 - m, 0)

obj_huber <- function(b0, beta, x, y, delta, lambda, lam2 = 0)
  mean(huber_loss(y - b0 - drop(x %*% beta), delta)) + lambda * sum(abs(beta)) + lam2 / 2 * sum(beta^2)
obj_qr <- function(b0, beta, x, y, tau, lambda, lam2 = 0)
  mean(check_loss(y - b0 - drop(x %*% beta), tau)) + lambda * sum(abs(beta)) + lam2 / 2 * sum(beta^2)
obj_svm <- function(b0, beta, x, y, lambda, lam2 = 0)
  mean(hinge_loss(y * (b0 + drop(x %*% beta)))) + lambda * sum(abs(beta)) + lam2 / 2 * sum(beta^2)

## ---- timing ------------------------------------------------------------

## median elapsed time of `reps` evaluations of `expr` (evaluated in the
## caller's environment); the first evaluation is a warm-up when reps > 1
time_it <- function(expr, reps = 3) {
  e <- substitute(expr); env <- parent.frame()
  t <- vapply(seq_len(reps), function(i) system.time(eval(e, env))[["elapsed"]], numeric(1))
  median(t)
}

## ---- reporting ---------------------------------------------------------

platform_info <- function() {
  c(R = R.version.string,
    platform = R.version$platform,
    os = paste(Sys.info()[c("sysname", "release")], collapse = " "),
    cpu = tryCatch(trimws(sub("model name\\s*:", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1])), error = function(e) NA),
    BLAS = basename(extSoftVersion()[["BLAS"]]),
    compiler = tryCatch(trimws(system2("g++", "--version", stdout = TRUE)[1]), error = function(e) NA))
}

pkg_versions <- function(pkgs) {
  v <- vapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_), "")
  data.frame(package = pkgs, version = v, row.names = NULL)
}

write_md_table <- function(df, file, digits = 3, append = FALSE, caption = NULL) {
  fmt <- function(v) if (is.numeric(v)) formatC(v, digits = digits, format = "g") else as.character(v)
  cells <- as.data.frame(lapply(df, fmt), stringsAsFactors = FALSE)
  lines <- c(if (!is.null(caption)) c(caption, ""),
             paste0("| ", paste(names(df), collapse = " | "), " |"),
             paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
             apply(cells, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
  cat(lines, file = file, sep = "\n", append = append)
}
