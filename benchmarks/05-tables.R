## Markdown tables for the benchmark write-up (written to benchmarks/results/*.md).
##   Rscript benchmarks/05-tables.R
source("benchmarks/setup.R")
out <- "benchmarks/results/tables.md"
cat("", file = out)

## ---------------- speed: one wide table per model ----------------
sp <- read.csv("benchmarks/results/speed.csv", stringsAsFactors = FALSE)
sp$size <- paste0(sp$n, " x ", sp$p)
for (m in unique(sp$model)) {
  d <- sp[sp$model == m, ]
  sizes <- unique(d$size)
  meths <- unique(d$method)
  tab <- data.frame(method = meths, stringsAsFactors = FALSE)
  for (sz in sizes) {
    v <- d$seconds[match(paste(meths, sz), paste(d$method, d$size))]
    ref <- d$seconds[d$size == sz & d$package == "hdstats" & !grepl("is_exact", d$method)][1]
    tab[[sz]] <- ifelse(is.na(v), "", sprintf("%.2f s (%.1fx)", v, v / ref))
  }
  names(tab)[1] <- "method (seconds; ratio to hdstats default)"
  write_md_table(tab, out, append = TRUE, caption = paste0("**", m, "**: median elapsed seconds of ", 3, " runs; the ratio is time / time of hdstats' default solver at that size."))
}

## ---------------- optimization quality ----------------
q <- read.csv("benchmarks/results/quality.csv", stringsAsFactors = FALSE)
for (rg in unique(q$regime)) for (m in unique(q$model)) {
  d <- q[q$regime == rg & q$model == m, ]
  lams <- sort(unique(d$lambda), decreasing = TRUE)
  solvers <- unique(d$solver)
  tab <- data.frame(solver = solvers, stringsAsFactors = FALSE)
  for (l in lams) {
    v <- d$rel_gap[match(paste(solvers, l), paste(d$solver, d$lambda))]
    dfv <- d$df[match(paste(solvers, l), paste(d$solver, d$lambda))]
    tab[[sprintf("lambda = %.3g", l)]] <- ifelse(is.na(v), "", sprintf("%.2g (df %d)", v, dfv))
  }
  tab$seconds <- sprintf("%.2f", d$seconds[match(solvers, d$solver)])
  write_md_table(tab, out, append = TRUE,
                 caption = paste0("**", m, ", ", rg, "**: relative objective gap (objective - best) / best at five lambda values; df = number of nonzero coefficients; seconds = time for the five fits."))
}

## ---------------- accuracy ----------------
a <- do.call(rbind, lapply(list.files("benchmarks/results", "^accuracy-raw", full.names = TRUE), read.csv, stringsAsFactors = FALSE))
se <- function(v) sd(v) / sqrt(length(v))
fmt <- function(m, s) sprintf("%.3f (%.3f)", m, s)
for (st in unique(a$setting)) for (e in unique(a$error[a$setting == st])) {
  d <- a[a$setting == st & a$error == e, ]
  g <- split(d, d$method)
  tab <- do.call(rbind, lapply(g, function(x) data.frame(
    method = x$method[1],
    `L2 error` = if (st == "regression") fmt(mean(x$l2), se(x$l2)) else "",
    `test error` = fmt(mean(x$test_mae), se(x$test_mae)),
    TPR = fmt(mean(x$tpr), se(x$tpr)), FPR = fmt(mean(x$fpr), se(x$fpr)),
    df = sprintf("%.1f", mean(x$df)), `CV seconds` = sprintf("%.1f", mean(x$seconds)),
    check.names = FALSE, stringsAsFactors = FALSE)))
  if (st == "classification") tab$`L2 error` <- NULL
  tab <- tab[order(if (st == "regression") as.numeric(sub(" .*", "", tab$`L2 error`)) else as.numeric(sub(" .*", "", tab$`test error`))), ]
  rownames(tab) <- NULL
  write_md_table(tab, out, append = TRUE,
                 caption = paste0("**", st, ", ", e, "** (n = 200, p = 500, s = 10, ", length(unique(d$rep)), " replications; mean (s.e.); test error = ",
                                  if (st == "regression") "mean absolute error on 1000 test observations" else "misclassification rate on 1000 test observations", ")."))
}

## ---------------- platform ----------------
pi <- platform_info()
pv <- pkg_versions(c("hdstats", "hqreg", "sparseSVM", "conquer", "quantreg", "gcdnet", "LiblineaR", "Rfit", "glmnet"))
cat("**Platform**: ", paste(names(pi), pi, sep = ": ", collapse = "; "), "\n\n", file = out, append = TRUE)
write_md_table(pv, out, append = TRUE, caption = "**Package versions**")
cat("tables written to", out, "\n")
