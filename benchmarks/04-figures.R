## Figures for the benchmark write-up (PNG, base graphics).
##   Rscript benchmarks/04-figures.R
## Colors follow the package (fixed slots of a colorblind-validated palette),
## never the rank of a series; hdstats is always slot 1.
pal <- c(hdstats = "#2a78d6", hqreg = "#eb6834", conquer = "#1baf7a", quantreg = "#eda100",
         sparseSVM = "#e87ba4", gcdnet = "#008300", LiblineaR = "#4a3aa7", Rfit = "#e34948",
         glmnet = "#4a3aa7")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid <- "#e6e5e1"; surface <- "#fcfcfb"
dir.create("benchmarks/figures", showWarnings = FALSE)

## ---------------- speed ----------------
sp <- read.csv("benchmarks/results/speed.csv", stringsAsFactors = FALSE)
sp$size <- paste0(sp$n, "x", sp$p)
sp$series <- ifelse(grepl("is_exact", sp$method), paste0(sp$package, " (is_exact)"), sp$package)
panels <- c("Huber", "Quantile", "SVM")
sizes <- unique(sp$size[sp$model %in% panels])          # in the order they were run
png("benchmarks/figures/speed.png", width = 2100, height = 900, res = 160, bg = surface)
par(mfrow = c(1, 3), mar = c(4.5, 9, 3, 1), family = "sans", col.axis = ink2, col.lab = ink2, fg = ink2)
for (m in panels) {
  d <- sp[sp$model == m & !grepl("10 lambda|single", sp$note), ]
  d$yi <- match(d$size, sizes)
  series <- unique(d$series)
  off <- seq(-0.28, 0.28, length.out = length(series))
  plot(NA, xlim = range(d$seconds) * c(0.6, 2), ylim = c(0.5, length(sizes) + 1.6), log = "x", axes = FALSE,
       xlab = "seconds for a 100-value path (log scale)", ylab = "", main = m, col.main = ink)
  abline(v = axTicks(1), col = grid); abline(h = seq_len(length(sizes) - 1) + 0.5, col = grid)
  axis(1, at = axTicks(1), labels = format(axTicks(1), drop0trailing = TRUE, scientific = FALSE), tick = FALSE)
  axis(2, at = seq_along(sizes), labels = paste0("n x p = ", sizes), las = 1, tick = FALSE, cex.axis = 0.9)
  for (k in seq_along(series)) {
    dd <- d[d$series == series[k], ]
    col <- pal[[sub(" .*", "", series[k])]]; pch <- if (grepl("is_exact", series[k])) 24 else 21
    points(dd$seconds, dd$yi + off[k], pch = pch, bg = col, col = surface, cex = 1.7, lwd = 1.5)
    if (length(dd$yi) > 1) segments(dd$seconds, dd$yi + off[k], dd$seconds, dd$yi + off[k])
  }
  legend("top", legend = series, pch = ifelse(grepl("is_exact", series), 24, 21),
         pt.bg = pal[sub(" .*", "", series)], col = surface, pt.cex = 1.4, bty = "n", ncol = 2,
         cex = 0.85, text.col = ink2)
}
dev.off()

## ---------------- accuracy ----------------
ac <- do.call(rbind, lapply(list.files("benchmarks/results", "^accuracy-raw", full.names = TRUE), read.csv, stringsAsFactors = FALSE))
dotplot <- function(d, value, xlab, main) {
  agg <- aggregate(d[[value]], list(method = d$method), function(v) c(mean = mean(v), se = if (length(v) > 1) sd(v) / sqrt(length(v)) else 0))
  agg <- data.frame(method = agg$method, mean = agg$x[, "mean"], se = agg$x[, "se"])
  agg <- agg[order(-agg$mean), ]
  y <- seq_len(nrow(agg))
  is_hd <- grepl("hdstats", agg$method)
  par(mar = c(4.5, 15, 3, 1))
  plot(NA, xlim = range(c(agg$mean - agg$se, agg$mean + agg$se)) * c(0.9, 1.05), ylim = c(0.5, nrow(agg) + 0.5),
       axes = FALSE, xlab = xlab, ylab = "", main = main, col.main = ink)
  abline(v = axTicks(1), col = grid); axis(1, tick = FALSE); axis(2, at = y, labels = agg$method, las = 1, tick = FALSE, cex.axis = 0.8)
  segments(agg$mean - agg$se, y, agg$mean + agg$se, y, col = ifelse(is_hd, pal[["hdstats"]], ink2), lwd = 2)
  points(agg$mean, y, pch = 21, bg = ifelse(is_hd, pal[["hdstats"]], ink2), col = surface, cex = 1.5, lwd = 1.5)
}
reg <- ac[ac$setting == "regression", ]
errs <- unique(reg$error)
png("benchmarks/figures/accuracy-regression.png", width = 2400, height = 800, res = 160, bg = surface)
par(mfrow = c(1, length(errs)), family = "sans", col.axis = ink2, col.lab = ink2, fg = ink2)
for (e in errs) dotplot(reg[reg$error == e, ], "l2", "L2 estimation error (mean +/- s.e.)", paste0(e, " errors"))
dev.off()
cl <- ac[ac$setting == "classification", ]
if (nrow(cl)) {
  png("benchmarks/figures/accuracy-classification.png", width = 2000, height = 700, res = 160, bg = surface)
  par(mfrow = c(1, 2), family = "sans", col.axis = ink2, col.lab = ink2, fg = ink2)
  dotplot(cl, "test_mae", "test misclassification rate (mean +/- s.e.)", "Test error")
  dotplot(cl, "fpr", "false positive rate (mean +/- s.e.)", "False positive rate")
  dev.off()
}
cat("figures written to benchmarks/figures/\n")
