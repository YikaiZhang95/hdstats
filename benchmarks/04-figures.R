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
png("benchmarks/figures/speed.png", width = 2100, height = 800, res = 160, bg = surface)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1), family = "sans", col.axis = ink2, col.lab = ink2, fg = ink2)
for (m in panels) {
  d <- sp[sp$model == m & !grepl("10 lambda|single", sp$note), ]
  sizes <- unique(sp$size[sp$model == m]); sizes <- sizes[order(as.numeric(sub("x.*", "", sizes)) * as.numeric(sub(".*x", "", sizes)))]
  d$xi <- match(d$size, sizes)
  plot(NA, xlim = c(0.8, length(sizes) + 0.9), ylim = range(d$seconds) * c(0.7, 1.6), log = "y", axes = FALSE,
       xlab = "problem size (n x p)", ylab = "seconds for a 100-value path (log scale)", main = m, col.main = ink)
  abline(h = axTicks(2), col = grid, lwd = 1); axis(1, at = seq_along(sizes), labels = sizes, tick = FALSE); axis(2, las = 1, tick = FALSE)
  for (s in unique(d$series)) {
    dd <- d[d$series == s, ]; dd <- dd[order(dd$xi), ]
    col <- pal[[sub(" .*", "", s)]]; lty <- if (grepl("is_exact", s)) 2 else 1
    lines(dd$xi, dd$seconds, col = col, lwd = 2, lty = lty)
    points(dd$xi, dd$seconds, col = col, bg = surface, pch = 21, cex = 1.3, lwd = 2)
    text(max(dd$xi) + 0.08, dd$seconds[which.max(dd$xi)], s, adj = c(0, 0.5), cex = 0.75, col = ink2, xpd = NA)
  }
}
dev.off()

## ---------------- accuracy ----------------
ac <- read.csv("benchmarks/results/accuracy-raw.csv", stringsAsFactors = FALSE)
dotplot <- function(d, value, xlab, main) {
  agg <- aggregate(d[[value]], list(method = d$method), function(v) c(mean = mean(v), se = sd(v) / sqrt(length(v))))
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
png("benchmarks/figures/accuracy-classification.png", width = 1400, height = 700, res = 160, bg = surface)
par(mfrow = c(1, 2), family = "sans", col.axis = ink2, col.lab = ink2, fg = ink2)
cl <- ac[ac$setting == "classification", ]
dotplot(cl, "test_mae", "test misclassification rate (mean +/- s.e.)", "Classification: test error")
dotplot(cl, "fpr", "false positive rate (mean +/- s.e.)", "Classification: false positives")
dev.off()
cat("figures written to benchmarks/figures/\n")
