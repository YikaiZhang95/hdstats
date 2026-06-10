## ---------------------------------------------------------------------------
## Figure: 100-lambda path runtime, hdstats vs peers, across (n, p) settings.
## Reads results/bench_speed.csv (produced by 02_speed.R).
## ---------------------------------------------------------------------------
sp <- read.csv("results/bench_speed.csv")
sp <- sp[order(sp$n, sp$p), ]
sp$lab <- sprintf("n=%d\np=%d", sp$n, sp$p)
blue <- "#2c7fb8"; orange <- "#fd8d3c"; grey <- "#bdbdbd"

png("results/hdstats_benchmark.png", width = 1500, height = 540, res = 110)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3.2, 1), cex.axis = 0.9)

bars <- function(mat, cols, leg, main) {
  barplot(t(mat), beside = TRUE, col = cols, names.arg = sp$lab, border = NA,
          main = main, ylab = "path time (s)", ylim = c(0, max(mat) * 1.18),
          cex.names = 0.85)
  legend("topleft", legend = leg, fill = cols, bty = "n", border = NA, cex = 0.95)
}

bars(as.matrix(sp[, c("svm_hdstats", "svm_sparseSVM")]),
     c(blue, grey), c("hdsvm", "sparseSVM"), "SVM (hinge) path")

m <- as.matrix(sp[, c("qr_hdstats", "qr_hqreg", "qr_conquer")])
bp <- barplot(t(m), beside = TRUE, col = c(blue, orange, grey), names.arg = sp$lab,
              border = NA, main = "Quantile regression path", ylab = "path time (s)",
              ylim = c(0, max(m) * 1.18), cex.names = 0.85)
legend("topleft", legend = c("hdqr", "hqreg", "conquer"),
       fill = c(blue, orange, grey), bty = "n", border = NA, cex = 0.95)
text(bp[1, ], m[, "qr_hdstats"], pos = 3, offset = 0.3, cex = 0.8, font = 2, col = blue,
     labels = sprintf("%.0f/%.0fx", sp$qr_hqreg / sp$qr_hdstats, sp$qr_conquer / sp$qr_hdstats))

bars(as.matrix(sp[, c("huber_hdstats", "huber_hqreg")]),
     c(blue, grey), c("hdhuber", "hqreg"), "Huber regression path")

mtext("hdstats vs mainstream packages: 100-lambda path time (matched tolerance, median of 5 runs)",
      outer = TRUE, line = -1.4, cex = 0.95, font = 2)
dev.off()
cat("Figure written: results/hdstats_benchmark.png\n")
