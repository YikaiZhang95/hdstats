## ---------------------------------------------------------------------------
## Install the packages needed to reproduce the hdstats benchmark.
##
## hdstats itself is compared against the closest mainstream high-dimensional
## path solvers:
##   sparseSVM  -- sparse penalized linear SVM  (path)         vs hdsvm
##   hqreg      -- penalized Huber & quantile regression (path) vs hdhuber, hdqr
##   conquer    -- convolution-smoothed penalized quantile reg. vs hdqr
## and against an exact convex reference for every model:
##   CVXR       -- disciplined convex programming (ECOS solver) vs all four
## quantreg / glmnet / e1071 / MASS are optional classical references.
## ---------------------------------------------------------------------------
pkgs <- c("Rcpp", "Matrix", "sparseSVM", "hqreg", "conquer", "CVXR",
          "quantreg", "glmnet", "e1071", "MASS")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) install.packages(need, repos = "https://cloud.r-project.org")

## install hdstats from the package root (run from the repo top level), e.g.:
##   R CMD INSTALL .
## or
##   remotes::install_github("YikaiZhang95/hdstats")

invisible(lapply(pkgs, function(p)
  cat(sprintf("%-12s %s\n", p,
      if (requireNamespace(p, quietly = TRUE))
        as.character(packageVersion(p)) else "MISSING"))))
