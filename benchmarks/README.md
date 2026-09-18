# Benchmarks: hdstats versus mainstream packages

Reproducible comparison of **hdstats** with the mainstream R packages that fit
the same or closely related models. All scripts run from the package root and
write their results to `results/` and figures to `figures/`:

| Script | What it measures |
|---|---|
| `01-speed.R` | Time to compute a 100-value solution path on identical data, single-threaded, median of 3 runs. |
| `02-quality.R` | Optimization accuracy: the value of the penalized objective reached by each solver at the same `lambda` values (no standardization), relative to the exact optimum. |
| `03-accuracy.R` | Statistical performance at the cross-validation-selected `lambda` of each package: estimation error, variable selection, test error, over replicated simulations. |
| `04-figures.R`, `05-tables.R` | Figures and markdown tables from the CSV results. |

Run with `HDSTATS_BENCH_QUICK=1` for a smoke test; `HDSTATS_BENCH_REPS` sets
the number of timing repetitions (speed) or simulation replications (accuracy).

## Competitors and how the problems were matched

| Model | hdstats | Competitors | Notes on matching |
|---|---|---|---|
| Huber regression | `hdhuber()` | **hqreg** (semismooth Newton coordinate descent) | hqreg's Huber loss is ours divided by `gamma`, so `lambda_hqreg = lambda / delta` reproduces the same problem. |
| Quantile regression | `hdqr()` | **hqreg** (`method = "quantile"`), **conquer** (convolution-smoothed QR, LAMM), **quantreg** (`rq(method = "lasso")`, interior-point LP) | Same check loss and `(1/n) sum + lambda * L1` scaling in hqreg and conquer. quantreg's lasso penalizes the summed loss with `lambda/2 * L1`, so `lambda_rq = 2 n lambda`. |
| Linear SVM | `hdsvm()` | **sparseSVM** (huberized hinge, semismooth Newton CD), **gcdnet** (`hhsvm`, Huberized *squared* hinge), **LiblineaR** (`type = 5`, L1-regularized *squared* hinge, one cost value), **quantreg** LP (hinge = check loss at `tau = 1`) | Only sparseSVM targets the same hinge-loss objective; gcdnet and LiblineaR fit different losses and are reported for speed only. sparseSVM codes the first level of `y` as +1, so its weights are sign-flipped before evaluating the objective. |
| Wilcoxon rank regression | `hdrr()` | **Rfit** (unpenalized), **quantreg** on all pairwise differences (unpenalized LP) | No mainstream package fits penalized rank regression; `hdrr()` computes a 100-value path while the competitors compute one unpenalized fit. |
| Least squares (reference) | | **glmnet** | Reference for the statistical comparison. |

`rqPen` was not included: its `alg = "huber"` is hqreg and its `alg = "br"` is
quantreg, both covered above.

## Results

RESULTS_PLACEHOLDER
