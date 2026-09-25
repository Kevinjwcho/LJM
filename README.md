# LJM

**Local Joint Models for Dynamic Survival Prediction**

`LJM` predicts survival risk dynamically from noisy longitudinal biomarkers.
Around a landmark time `s`, each biomarker trajectory is summarised by a local
linear approximation — its recent **level** (intercept) and recent **change**
(slope) — and a flexible transformation survival model relates these to the
conditional risk. The associations can be constant or **time-varying** (via
B-splines), and multiple biomarkers are handled jointly with a fast EM
algorithm.

## Installation

```r
# from GitHub
# install.packages("remotes")
remotes::install_github("kevinjwcho/LJM")

# or from a local source tree
# install.packages("path/to/LJM", repos = NULL, type = "source")
```

## Quick start

```r
library(LJM)

data("pbc2", package = "LJM")
d <- pbc2
d$id  <- as.numeric(d$id)
d$Y.1 <- log(d$serBilir)   # biomarker 1
d$Y.2 <- d$albumin         # biomarker 2
y_vars   <- c("Y.1", "Y.2")
var_list <- list(id = "id", time = "year", EvTime = "years", event = "status2")

## 1) landmarked dataset at s = 5 (complete-case filter on by default)
td <- LJM_dat(d, s = 5, var_list = var_list, h = 4, y_vars = y_vars)

## 2) fit (Epanechnikov kernel, fixed bandwidth h = 4)
fit <- LJM(td, y_vars = y_vars, s = 5, h = 4,
           base_terms = "drug", ker = "epanechnikov")
fit$coefficients$eta

## 3) predict conditional risk at horizon tau = 2, then evaluate
pr <- predict(fit, testdat = d, tau = 2)
ev <- list(id = "id", EvTime = "years", event = "status2")
AUCdyn(pr, data = d, landmarks = 5, tau = 2, var_list = ev)
PEdyn (pr, data = d, landmarks = 5, tau = 2, var_list = ev)
```

### Time-varying associations

```r
spl    <- list(df = 3, degree = 1, knots = NULL, Bknots = NULL)
fit_tv <- LJM(td, y_vars = y_vars, s = 5, h = 4, base_terms = "drug",
              ker = "epanechnikov", Bs = list(spl, list(NULL)))
cb <- confBands.LJM(fit_tv, K = 1)   # pointwise confidence bands, marker 1
```

## Main functions

| Function | Purpose |
|---|---|
| `LJM_dat()` | Build the landmarked longitudinal + survival dataset at time `s`. |
| `LJM()` | Fit the model (fixed or time-varying via `Bs`). |
| `predict()` | Conditional survival / risk at horizon `tau`. |
| `AUCdyn()`, `PEdyn()` | Time-dependent AUC and prediction error. |
| `confBands.LJM()` | Confidence bands for time-varying coefficient functions. |
| `select_h_longitudinal()` | Cross-validated bandwidth selection. |

## Reference

Cho, J. and Chen, K. (2026). *Survival Prediction with Multiple Longitudinal
Biomarkers: Transformation Models and a Local Likelihood Approach.* Manuscript.

`LJM` was developed under the name `JEL` (versions up to 2.5.0); see `NEWS.md`.

## License

GPL-3
