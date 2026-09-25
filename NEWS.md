# LJM 1.0.0

* The package is renamed from `JEL` to `LJM` (local joint model), matching the
  method name in the manuscript. The code is that of `JEL` 2.5.0; estimates and
  predictions are unchanged.
* User-facing names: `JEL()` -> `LJM()`, `JEL_dat()` -> `LJM_dat()`, class
  `"JEL"` -> `"LJM"` (so `predict.JEL()` -> `predict.LJM()` and
  `confBands.JEL()` -> `confBands.LJM()`). `AUCdyn()`, `PEdyn()` and
  `select_h_longitudinal()` keep their names. The old names are not kept as
  aliases; scripts written for `JEL` run unchanged with `JEL` installed, and the
  two packages can be installed side by side.
* Internal helpers keep their historical names (e.g. `InitVal_LLAJEL()`,
  `fitLLAJEL()`); they are not exported.
* The history below refers to the package under its former name.

# JEL 2.5.0

## Transformation models in the local joint model

* `JEL()` gains `rho = 0`, the transformation parameter of
  G(x) = log(1 + rho x)/rho (0: Cox, 1: proportional odds). Before 2.5 the
  local (LLA) estimation path was Cox-only: the rho code in `EMiterJEL`,
  `LHGeneric_JEL` and `LambGeneric_JEL` belonged to the univariate global JEL,
  and `fitMultiJEL` had fixed `rho = 0`. Only `pred_surv_compute_LLA` accepted rho.
* C++: `ll_lla`, `gradll_lla`, `sdll_lla` use the G survival density
  Delta{log l0 + lp - log(1 + rho u)} - G(u), u = Lambda0(V) e^lp, with
  gradient eta (Delta - u)/(1 + rho u) and Hessian
  -u (1 + rho Delta)/(1 + rho u)^2 eta eta'. `sdll_lla` now takes `Delta`.
  `Esurv`/`Esurv_exp` weight the 1-D quadrature by the G density, and
  `Esurv_exp` returns E[xi e^zeta | O] with E(xi | b, O) = (1 + Delta rho)/(1 + rho u)
  (Zeng and Lin 2007), so the Lambda0 update and the (eta, phi) Newton step
  (`Eetaphi`/`Setaphi`/`Hetaphi`) follow automatically. The local variable
  `rho` (number of quadrature nodes) in `Esurv`/`Esurv_exp` was renamed `nq`.
* `rho` is threaded through `RefinedfastEM_LLA`, the PRES standard errors
  (`score_JEL`, `PRES_hessian`) and `predict()`, which now defaults to the
  rho stored in the fit. `rho != 0` stops for the time-varying model.
* Every rho = 0 code path is the unchanged 2.4 code, so results are
  bit-identical to 2.4.

## Verification (2026-09-23)

* PBC, s = 3 and 5, h = 5, Vcov = TRUE: 2.5 (rho = 0) vs 2.4 identical in all
  coefficients, Vcov, b-hat, predictions and iteration counts.
* Same inputs, rho = 1e-8 vs 0: differences ~1e-9 in ll, gradient and Hessian.
* rho = 0.5, 1: analytic gradient/Hessian match numerical derivatives
  (relative error 1e-10 to 1e-11); quadrature matches numerical integration
  (1e-10 to 1e-12).
* rho = 0 on the general G path (rho = 1e-8) costs the same per EM iteration.
* Known, not addressed: on PBC the PRES/ridge standard errors are unstable
  when the information matrix is near singular (s = 5: serBilir intercept SE
  0.598 at rho = 0 vs 2.444 at rho = 1e-8 with near-identical estimates).
* `fit$coefficients$beta` holds base R's `beta()` function (harmless, cosmetic).

# JEL 2.4.0

## Changed

* `Predict_BLUP_LLA()` now computes the MAP/BLUP predictor in the precision form
  `b_hat = (D^-1 + Z' Omega Z)^-1 (D^-1 c + Z' Omega y)` with
  `Omega = diag(K_h(t-s) / sigma_k^2)`, matching Section 3.1 of the manuscript.
  The previous covariance form `c + D Z'(Omega^-1 + Z D Z')^-1 (y - Z c)` is
  algebraically identical (Woodbury) but inverted an `N_i x N_i` matrix and had
  to guard `Omega^-1` with `w + 1e-12`, which put entries of order `1e12` on the
  diagonal for zero-weight observations. Verified on PBC at s in {3,5} and
  h in {3,4,5}: `b_hat` identical to machine zero, predicted survival within
  1.4e-12.

* Initial values: the pooled kernel linear mixed model of Step 1 is now fitted at
  the SAME bandwidth `h` as the EM (pooling across subjects makes the wider
  window unnecessary). If that fit fails or returns a zero variance component,
  it falls back to the full-record fit of 2.3, and then to per-subject WLS.
  `fit$init_h` and `fit$init_mode` report what was used.

# JEL 2.3.0

## Bug fixes

* `PRES_hessian()` (the finite-difference Hessian behind `Vcov = TRUE`) now
  scales its step to the parameter magnitude and, for the entries of `D`,
  halves it until every stencil point keeps `D` positive-definite. With the
  previous absolute step of `1e-4`, a near-zero random-slope variance put a
  stencil point at a singular `D`, giving Hessian entries of order `1e15` and,
  through the ridged inverse, standard errors of `eta` near `1000`. Point
  estimates are unaffected.

* `RefinedfastEM_LLA()`: the (eta, phi) Newton update is now guarded. A
  non-finite score/Hessian (exp overflow on small risk sets or badly scaled
  markers) restarts (eta, phi) from zero once, later occurrences skip the
  update, and a step larger than 10 in any coordinate is damped. The number of
  interventions is returned as `fit$newton_guard` (0 for all manuscript fits).
* `InitVal_LLAJEL()`: if the Cox-based initial (eta, phi) give a centred linear
  predictor wider than 10, they are shrunk proportionally before the EM starts;
  the factor is returned as `fit$init_shrink` (1 = no intervention).

## Breaking changes

* `JEL()` no longer takes `h_init`. Initial values (the pooled kernel linear
  mixed model giving `c`, `Sigma_b(s)`, `sigma^2` and the BLUPs of `b_i(s)`,
  and the Cox fit seeding `lambda_0`, `phi`, `eta`) are now always built from
  the **full longitudinal record** of each retained subject, using the same
  kernel `ker` as the EM fit at a bandwidth wide enough to include every
  observation (`full_range_h()`). The kernel shape is preserved -- weights are
  largest at the landmark and taper towards the edges -- but nothing is
  excluded, and the initialization bandwidth is no longer user-tunable.
* Previously `h_init = NULL` (the default) reused the fitting bandwidth `h` for
  initialization, so the initial values were built from a localized sub-window
  and observations outside `|t - s| <= h` were discarded. Results from 2.2.x
  are therefore **not** reproduced by 2.3.0; anything fitted with the old
  default has to be refitted.

# JEL 2.2.1

## Bug fixes

* `predict.JEL()` now drops longitudinal rows of `testdat` observed after the
  landmark time `s` before computing the BLUPs (new argument
  `newdata_past_only = TRUE`). Since the 2.0 refactor, `predict()` routed new
  data through the same `JEL_dat()` used for estimation, which keeps all
  observations for the symmetric kernel; as a result, held-out predictions
  could condition on post-landmark measurements. Dynamic predictions must
  condition only on the history up to `s` (JEL <= 1.2 enforced this via a hard
  `t <= s` cutoff in `JEL_dat()`). Estimation is unchanged. Set
  `newdata_past_only = FALSE` to reproduce the old behaviour for diagnostics.

# JEL 2.2

## Renamed

* The main fitting function `LLAJEL()` is renamed to `JEL()`, matching the
  package and method name; `JEL()`, `JEL_dat()`, and `predict.JEL()` now form a
  consistent interface. The old name `LLAJEL()` is no longer available.

## Removed

* Removed the legacy `JEL1()` entry point (the old split-API version). The
  current interface is `JEL()` together with `JEL_dat()`, as used throughout the
  manuscript; nothing in the package or analyses depended on `JEL1()`.

## Behaviour changes

* The population mean of the random effects, `c`, is now **always** re-estimated
  at each EM iteration via the M-step update `c_hat = mean(b_hat_i)`, matching
  the manuscript. This is an internal detail and is no longer exposed as a
  toggle (the former `update_c` argument has been removed).
* `JEL_dat()` applies the complete-case filter **by default**: `min_window_obs`
  now defaults to `2`, so when both `h` and `y_vars` are supplied a subject is
  kept only if it has at least `min_window_obs` non-missing observations of every
  marker in `y_vars` inside the kernel window `|t - s| <= h`. The filter is
  skipped when `h` or `y_vars` is absent, and can be disabled with
  `min_window_obs = NULL`.

## New features

* `JEL()` gains `h_init`: an optional bandwidth used only to build the initial
  local-linear BLUPs of the random effects, decoupled from the fitting bandwidth
  `h`. Defaults to `h`.
* `JEL_dat()` gains `h` and `y_vars`, used together with `min_window_obs` to
  define the complete-case filter described above.

## Documentation

* Documented the previously undocumented `Bs` argument of `JEL()` and the new
  arguments above; regenerated the corresponding help pages.

## Tests

* Added `test-updates.R` covering the default and explicit complete-case filter
  in `JEL_dat()` and the M-step re-estimation of `c` in `JEL()`.
