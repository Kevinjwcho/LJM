# JEL 2.2

## New features

* `LLAJEL()` gains `update_c` (default `TRUE`): the population mean of the
  random effects, `c`, is now re-estimated at each EM iteration via the M-step
  update `c_hat = mean(b_hat_i)`, matching the manuscript. Set `update_c = FALSE`
  for the legacy behaviour that freezes `c` at its initial value.
* `LLAJEL()` gains `h_init`: an optional bandwidth used only to build the initial
  local-linear BLUPs of the random effects, decoupled from the fitting bandwidth
  `h`. Defaults to `h`.
* `JEL_dat()` gains `h`, `min_window_obs`, and `y_vars`: an optional
  complete-case filter that keeps a subject only if it has at least
  `min_window_obs` non-missing observations of every marker in `y_vars` inside
  the kernel window `|t - s| <= h`.

## Documentation

* Documented the previously undocumented `Bs` argument of `LLAJEL()` and the new
  arguments above; regenerated the corresponding help pages.

## Tests

* Added `test-updates.R` covering the `min_window_obs` complete-case filter and
  the `update_c` M-step re-estimation.
