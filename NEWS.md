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
