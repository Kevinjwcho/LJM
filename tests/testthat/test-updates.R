test_that("LJM_dat complete-case filter (default and explicit)", {
  data("pbc2")
  d <- pbc2
  d$id <- as.numeric(d$id)
  d$Y1 <- log(d$serBilir)
  vl <- list(id = "id", time = "year", EvTime = "years", event = "status2")

  # No h/y_vars -> filter cannot apply -> all at-risk subjects kept.
  base <- LJM_dat(d, s = 5, var_list = vl)
  # Explicit filter.
  filt <- LJM_dat(d, s = 5, var_list = vl,
                  h = 4, min_window_obs = 2, y_vars = "Y1")
  # min_window_obs defaults to 2, so supplying h + y_vars filters by default.
  deflt <- LJM_dat(d, s = 5, var_list = vl, h = 4, y_vars = "Y1")

  base_ids  <- unique(base$Surv_dat$id)
  filt_ids  <- unique(filt$Surv_dat$id)
  deflt_ids <- unique(deflt$Surv_dat$id)

  # filtered cohort is a subset, no larger than the unfiltered one
  expect_true(all(filt_ids %in% base_ids))
  expect_lte(length(filt_ids), length(base_ids))
  # the default (min_window_obs = 2) matches the explicit min_window_obs = 2
  expect_setequal(deflt_ids, filt_ids)

  # every retained subject really has >= 2 non-missing Y1 obs inside |year - s| <= h
  cnt <- vapply(filt_ids, function(i) {
    sub <- d[d$id == i & abs(d$year - 5) <= 4, , drop = FALSE]
    sum(!is.na(sub$Y1))
  }, numeric(1))
  expect_true(all(cnt >= 2))
})

test_that("LJM re-estimates the random-effect mean c (M-step, always on)", {
  data("pbc2")
  d <- pbc2
  d$id <- as.numeric(d$id)
  d$Y1 <- log(d$serBilir)
  vl <- list(id = "id", time = "year", EvTime = "years", event = "status2")
  # complete-case is applied by default here (h + y_vars supplied)
  td <- LJM_dat(d, s = 5, var_list = vl, h = 4, y_vars = "Y1")

  fit <- LJM(td, y_vars = "Y1", s = 5, h = 4, base_terms = c("drug"),
             ker = "epanechnikov", max.iter = 50, Vcov = FALSE, verbose = FALSE)

  expect_s3_class(fit, "LJM")

  # the EM always re-estimates c = mean(b_hat); check it is present and sane
  c_hat <- as.numeric(fit$coefficients$c)
  expect_false(is.null(c_hat))
  expect_true(all(is.finite(c_hat)))
  expect_length(c_hat, 2L)   # intercept + slope for a single marker
})

test_that("AUCdyn/PEdyn fit the censoring KM with one row per subject", {
  data("pbc2")
  d <- pbc2
  d$id <- as.numeric(d$id)
  one <- d[!duplicated(d$id), ]
  ids <- one$id[one$years > 3]
  set.seed(1)
  pr <- data.frame(id = ids, tau = 2, pred_surv = runif(length(ids), 0.5, 1))
  ev <- list(id = "id", EvTime = "years", event = "status2")
  # long format (several rows per subject) must give the same result as one row per subject
  expect_equal(AUCdyn(pr, data = d,   landmarks = 3, tau = 2, var_list = ev),
               AUCdyn(pr, data = one, landmarks = 3, tau = 2, var_list = ev))
  expect_equal(PEdyn(pr, data = d,   landmarks = 3, tau = 2, var_list = ev),
               PEdyn(pr, data = one, landmarks = 3, tau = 2, var_list = ev))
})
