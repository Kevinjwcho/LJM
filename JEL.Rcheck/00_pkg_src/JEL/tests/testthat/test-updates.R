test_that("JEL_dat min_window_obs applies a complete-case filter", {
  data("pbc2")
  d <- pbc2
  d$id <- as.numeric(d$id)
  d$Y1 <- log(d$serBilir)
  vl <- list(id = "id", time = "year", EvTime = "years", event = "status2")

  base <- JEL_dat(d, land_time = 5, var_list = vl)
  filt <- JEL_dat(d, land_time = 5, var_list = vl,
                  h = 4, min_window_obs = 2, y_vars = "Y1")

  base_ids <- unique(base$Surv_dat$id)
  filt_ids <- unique(filt$Surv_dat$id)

  # the filtered cohort is a subset and no larger than the unfiltered one
  expect_true(all(filt_ids %in% base_ids))
  expect_lte(length(filt_ids), length(base_ids))

  # every retained subject really has >= 2 non-missing Y1 obs inside |year - s| <= h
  cnt <- vapply(filt_ids, function(i) {
    sub <- d[d$id == i & abs(d$year - 5) <= 4, , drop = FALSE]
    sum(!is.na(sub$Y1))
  }, numeric(1))
  expect_true(all(cnt >= 2))
})

test_that("update_c re-estimates the random-effect mean c (M-step)", {
  data("pbc2")
  d <- pbc2
  d$id <- as.numeric(d$id)
  d$Y1 <- log(d$serBilir)
  vl <- list(id = "id", time = "year", EvTime = "years", event = "status2")
  td <- JEL_dat(d, land_time = 5, var_list = vl)

  fit_upd <- JEL(td, y_vars = "Y1", s = 5, h = 4, base_terms = c("drug"),
                    ker = "epanechnikov", max.iter = 50, Vcov = FALSE,
                    verbose = FALSE, update_c = TRUE)
  fit_frz <- JEL(td, y_vars = "Y1", s = 5, h = 4, base_terms = c("drug"),
                    ker = "epanechnikov", max.iter = 50, Vcov = FALSE,
                    verbose = FALSE, update_c = FALSE)

  expect_s3_class(fit_upd, "JEL")
  expect_s3_class(fit_frz, "JEL")

  c_upd <- fit_upd$coefficients$c
  c_frz <- fit_frz$coefficients$c
  expect_false(is.null(c_upd))
  expect_false(is.null(c_frz))

  # re-estimating c should move it away from the frozen initial value
  expect_false(isTRUE(all.equal(as.numeric(c_upd), as.numeric(c_frz))))
})
