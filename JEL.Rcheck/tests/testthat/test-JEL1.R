test_that("JEL1 works well", {
  data("pbc2")
  str(pbc2)

  s = 5
  u = 3
  tau = 2

  # test the functions
  testJEL1 = JEL1(pbc2, land_time = s, u = u,
                  lme_fixed = list(
                    log(serBilir) ~ year
                  ),
                  lme_random = list(
                    ~ year | id
                  ),
                  cox_formula = Surv(years, status2) ~ drug,
                  Bs = NULL,
                  rho = 0
  )

  testJEL2 = JEL1(pbc2, land_time = s, u = u,
                  lme_fixed = list(
                    log(serBilir) ~ year
                  ),
                  lme_random = list(
                    ~ year | id
                  ),
                  cox_formula = Surv(years, status2) ~ drug,
                  Bs = NULL,
                  rho = 1
  )

  testJEL3 = JEL1(pbc2, land_time = s, u = u,
                  lme_fixed = list(
                    log(serBilir) ~ year
                  ),
                  lme_random = list(
                    ~ year | id
                  ),
                  cox_formula = Surv(years, status2) ~ drug,
                  Bs = list(df = 3, degree = 1, knots = NULL, Bknots = NULL),
                  rho = 0
  )


  # print(testJEL1)
  # print(testJEL2)

  prediction_result1 = predict(testJEL1, tau = tau)
  prediction_result1 = predict(testJEL1, tau = tau, CI = TRUE)
  prediction_result2 = predict(testJEL2)
  prediction_result3 = predict(testJEL3, tau = tau)

  # conint_tv = confBands.JEL(testJEL1, alpha = 0.05)
  # conint_tv = confBands.JEL(testJEL2, alpha = 0.05)
  conint_tv = confBands.JEL(testJEL3, alpha = 0.05)

  AUCtest1 = AUCdyn(pred_surv_result = prediction_result1,
                    data = pbc2, landmarks = s, tau = tau,
                    var_list = list(id = "id", EvTime = "years", event = "status2")
  )

  AUCtest2 = AUCdyn(pred_surv_result = prediction_result2,
                    data = pbc2, landmarks = s, tau = tau,
                    var_list = list(id = "id", EvTime = "years", event = "status2")
  )

  AUCtest3 = AUCdyn(pred_surv_result = prediction_result3,
                    data = pbc2, landmarks = s, tau = tau,
                    var_list = list(id = "id", EvTime = "years", event = "status2")
  )

  PEtest1 = PEdyn(pred_surv_result = prediction_result1,
                  data = pbc2, landmarks = s, tau = tau,
                  var_list = list(id = "id", EvTime = "years", event = "status2")
  )

  PEtest2 = PEdyn(pred_surv_result = prediction_result2,
                  data = pbc2, landmarks = s, tau = tau,
                  var_list = list(id = "id", EvTime = "years", event = "status2")
  )

  PEtest3 = PEdyn(pred_surv_result = prediction_result3,
                  data = pbc2, landmarks = s, tau = tau,
                  var_list = list(id = "id", EvTime = "years", event = "status2")
  )
})


