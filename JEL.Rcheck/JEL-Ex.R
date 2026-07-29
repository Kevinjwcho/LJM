pkgname <- "JEL"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('JEL')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("JEL1")
### * JEL1

flush(stderr()); flush(stdout())

### Name: JEL1
### Title: Jointly Estimated Landmarking (JEL) for a Singular Longitudinal
###   Process
### Aliases: JEL1

### ** Examples

# Load example dataset
data("pbc2")
str(pbc2)

# Specify landmarking parameters
s <- 5
u <- 3
tau <- 2

# Example 1: JEL with rho = 0 (Cox proportional hazards model)
testJEL1 <- JEL1(
  pbc2, land_time = s, u = u,
  lme_fixed = list(log(serBilir) ~ year),
  lme_random = list(~ year | id),
  cox_formula = Surv(years, status2) ~ drug,
  Bs = NULL,
  rho = 0
)

# Example 2: JEL with rho = 1 (Proportional odds model)
testJEL2 <- JEL1(
  pbc2, land_time = s, u = u,
  lme_fixed = list(log(serBilir) ~ year),
  lme_random = list(~ year | id),
  cox_formula = Surv(years, status2) ~ drug,
  Bs = NULL,
  rho = 1
)

# Example 3: JEL using time-varying effects with B-splines
testJEL3 <- JEL1(
  pbc2, land_time = s, u = u,
  lme_fixed = list(log(serBilir) ~ year),
  lme_random = list(~ year | id),
  cox_formula = Surv(years, status2) ~ drug,
  Bs = list(df = 3, degree = 1, knots = NULL, Bknots = NULL),
  rho = 0
)

# Predictions of conditional survival probabilities, Pr(T > tau + s | T > s)
prediction_result1 <- predict(testJEL1, tau = tau)
prediction_result2 <- predict(testJEL2)
prediction_result3 <- predict(testJEL3, tau = tau)

# Confidence bands for \eta_0(t) and \eta_1(t) from time-varying model
conint_tv <- confBands.JEL(testJEL3, alpha = 0.05)

# AUC calculations
AUCtest1 <- AUCdyn(
  pred_surv_result = prediction_result1,
  data = pbc2, landmarks = s, tau = tau,
  var_list = list(id = "id", EvTime = "years", event = "status2")
)

# PE calculations
PEtest1 <- PEdyn(
  pred_surv_result = prediction_result1,
  data = pbc2, landmarks = s, tau = tau,
  var_list = list(id = "id", EvTime = "years", event = "status2")
)





cleanEx()
nameEx("JEL_dat")
### * JEL_dat

flush(stderr()); flush(stdout())

### Name: JEL_dat
### Title: Create landmarked longitudinal and survival datasets for
###   JEL/LLAJEL
### Aliases: JEL_dat

### ** Examples

## Not run: 
##D dat_split <- JEL_dat(
##D   data = dat,
##D   land_time = 2,
##D   negative = FALSE,
##D   var_list = list(id="id", time="time", EvTime="Time", event="event")
##D )
##D head(dat_split$LMM_dat)
##D head(dat_split$Surv_dat)
## End(Not run)




cleanEx()
nameEx("LLAJEL")
### * LLAJEL

flush(stderr()); flush(stdout())

### Name: LLAJEL
### Title: Fit an LLA-based Jointly Estimated Landmarking (LLAJEL) model
### Aliases: LLAJEL

### ** Examples

## Not run: 
##D fit <- JEL(
##D   train_dataset = train_dataset,
##D   y_vars = c("Y.1", "Y.2", "Y.3"),
##D   s = 2, h = 2,
##D   base_terms = c("age", "sex")
##D )
##D print(fit$coefficients)
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
