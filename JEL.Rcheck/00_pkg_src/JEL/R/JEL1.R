#' Jointly Estimated Landmarking (JEL) for a Singular Longitudinal Process
#'
#' The \code{JEL1} function implements a Jointly Estimated Landmarking approach to predict survival probabilities at specified landmark times.
#' It integrates longitudinal mixed-effects modeling and Cox proportional hazards modeling, providing dynamic predictions
#' for survival outcomes based on past longitudinal data.
#' @param data A data frame containing both longitudinal and survival data. The data must include the following columns: id, time variable for the longitudinal data, survival event time, and survival event indicator (censor). All these variables must be referenced in \code{lme_fixed}, \code{lme_random}, and \code{cox_formula}. Refer to \code{pbc2} as an example.
#' @param land_time A numeric value specifying the landmarking time point for dynamic survival prediction.
#' @param u A numeric value specifying the interval length of past observations to include in the modeling.
#' @param lme_fixed A list containing fixed-effects formulas for the longitudinal mixed-effects models.
#'
#'   - Must contain exactly one formula of fixed effect of linear mixed model.
#' @param lme_random A list containing random-effects formulas for the longitudinal mixed-effects models.
#'
#'   - Must contain exactly one formula of the form \code{~ time | id}.
#' @param cox_formula A formula for the Cox proportional hazards model, specifying the survival outcome and covariates.
#'
#'   - Must follow the format \code{Surv(time, event) ~ covariates}.
#' @param Bs An optional list specifying B-spline parameters for time-varying effects, with components:
#'
#'   - \code{df}: Degrees of freedom.
#'
#'   - \code{degree}: Degree of the spline.
#'
#'   - \code{knots}: Locations of internal knots.
#'
#'   - \code{Bknots}: Boundary knots.
#'
#'   If \code{NULL}, the landmarking model is fit without time-varying effects (default: \code{NULL}).
#' @param rho A numeric value specifying the penalty parameter for the landmarking model (default: \code{0}).
#'
#' @details
#' The Jointly Estimated Landmarking (JEL) framework dynamically combines longitudinal and survival data to estimate survival probabilities
#' at specified landmark times. It involves preparing data, fitting longitudinal mixed-effects models and Cox proportional hazards models,
#' and integrating these into a unified framework.
#'
#' @return A list of class \code{"JEL"} containing:
#' \describe{
#'   \item{\code{fitLME_list}}{A list of fitted longitudinal mixed-effects models, with estimated fixed and random effects.}
#'   \item{\code{fitCOX}}{The fitted Cox proportional hazards model, including coefficients, baseline hazard estimates, and diagnostics.}
#'   \item{\code{fitJEL}}{The fitted landmarking model that combines results from the longitudinal and survival models, with sub-items:
#'     \describe{
#'       \item{\code{coefficients}}{A list containing:
#'         \describe{
#'           \item{\code{beta}}{Estimated fixed effects for the longitudinal process.}
#'           \item{\code{phi}}{Estimated baseline survival covariates.}
#'           \item{\code{eta}}{Estimated association paramters between longitudinal and survival processes.}
#'           \item{\code{Ysigma}}{Residual variance of the longitudinal process.}
#'           \item{\code{Bsigma}}{Random effects variance-covariance matrix.}
#'           \item{\code{lamb}}{Estimated baseline hazard function over time.}
#'         }
#'       }
#'       \item{\code{logLik}}{Log-likelihood of the fitted landmarking model.}
#'       \item{\code{numIter}}{Number of EM algorithm iterations performed.}
#'       \item{\code{Vcov}}{Variance-covariance matrix of the parameter estimates.}
#'       \item{\code{est.bi}}{Predicted random effects for each subject.}
#'       \item{\code{convergence}}{Indicator of convergence status ("success" or "failure").}
#'       \item{\code{control}}{Control parameters used for model fitting, including tolerance and maximum iterations.}
#'       \item{\code{time.SE}}{Computation time for standard error estimation.}
#'       \item{\code{N}}{Total number of observations in the longitudinal dataset.}
#'       \item{\code{n}}{Number of unique subjects in the dataset.}
#'       \item{\code{d}}{Event indicator vector for survival data (1 if the event occurred, 0 otherwise).}
#'       \item{\code{rho}}{Penalty parameter used in the landmarking model.}
#'       \item{\code{dataMat}}{A list summarizing the prepared data, including:
#'         \describe{
#'           \item{\code{B}}{Basis matrix for time-varying effects (if applicable).}
#'           \item{\code{ID}}{Subject ID vector.}
#'         }
#'       }
#'     }
#'   }
#'   \item{\code{dataset}}{A list of prepared datasets used in the JEL model:
#'     \describe{
#'       \item{\code{LMM_dat}}{The dataset formatted for the longitudinal mixed-effects models.
#'        And the data is transformed the time variable given \code{s} and \code{u}.}
#'       \item{\code{Surv_dat}}{The dataset formatted for the Cox proportional hazards model. And the data involves start, stop, and event variable for JEL model.}
#'       }
#'    }
#' }
#'
#' @examples
#' # Load example dataset
#' data("pbc2")
#' str(pbc2)
#'
#' # Specify landmarking parameters
#' s <- 5
#' u <- 3
#' tau <- 2
#'
#' # Example 1: JEL with rho = 0 (Cox proportional hazards model)
#' testJEL1 <- JEL1(
#'   pbc2, land_time = s, u = u,
#'   lme_fixed = list(log(serBilir) ~ year),
#'   lme_random = list(~ year | id),
#'   cox_formula = Surv(years, status2) ~ drug,
#'   Bs = NULL,
#'   rho = 0
#' )
#'
#' # Example 2: JEL with rho = 1 (Proportional odds model)
#' testJEL2 <- JEL1(
#'   pbc2, land_time = s, u = u,
#'   lme_fixed = list(log(serBilir) ~ year),
#'   lme_random = list(~ year | id),
#'   cox_formula = Surv(years, status2) ~ drug,
#'   Bs = NULL,
#'   rho = 1
#' )
#'
#' # Example 3: JEL using time-varying effects with B-splines
#' testJEL3 <- JEL1(
#'   pbc2, land_time = s, u = u,
#'   lme_fixed = list(log(serBilir) ~ year),
#'   lme_random = list(~ year | id),
#'   cox_formula = Surv(years, status2) ~ drug,
#'   Bs = list(df = 3, degree = 1, knots = NULL, Bknots = NULL),
#'   rho = 0
#' )
#'
#' # Predictions of conditional survival probabilities, Pr(T > tau + s | T > s)
#' prediction_result1 <- predict(testJEL1, tau = tau)
#' prediction_result2 <- predict(testJEL2)
#' prediction_result3 <- predict(testJEL3, tau = tau)
#'
#' # Confidence bands for \eta_0(t) and \eta_1(t) from time-varying model
#' conint_tv <- confBands.JEL(testJEL3, alpha = 0.05)
#'
#' # AUC calculations
#' AUCtest1 <- AUCdyn(
#'   pred_surv_result = prediction_result1,
#'   data = pbc2, landmarks = s, tau = tau,
#'   var_list = list(id = "id", EvTime = "years", event = "status2")
#' )
#'
#' # PE calculations
#' PEtest1 <- PEdyn(
#'   pred_surv_result = prediction_result1,
#'   data = pbc2, landmarks = s, tau = tau,
#'   var_list = list(id = "id", EvTime = "years", event = "status2")
#' )
#'
#'
#' @import dplyr
#' @importFrom nlme lme ranef fixef
#' @importFrom survival coxph Surv survfit
#' @import splines
#' @import statmod
#' @import ucminf
#' @importFrom Rcpp sourceCpp
#' @useDynLib JEL, .registration = TRUE
#' @export
JEL1 <- function(data, land_time, u, lme_fixed, lme_random, cox_formula, Bs = NULL, rho = 0, time_var) {

  # Ensure lme_fixed and lme_random have the same length
  if (!(length(lme_fixed) == length(lme_random) & length(lme_random) == 1)) {
    stop("The number of fixed effects formulas (lme_fixed) and random effects formulas (lme_random) must be one.")
  }

  if(is.null(time_var)){
    stop("Please specify the time variable for the longitudinal data.")
  }

  # Extract id_var and time_var from lme_random
  id_var <- NULL
  if (length(lme_random) > 0) {
    first_random <- lme_random[[1]]
    random_parts <- nlme::splitFormula(first_random, "|")
    if (length(random_parts) == 2) {
      id_var <- as.character(random_parts[[2]])[2]
    } else {
      stop("Invalid random effects formula. Ensure the formula uses '|' (e.g., ~ time | id).")
    }
  }

  # Extract event_var and surv_time_var from cox_formula
  event_var <- NULL
  surv_time_var <- NULL
  if (!is.null(cox_formula)) {
    surv_vars <- all.vars(cox_formula[[2]])
    if (length(surv_vars) >= 3) {
      surv_time_var <- surv_vars[1:2]
      event_var <- surv_vars[3]
    } else if (length(surv_vars) == 2) {
      surv_time_var <- surv_vars[1]
      event_var <- surv_vars[2]
    }
  }

  # Data Preparation
  dataset <- JEL_dat(
    data = data,
    land_time = land_time,
    var_list = list(
      id = id_var, time = time_var,
      EvTime = surv_time_var, event = event_var
    )
  )

  # Fit Longitudinal Mixed-Effects Models
  fitLME_list <- mapply(function(fixed, random) {
    lme(
      fixed,
      data = dataset$LMM_dat,
      random = list(as.formula(random)),
      method = "REML",
      control = list(opt = "optim"),
      na.action = "na.exclude"
    )
  }, lme_fixed, lme_random, SIMPLIFY = FALSE)

  # fitLME <- lme(
  #   fixed = lme_fixed[[1]],
  #   data = dataset$LMM_dat,
  #   random = list(as.formula(lme_random[[1]])),
  #   method = "REML",
  #   control = list(opt = "optim"),
  #   na.action = "na.exclude"
  # )

  # Build the Cox model formula dynamically
  rhs <- as.character(cox_formula)[3]  # Extract right-hand side of cox_formulas
  surv_formula <- as.formula(paste("Surv(start, stop, event) ~", rhs))

  # Fit Survival Model
  fitCOX <- coxph(
    surv_formula,
    data = dataset$Surv_dat,
    x = TRUE,
    model = TRUE,
    id = dataset$Surv_dat[[id_var]],      # Use id_var to extract the column
    cluster = dataset$Surv_dat[[id_var]] # Use id_var for clustering
  )

  # Fit Multi-Joint EL Model (Standard or Time-Varying)
  if (is.null(Bs)) {
    fit.JEL = fitJEL(
      fitLME_list[[1]],
      fitCOX,
      dataset$Surv_dat,
      model = "Fixed",
      rho = rho,
      control = list(SE.method = "PRES")
    )
  } else {
    df = Bs$df
    degree = Bs$degree
    knots = Bs$knots
    Bknots = Bs$Bknots

    fit.JEL <- fitJEL(
      fitLME_list[[1]],
      fitCOX,
      dataset$Surv_dat,
      df = df,
      degree = degree,
      knots = knots,
      Bknots = Bknots,
      model = "TimeVar",
      rho = rho,
      control = list(SE.method = "PRES")
    )
  }

  # Add class and return result
  result <- list(
    fitLME_list = fitLME_list,
    fitCOX = fitCOX,
    fitJEL = fit.JEL,
    dataset = dataset
  )
  class(result) <- "JEL"
  return(result)
}
