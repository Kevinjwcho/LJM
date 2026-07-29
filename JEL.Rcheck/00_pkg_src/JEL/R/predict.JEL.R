#' Predict Conditional Survival Probabilities Using JEL Models
#'
#' This function predicts conditional survival probabilities using JEL (Joint Estimation with Landmarks) models.
#'
#' @param JEL_models A list from the \code{JEL1} or \code{JEL2}, including:
#'
#'   - \code{fitJEL}: The fitted Joint Estimation Landmarking (JEL) model.
#'
#'   - \code{fitLME_list}: A list of fitted longitudinal linear mixed-effects (LME) models.
#'
#'   - \code{fitCOX}: The fitted Cox proportional hazards model.
#'
#'   - \code{dataset}: A list of training and testing datasets, with the following components:
#'      \itemize{
#'          \item \code{LMM_dat}: Data used for the longitudinal mixed-effects models.
#'
#'          \item \code{Surv_dat}: Data used for the survival Cox model.
#'      }
#'
#'     - \code{s}: Landmark time used for model training.
#'
#'     - \code{u}: Interval length of past observations included in the modeling.
#'
#'     - \code{var_list}: A list of variable names for identification (\code{id}), time (\code{time}), survival event time (\code{EvTime}), and event indicator (\code{event}).
#' @param testdat An optional data frame containing new data for prediction. If \code{NULL}, the training dataset from \code{JEL_models} is used (default: \code{NULL}).
#' @param tau A numeric value specifying the prediction time horizon. If \code{NULL}, the full conditional survival probabilities are returned (default: \code{NULL}).
#' @param CI A logical value indicating whether to compute confidence intervals for the predicted survival probabilities (default: \code{FALSE}). The confidence interval is obtained using Monte Carlo approach.
#' @param MC An integer specifying the number of Monte Carlo simulations for confidence interval estimation (default: 100).
#' @param alpha A numeric value specifying the significance level for confidence interval estimation (default: 0.05).
#'
#' @details
#' This function computes conditional survival probabilities for individuals using the components of a JEL model. It estimates the Best Linear Unbiased Predictors (BLUPs) for random effects based on the new or training dataset and combines them with survival predictions from the fitted Cox model. If a prediction time horizon (`tau`) is specified, the function filters the results to include survival probabilities at or before the specified time horizon.
#'
#' The workflow involves:
#'
#' 1. BLUP Estimation: Predicts the random effects for new data using the fitted JEL and LME models.
#'
#' 2. Survival Probability Computation: Uses the Cox model and the estimated BLUPs to calculate conditional survival probabilities.
#'
#' 3. Filtering by \code{tau} (if provided): Returns survival probabilities up to the specified time horizon.
#'
#'
#' @return A data frame or list containing the predicted survival probabilities. The structure of the output depends on whether \code{tau} is specified:
#'
#'   - If \code{tau} is \code{NULL}, the function returns a full prediction for all time points.
#'
#'   - If \code{tau} is specified, the function returns predictions filtered at or before the given time horizon.
#'
#' @export

predict.JEL <- function(fitLLAJEL,
                           testdat = NULL,
                           tau = NULL,
                           CI = FALSE,
                           MC = 100,
                           alpha = 0.05,
                           rho = 0,
                           bLLA_new = NULL) {
  
  # -----------------------------
  # Checks
  # -----------------------------
  if (is.null(fitLLAJEL$fitCOX)) stop("fitLLAJEL must contain $fitCOX.")
  if (is.null(fitLLAJEL$prep)) stop("fitLLAJEL must contain $prep (from LLAJEL).")
  if (is.null(fitLLAJEL$dataset)) {
    # 너 fitLLAJEL 리턴에서 dataset 저장 안 했으면 여기서 안내
    # (JEL() 리턴에 dataset=train_dataset 넣으면 가장 깔끔)
    if (is.null(testdat)) stop("fitLLAJEL$dataset is missing. Provide testdat explicitly.")
  }
  
  # -----------------------------
  # Build Cox test dataset
  # -----------------------------
  if (is.null(testdat)) {
    Cox_dat.te <- fitLLAJEL$dataset$Surv_dat
    LMM_dat.te <- fitLLAJEL$dataset$LMM_dat
    bLLA_new <- fitLLAJEL$est.bi
  } else {
    # testdat = raw long+surv merged data.frame
    if (is.null(fitLLAJEL$dataset$var_list)) stop("fitLLAJEL$dataset$var_list is missing.")
    if (is.null(fitLLAJEL$dataset$s)) stop("fitLLAJEL$dataset$s is missing (landmark).")
    
    test_ds <- JEL_dat(
      data = testdat,
      land_time = fitLLAJEL$dataset$s,
      var_list = fitLLAJEL$dataset$var_list
    )
    Cox_dat.te <- test_ds$Surv_dat
    LMM_dat.te <- test_ds$LMM_dat
  }
  
  # -----------------------------
  # Build bBLUP_new if not provided
  # -----------------------------
  if (is.null(bLLA_new)) {
    
    s0  <- fitLLAJEL$prep$s0
    h   <- fitLLAJEL$prep$h
    ker <- fitLLAJEL$prep$ker
    
    prep_te <- prep_LLA_landmark(
      LMM_dat  = LMM_dat.te,
      Surv_dat = Cox_dat.te,
      y_vars   = fitLLAJEL$prep$y_vars,
      s0       = s0,
      h        = h,
      ker      = ker,
      var_list = fitLLAJEL$dataset$var_list
    )
    
    # --- kernel-aware BLUP (posterior mean) ---
    blup_te <- Predict_BLUP_LLA(prep_te = prep_te, fitLLAJEL = fitLLAJEL)
    
    # blup_te$mean: (n_subject x 2*nK) matrix
    bLLA_new <- blup_te$mean
    
    # rownames를 subject id로 강제 (뒤에서 매칭하려고)
    if (is.null(rownames(bLLA_new))) {
      rownames(bLLA_new) <- as.character(sort(unique(Cox_dat.te$id)))
    }
  }
  
  # -----------------------------
  # Predict conditional survival
  # -----------------------------
  pred <- pred_surv_compute_LLA(
    fitLLAJEL = fitLLAJEL,
    Cox_dat.te = Cox_dat.te,
    bLLA_new = bLLA_new,
    tau = tau,
    rho = rho,
    CI = CI,
    MC = MC,
    alpha = alpha
  )
  
  pred
}