#' Confidence Bands for Time-Varying Coefficients in LJM Models
#'
#' This function constructs confidence bands for the time-varying coefficients of a local joint model (LJM) model.
#' It requires a B-spline configuration (\code{Bs}) to compute the confidence intervals for both slope and intercept terms.
#'
#' @param object A fitted time-varying LJM model object (class \code{"LJM"}), i.e. one fit with a B-spline basis (\code{Bs}).
#' @param K A numeric value indicating the index of the longitudinal process for which confidence bands are computed.
#'   - Default is 1, which applies to single-process models.
#' @param alpha A numeric value specifying the significance level for the confidence bands. Default is 0.05 (95\% confidence level).
#'
#' @details
#' The \code{confBands.LJM} function computes confidence bands for the time-varying slope and intercept coefficients in an LJM model using the following steps:
#'
#' 1. B-Spline Configuration:
#'
#'    - Extracts the B-spline configuration (\code{Bs}) from the fitted LJM model to construct the design matrix for the time variable.
#'
#'    - Ensures the spline configuration is properly indexed for multi-process models.
#'
#' 2. Variance-Covariance Matrix Adjustment:
#'
#'    - Performs eigen decomposition of the variance-covariance matrix to ensure it is positive semi-definite.
#'
#'    - Sets any negative eigenvalues to zero before recomputing the matrix.
#'
#' 3. Confidence Band Computation:
#'
#'    - Computes point estimates for the time-varying coefficients (slope and intercept) based on the B-spline design matrix.
#'
#'    - Uses the adjusted variance-covariance matrix to calculate the standard errors.
#'
#'    - Constructs confidence intervals for each coefficient using the chi-square distribution.
#'
#'    - Confidence bands are calculated based on the chi-square distribution, with critical values determined by the specified significance level (\code{alpha}) and the degrees of freedom (\code{df}) derived from the B-spline configuration.
#'
#' @return A list containing two data frames:
#'
#'   - \code{intercept}: A data frame with columns \code{tau} (prediction time after the landmarking time point), \code{est} (point estimate), \code{lower} (lower confidence limit), and \code{upper} (upper confidence limit) for the intercept term.
#'
#'   - \code{slope}: A data frame with the same structure as \code{intercept}, providing the confidence bands for the slope term.
#'
#' @examples
#' \donttest{
#' data("pbc2", package = "LJM")
#' d <- pbc2; d$id <- as.numeric(d$id)
#' d$Y.1 <- log(d$serBilir); d$Y.2 <- d$albumin
#' vl <- list(id = "id", time = "year", EvTime = "years", event = "status2")
#' td  <- LJM_dat(d, s = 5, var_list = vl, h = 4, y_vars = c("Y.1", "Y.2"))
#'
#' ## time-varying association for marker 1 (B-spline basis), constant for marker 2
#' spl    <- list(df = 3, degree = 1, knots = NULL, Bknots = NULL)
#' fit_tv <- LJM(td, y_vars = c("Y.1", "Y.2"), s = 5, h = 4, base_terms = "drug",
#'               ker = "epanechnikov", Bs = list(spl, list(NULL)))
#' cb <- confBands.LJM(fit_tv, K = 1)
#' str(cb, max.level = 1)
#' }
#'
#' @seealso \code{\link{LJM}} (fit with \code{Bs} for a time-varying model)
#'
#' @export
confBands.LJM <- function(object, K = 1, alpha = 0.05){

  # New LJM() returns a flat object; the legacy JEL1/JEL2 nested it under $fitJEL.
  fit.JEL <- if (!is.null(object$fitJEL)) object$fitJEL else object
  fit.COX <- object$fitCOX
  B_conf_lst <- fit.JEL$dataMat$B
  if(inherits(B_conf_lst, "bs")){
    B_conf = B_conf_lst
    single = TRUE
  }else{
    B_conf <- fit.JEL$dataMat$B[[K]]
    single = FALSE
  }
  if(is.null(B_conf)){
    stop("B spline configuration is required. \n")
  }else{
    df <- ncol(B_conf)
    B <- predict(B_conf, fit.COX$y[,2])

    eta = fit.JEL$coefficients$eta

    if(single){
      slope_n = paste0("b_1")
      intercept_n = paste0("b_0")
    }else{
      slope_n = paste0("b", K, "_1")
      intercept_n = paste0("b", K, "_0")
    }
    # Slope -------------------------------------------------------------------
    ind_eta <- grep(slope_n, names(eta))
    eta_b1 = eta[ind_eta]
    eta_b1B = B%*%eta_b1
    ind_rn <- grep(slope_n, rownames(fit.JEL$Vcov))
    # if eigen value is negative then set to 0
    e.dcom <- eigen(fit.JEL$Vcov)
    e.dcom$values[e.dcom$values < 0] <- 0
    new_Vmat <- e.dcom$vectors %*% diag(e.dcom$values) %*% t(e.dcom$vectors)
    Vmat <- new_Vmat[ind_rn, ind_rn]
    var_b1 <- diag(B%*%Vmat%*%t(B))
    se_b1 <- sqrt(var_b1)
    level = 1-alpha
    conf_B <- cbind(eta_b1B, eta_b1B - sqrt(qchisq(level, df))*se_b1, eta_b1B + sqrt(qchisq(level, df))*se_b1)
    colnames(conf_B) <- c("est", "lower", "upper")

    slope_dat = data.frame(tau = fit.COX$y[,2], conf_B)
    # arange data by tau
    slope_dat = slope_dat[order(slope_dat$tau),]
    rownames(slope_dat) <- NULL

    # Intercept ---------------------------------------------------------------

    ind_eta <- grep(intercept_n, names(eta))
    eta_b0 = eta[ind_eta]
    eta_b0B = B%*%eta_b0
    ind_rn <- grep(intercept_n, rownames(fit.JEL$Vcov))
    # if eigen value is negative then set to 0
    e.dcom <- eigen(fit.JEL$Vcov)
    e.dcom$values[e.dcom$values < 0] <- 0
    new_Vmat <- e.dcom$vectors %*% diag(e.dcom$values) %*% t(e.dcom$vectors)
    Vmat <- new_Vmat[ind_rn, ind_rn]
    var_b0 <- diag(B%*%Vmat%*%t(B))
    #take 1/2 power of var_b1
    se_b0 <- sqrt(var_b0)
    conf_B <- cbind(eta_b0B, eta_b0B - sqrt(qchisq(level, df))*se_b0, eta_b0B + sqrt(qchisq(level, df))*se_b0)
    colnames(conf_B) <- c("est", "lower", "upper")


    intercept_dat = data.frame(tau = fit.COX$y[,2], conf_B)
    # arange data by tau
    intercept_dat = intercept_dat[order(intercept_dat$tau),]
    rownames(intercept_dat) <- NULL

  }
  return(list(intercept = intercept_dat, slope = slope_dat))
}
