#' Dynamic Prediction Error (PE) Calculation
#'
#' This function computes the dynamic prediction error (PE) for survival predictions at a specified landmark time and prediction horizon.
#'
#' @param pred_surv_result The result from the \code{predict.JEL} function. It can be a \code{data.frame} or a \code{list}:
#'
#'   - If a \code{data.frame}, it should include columns for \code{id}, \code{tau}, and survival probabilities.
#'
#'   - If a \code{list}, it should include multiple survival prediction outputs.
#' @param data A \code{data.frame} containing the survival data, including subject IDs, event times, and event statuses.
#' @param landmarks A numeric value specifying the landmark time for dynamic prediction.
#' @param tau A numeric value specifying the prediction horizon.
#' @param var_list A \code{list} specifying the column names in the \code{data}:
#'
#'   - \code{id}: Column name for the subject IDs (default: \code{"id"}).
#'
#'   - \code{EvTime}: Column name for the event times (default: \code{"EvTime"}).
#'
#'   - \code{event}: Column name for the event status (default: \code{"event"}).
#'
#' @details
#' The function calculates the dynamic PE by comparing the predicted survival probabilities from \code{predict.JEL} against actual survival outcomes at a specified landmark time and prediction horizon.
#' The prediction error is computed using Kaplan-Meier estimates to adjust for censored cases.
#'
#' @return A numeric value representing the dynamic PE.
#'
#' @examples
#' \dontrun{
#' ## `pr` is the output of predict() on a fitted JEL object
#' ev <- list(id = "id", EvTime = "years", event = "status2")
#' PEdyn(pr, data = mydata, landmarks = 5, tau = 2, var_list = ev)
#' }
#'
#' @seealso \code{\link{AUCdyn}}, \code{\link[=predict.JEL]{predict}}, \code{\link{JEL}}
#'
#'
#' @export
PEdyn <- function(pred_surv_result, data, landmarks, tau,
                   var_list = list(id = "id", EvTime = "EvTime", event = "event")){

#   var_list = list(id = "id", EvTime = "years", event = "status2")
#   data = pbc2
#   landmarks = 3
#   tau = 2
  if(is.data.frame(pred_surv_result)){
    pred_result <- pred_surv_result[, c('id', 'tau', 'pred_surv')]
  }else if(is.list(pred_surv_result)){
    sub_pred <- lapply(pred_surv_result, function(x) x[which.max(x$tau[x$tau <= tau]), c('id', 'tau', 'pred_surv')])
    pred_result <- data.frame(id = as.numeric(names(sub_pred)), do.call('rbind', sub_pred))
  }else {
    stop("pred_surv_dat should be list or data.frame")
  }

  # KM_est compute
  KM_est <- survfit(Surv(data[[var_list$EvTime]], data[[var_list$event]]) ~ 1)

  # test data
  survdata <- data %>% dplyr::select(var_list$id, var_list$EvTime, var_list$event) %>% unique() %>% filter(.[[var_list$EvTime]] >= landmarks)

  # Filter survdata to only include IDs present in pred_result
  filter_test.id <- survdata[which(survdata[[var_list$id]] %in% pred_result$id), ]

  # Check if all IDs in pred_result are present in survdata
  if (!all(pred_result$id %in% filter_test.id[[var_list$id]])) {
    stop("Not all IDs in pred_result are present in survdata. Please ensure both datasets have matching IDs.")
  }

  # Check if the order of IDs matches
  if (!all(pred_result$id == filter_test.id[[var_list$id]])) {
    # Reorder survdata to match the order of IDs in pred_result
    filter_test.id <- filter_test.id[match(pred_result$id, filter_test.id[[var_list$id]]), ]
  }

  T_surv <- filter_test.id[[var_list$EvTime]]
  event <- filter_test.id[[var_list$event]]

  tau0 = tau+landmarks

  n_t <- nrow(filter_test.id) # total survivals
  ind1 <- which(T_surv > tau0) # total # of survivals after tau0
  ind2 <- which(T_surv <= tau0 & event == 1) # total # of events before tau0
  ind3 <- which(T_surv <= tau0 & event == 0) # total # of censored before tau0

  if(any(ind1, na.rm = T)){
    pi_set1 <- pred_result[ind1, -c(1,2)] # remove id and tau
    L1 <- sum((1-pi_set1)^2, na.rm = T)
  }else{
    L1 = 0
  }

  if(any(ind2, na.rm = T)){
    pi_set2 <- pred_result[ind2, -c(1,2)] # remove id and tau
    L2 <- sum((0-pi_set2)^2, na.rm = T)
  }else{
    L2 = 0
  }

  if(any(ind3, na.rm = T)){
    Tsurv_ind3 <- T_surv[ind3]

    # L3 <- matrix(0, nrow = length(Tsurv_ind3), ncol = 1)
    L3 <- numeric(length(Tsurv_ind3))
    KM_surv <- data.frame(time = KM_est$time, surv = KM_est$surv) %>% dplyr::filter(time >= landmarks)
    num_ind <- max(which(KM_surv$time <= tau0))
    pi_set3 <- pred_result[ind3, -c(1,2)] # remove id and tau

    for(j in 1:length(Tsurv_ind3)){
      T_ind <- max(which(KM_surv$time <= Tsurv_ind3[j]))
      weight <- (KM_surv$surv[num_ind])/(KM_surv$surv[T_ind])
      L3[j] <- weight*(1-pi_set3[j])^2 + (1-weight)*(0-pi_set3[j])^2
    }
    L3 <- sum(L3)
  }else{
    L3 = 0
  }

  pe <- n_t^(-1)*(L1+L2+L3)

  return(pe)
}
