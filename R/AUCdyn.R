#' Dynamic Area Under the Curve (AUC) Calculation
#'
#' This function computes the dynamic area under the curve (AUC) for survival predictions at a specified landmark time and prediction horizon.
#'
#' @param pred_surv_result The result from the \code{predict.LJM} function. It can be a \code{data.frame} or a \code{list}:
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
#' The function calculates the dynamic AUC by comparing predicted survival probabilities from \code{predict.LJM} against actual survival outcomes at the specified landmark and prediction horizon. It uses Kaplan-Meier estimates to adjust survival probabilities for specific cases.
#'
#' @return A numeric value representing the dynamic AUC.
#'
#' @examples
#' data("pbc2", package = "LJM")
#' d <- pbc2; d$id <- as.numeric(d$id)
#'
#' ## a simple illustrative predictor for the subjects at risk at s = 5:
#' ## survival probability decreasing in the last observed log bilirubin
#' last <- d[d$year <= 5 & d$years > 5, ]
#' last <- last[!duplicated(last$id, fromLast = TRUE), ]
#' pr <- data.frame(id = last$id, tau = 2,
#'                  pred_surv = exp(-0.1 * last$serBilir))
#'
#' ev <- list(id = "id", EvTime = "years", event = "status2")
#' AUCdyn(pr, data = d, landmarks = 5, tau = 2, var_list = ev)
#'
#' @seealso \code{\link{PEdyn}}, \code{\link[=predict.LJM]{predict}}, \code{\link{LJM}}
#'
#'
#' @export
AUCdyn <- function(pred_surv_result, data, landmarks, tau,
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

  # KM_est compute: one row per subject. `data` is usually in long format (one
  # row per visit); fitting the KM on it would count each subject once per visit.
  subj_dat <- data[!duplicated(data[[var_list$id]]), , drop = FALSE]
  KM_est <- survfit(Surv(subj_dat[[var_list$EvTime]], subj_dat[[var_list$event]]) ~ 1)

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

  pairs <- combn(1:nrow(filter_test.id), 2)
  new_pairs <- rbind(pairs[2, ], pairs[1, ])
  pairs <- cbind(pairs, new_pairs)

  T_surv <- filter_test.id[[var_list$EvTime]]
  event <- filter_test.id[[var_list$event]]

  Ti <- T_surv[pairs[1, ]]
  Tj <- T_surv[pairs[2, ]]

  di <- event[pairs[1, ]]
  dj <- event[pairs[2, ]]


  pi.u.t.i <- pred_result[pairs[1, ], -c(1,2)] # remove id and tau
  pi.u.t.j <- pred_result[pairs[2, ], -c(1,2)] # remove id and tau

  concord <- (pi.u.t.i < pi.u.t.j) # model's concord

  tau0 = tau+landmarks

  ind1 <- (Ti <= tau0 & di == 1) & Tj > tau0
  ind2 <- (Ti <= tau0 & di == 0) & Tj > tau0
  ind3 <- (Ti <= tau0 & di == 1) & (Tj <= tau0 & dj == 0)
  ind4 <- (Ti <= tau0 & di == 0) & (Tj <= tau0 & dj == 0)

  ind <- ind1 | ind2 | ind3 | ind4

  KM_surv <- data.frame(time = KM_est$time, surv = KM_est$surv) %>% dplyr::filter(time >= landmarks)
  num_ind <- max(which(KM_surv$time <= tau0))

  if(any(ind2, na.rm = T)){
    Ti_ind2_denum <- sapply(Ti[ind2], function(x)max(which(KM_surv$time <= x))) #pi(T_i)
    auc2_v <- sapply(1:length(Ti[ind2]), function(x) 1-(KM_surv$surv[num_ind])/(KM_surv$surv[Ti_ind2_denum[x]])) # 1-pi(tau)/pi(T_i)
    ind[ind2] <- ind[ind2]*auc2_v
  }

  if(any(ind3, na.rm = T)){
    Tj_ind3_denum <- sapply(Tj[ind3], function(x)max(which(KM_surv$time <= x))) #pi(T_j)
    auc3_v <- sapply(1:length(Tj[ind3]), function(x) (KM_surv$surv[num_ind])/(KM_surv$surv[Tj_ind3_denum[x]])) # pi(tau)/pi(T_j)
    ind[ind3] <- ind[ind3]*auc3_v
  }

  if(any(ind4, na.rm = T)){
    Ti_ind4_denum <- sapply(Ti[ind4], function(x)max(which(KM_surv$time <= x))) #pi(T_i)
    Tj_ind4_denum <- sapply(Tj[ind4], function(x)max(which(KM_surv$time <= x))) #pi(T_j)
    auc4_vi <- sapply(1:length(Ti[ind4]), function(x) 1-(KM_surv$surv[num_ind])/(KM_surv$surv[Ti_ind4_denum[x]])) # 1-pi(tau)/pi(T_i)
    auc4_vj <- sapply(1:length(Tj[ind4]), function(x) (KM_surv$surv[num_ind])/(KM_surv$surv[Tj_ind4_denum[x]])) # pi(tau)/pi(T_j)
    auc4_v <- auc4_vi*auc4_vj
    ind[ind4] <- ind[ind4]*auc4_v
  }

  auc <- sum(c(ind) * concord, na.rm = T)/sum(ind, na.rm = T)

  return(auc)
}
