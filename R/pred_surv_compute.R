pred_surv_compute <- function(fit.JEL, fit.COX, Cox_dat.tr, Cox_dat.te, BLUP_new, rho, CI, MC, alpha){

  #### For test
  # fit.JEL = JEL_models[[1]]$JELmodels[[1]]

  lamb <- fit.JEL$coefficients$lamb
  surv_t_lamb <- lamb$time

  eta <- fit.JEL$coefficients$eta %>% as.matrix()

  phi <- fit.JEL$coefficients$phi
  if(phi %>% length == 0) phi <- 0

  # BLUP_new
  ID1_surv <- as.vector(Cox_dat.te$id) # Subjects
  uniqueID_surv <- !duplicated(ID1_surv) #
  tempID_surv <- which(uniqueID_surv)
  tempID_surv <- c(tempID_surv, length(ID1_surv) + 1)
  ni_surv <- diff(tempID_surv)
  ID_surv <- rep(1:sum(uniqueID_surv), times = ni_surv) # ID variable

  B <- fit.JEL$dataMat$B
  if(inherits(B, "bs")){
    B <- list(B)
  }


  formSurv <- formula(fit.COX)
  TermsSurv <- fit.COX$terms # terms import
  mfSurv <- model.frame(TermsSurv, Cox_dat.te)
  W_new <- as.matrix(model.matrix(formSurv, mfSurv)[,-1]) # remove intercept part

  n_surv = length(ID_surv)
  bashaz = lamb$bashaz


  # rho <- 0
  trans_fun <- function(x){
    if(rho == 0){
      x
    }else{
      log(1 + rho*x)/rho
    }
  }

  if(CI == FALSE){

    if(is.list(B)){

      eta.inds <- sapply(B, function(x) 2*ncol(x))
      BBc = do.call(cbind, B)
      # Duplicate each column and place duplicates next to originals
      BBcd <- do.call(cbind, lapply(seq_len(ncol(BBc)), function(i) {
        cbind(BBc[, i], BBc[, i])
      }))

      etamat = matrix(eta, nrow = nrow(BBcd), ncol = ncol(BBcd), byrow = TRUE)
      BBeta = BBcd * etamat

      # Function to add odd and even columns
      add_odd_even_columns <- function(matrix_data) {
        odd_indices <- seq(1, ncol(matrix_data), 2)   # Indices for odd columns
        even_indices <- seq(2, ncol(matrix_data), 2) # Indices for even columns

        # Sum odd and even columns
        odd_sum <- rowSums(matrix_data[, odd_indices, drop = FALSE])
        even_sum <- rowSums(matrix_data[, even_indices, drop = FALSE])

        # Combine results into a new matrix
        result_matrix <- cbind(odd_sum, even_sum)

        return(result_matrix)
      }

      start = 1
      end = 0
      etaBBs_mat = matrix(0, nrow = nrow(BBeta), ncol = 2*length(eta.inds))
      for(k in 1:length(eta.inds)){
        eta1 = eta.inds[k]
        end = end + eta1
        etaBB = BBeta[, start:end]
        etaBBs = add_odd_even_columns(etaBB)
        etaBBs_mat[, (2*k-1):(2*k)] = etaBBs
        start = end + 1
      }

      BLUP_eta <- lapply(1:n_surv, function(i) etaBBs_mat%*%BLUP_new$mean[,i])

      if(ncol(W_new) == 0){
        exp_h <- lapply(1:n_surv, function(i) exp(BLUP_eta[[i]]))
      }else{
        exp_h <- lapply(1:n_surv, function(i) exp(as.numeric(W_new[i,]%*%phi) +  BLUP_eta[[i]]))
      }

      cum_haz <- lapply(1:n_surv, function(i) apply(exp_h[[i]] * bashaz, 2, cumsum))

      surv_f <- lapply(1:n_surv, function(i) exp(-trans_fun(cum_haz[[i]])))

      surv_tab <- lapply(1:n_surv, function(i) data.frame(tau = surv_t_lamb,
                                                          pred_surv = surv_f[[i]])
      )

    }else{
      ## upper and lower bound of eta*b_i
      # chi_val <- qchisq((1-alpha), nrow(eta))

      if(is.null(ncol(BLUP_new$mean))){
        BLUP_eta <- lapply(1:length(BLUP_new$mean), function(i) eta*BLUP_new$mean[i])
      }else{
        BLUP_eta <- lapply(1:ncol(BLUP_new$mean), function(i) t(eta)%*%BLUP_new$mean[,i])
      }

      if(ncol(W_new) == 0){
        exp_h <- lapply(1:n_surv, function(i) exp(BLUP_eta[[i]]))
      }else{
        exp_h <- lapply(1:n_surv, function(i) exp(as.numeric(W_new[i,]%*%phi) +  BLUP_eta[[i]]))
      }
      cum_haz <- lapply(1:n_surv, function(i) apply(exp_h[[i]] %o% bashaz, 1, cumsum))

      surv_f <- lapply(1:n_surv, function(i) exp(-trans_fun(cum_haz[[i]])))

      surv_tab <- lapply(1:n_surv, function(i) data.frame(tau = surv_t_lamb,
                                                          pred_surv = surv_f[[i]]
                                                          # upper = upper_surv[[i]],
                                                          # lower = lower_surv[[i]]
      ))
    }
    names(surv_tab) = unique(ID1_surv)
  }else{
    if(fit.JEL$Vcov %>% is.null){
      stop("Variance-covariance matrix of eta is not available. Please set CI = FALSE.")
    }
    # index check
    Vcov = fit.JEL$Vcov
    etaVcov_ind = grepl("eta:eta", rownames(Vcov))
    phiVcov_ind = grepl("phi:", rownames(Vcov))

    eta_MC = matrix(0, nrow = length(eta), ncol = MC)
    for(i in 1:MC){
      eta_MC[,i] = MASS::mvrnorm(1, eta, Vcov[etaVcov_ind, etaVcov_ind])
    }

    if(ncol(W_new) == 0){
      phi_MC = matrix(0, nrow = 1, ncol = MC)
    }else{
      phi_MC = matrix(0, nrow = length(phi), ncol = MC)
      for(i in 1:MC){
        phi_MC[,i] = MASS::mvrnorm(1, phi, Vcov[phiVcov_ind, phiVcov_ind])
      }
    }

    phi_MC = matrix(0, nrow = length(phi), ncol = MC)
    for(i in 1:MC){
      phi_MC[,i] = MASS::mvrnorm(1, phi, Vcov[phiVcov_ind, phiVcov_ind])
    }

    if(is.list(B)){
      surv_tab_list = list()
      for(mm in 1:MC){
        eta = eta_MC[,mm]
        phi = phi_MC[,mm]
        eta.inds <- sapply(B, function(x) 2*ncol(x))
        BBc = do.call(cbind, B)
        # Duplicate each column and place duplicates next to originals
        BBcd <- do.call(cbind, lapply(seq_len(ncol(BBc)), function(i) {
          cbind(BBc[, i], BBc[, i])
        }))

        etamat = matrix(eta, nrow = nrow(BBcd), ncol = ncol(BBcd), byrow = TRUE)
        BBeta = BBcd * etamat

        # Function to add odd and even columns
        add_odd_even_columns <- function(matrix_data) {
          odd_indices <- seq(1, ncol(matrix_data), 2)   # Indices for odd columns
          even_indices <- seq(2, ncol(matrix_data), 2) # Indices for even columns

          # Sum odd and even columns
          odd_sum <- rowSums(matrix_data[, odd_indices, drop = FALSE])
          even_sum <- rowSums(matrix_data[, even_indices, drop = FALSE])

          # Combine results into a new matrix
          result_matrix <- cbind(odd_sum, even_sum)

          return(result_matrix)
        }

        start = 1
        end = 0
        etaBBs_mat = matrix(0, nrow = nrow(BBeta), ncol = 2*length(eta.inds))
        for(k in 1:length(eta.inds)){
          eta1 = eta.inds[k]
          end = end + eta1
          etaBB = BBeta[, start:end]
          etaBBs = add_odd_even_columns(etaBB)
          etaBBs_mat[, (2*k-1):(2*k)] = etaBBs
          start = end + 1
        }

        BLUP_eta <- lapply(1:n_surv, function(i) etaBBs_mat%*%BLUP_new$mean[,i])

        if(ncol(W_new) == 0){
          exp_h <- lapply(1:n_surv, function(i) exp(BLUP_eta[[i]]))
        }else{
          exp_h <- lapply(1:n_surv, function(i) exp(as.numeric(W_new[i,]%*%phi) +  BLUP_eta[[i]]))
        }

        cum_haz <- lapply(1:n_surv, function(i) apply(exp_h[[i]] * bashaz, 2, cumsum))

        surv_f <- lapply(1:n_surv, function(i) exp(-trans_fun(cum_haz[[i]])))

        surv_tab <- lapply(1:n_surv, function(i) data.frame(tau = surv_t_lamb,
                                                            pred_surv = surv_f[[i]])
        )
        surv_tab_list[[mm]] = surv_tab
      }
      pred_surv_list = lapply(surv_tab_list, function(x) do.call('cbind', x))
      pred_surv_mat = lapply(1:ncol(pred_surv_list[[1]]), function(i) do.call('cbind', lapply(pred_surv_list, function(x) x[,i])))
      pred_surv_mean = lapply(pred_surv_mat, function(x) apply(x, 1, mean))
      pred_surv_CI =lapply(pred_surv_mat, function(x)  apply(x, 1, quantile, probs = c(alpha/2, 1-alpha/2)))
      surv_tab = lapply(1:length(pred_surv_mean), function(i) data.frame(tau = surv_t_lamb, pred_surv = pred_surv_mean[[i]], lower = pred_surv_CI[[i]][1,], upper = pred_surv_CI[[i]][2,]))
    }else{

      surv_tab_list = list()
      for(mm in 1:MC){
        eta = eta_MC[,mm]
        phi = phi_MC[,mm]
        if(is.null(ncol(BLUP_new$mean))){
          BLUP_eta <- lapply(1:length(BLUP_new$mean), function(i) eta*BLUP_new$mean[i])
        }else{
          BLUP_eta <- lapply(1:ncol(BLUP_new$mean), function(i) t(eta)%*%BLUP_new$mean[,i])
        }
        if(ncol(W_new) == 0){
          exp_h <- lapply(1:n_surv, function(i) exp(BLUP_eta[[i]]))
        }else{
          exp_h <- lapply(1:n_surv, function(i) exp(as.numeric(W_new[i,]%*%phi) +  BLUP_eta[[i]]))
        }
        cum_haz <- lapply(1:n_surv, function(i) apply(exp_h[[i]] %o% bashaz, 1, cumsum))

        surv_f <- lapply(1:n_surv, function(i) exp(-trans_fun(cum_haz[[i]])))

        surv_tab <- lapply(1:n_surv, function(i) data.frame(pred_surv = surv_f[[i]]))
        surv_tab_list[[mm]] = surv_tab
      }
      pred_surv_list = lapply(surv_tab_list, function(x) do.call('cbind', x))
      pred_surv_mat = lapply(1:ncol(pred_surv_list[[1]]), function(i) do.call('cbind', lapply(pred_surv_list, function(x) x[,i])))
      pred_surv_mean = lapply(pred_surv_mat, function(x) apply(x, 1, mean))
      pred_surv_CI =lapply(pred_surv_mat, function(x)  apply(x, 1, quantile, probs = c(alpha/2, 1-alpha/2)))
      surv_tab = lapply(1:length(pred_surv_mean), function(i) data.frame(tau = surv_t_lamb, pred_surv = pred_surv_mean[[i]], lower = pred_surv_CI[[i]][1,], upper = pred_surv_CI[[i]][2,]))
    }
    names(surv_tab) = unique(ID1_surv)
  }
  return(surv_tab)
}
