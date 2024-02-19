# single stage Q-learning
qlearning <- function(K_list, A_list, P_list, Y){
  
  # using DTRlearn2
  suppressMessages(library(DTRlearn2))
  T <- length(A_list)
  N <- length(A_list[[1]])

  R_list <- list()
  for(t in 1:T){
    if(t==T){
      R_list[[t]] <- Y
    }else{
      R_list[[t]] <- rep(0, N)
    }
  }

  res <- try(ql(H = K_list, AA = A_list, RR = R_list, pi = P_list, K = T, lasso = TRUE), silent = TRUE)
  if(!is.list(res)){ stop('Q-learning estimation failed!') }

  beta <- list()
  beta0 <- list()

  for(t in 1:T){
    p <- dim(K_list[[t]])[2]
    keyword <- paste('stage', t, sep = '')
    beta0[[t]] <- res[[keyword]]$co[p+2]
    beta[[t]] <-  res[[keyword]]$co[(p+3):(2*p+2)]
  }

  return(list(beta = beta, beta0 = beta0))
}

# Bisection for Q-learning
br2_q <- function(H_list, A_list, P_list, P, Reward, Risk, tau, 
                  epsilon, eta, sigma_list = NULL){
  
  # Get number of stage and sample size
  suppressMessages(require(glmnet))
  suppressMessages(require(DTRlearn2))
  T <- length(A_list)
  N <- length(A_list[[1]])
  
  # Get kernel bandwidth
  if(kernel == 'gaussian' & is.null(sigma_list)){
    sigma_list <- list()
    for(t in 1:T){
      sigma_list[[t]] <- rbf_bandwidth(H_list[[t]], A_list[[t]])
    }
  }
  
  K_list <- K_matrix_list(H_list = H_list, kernel = kernel, sigma_list = sigma_list)
  
  # Initialize gamma_max and gamma_min
  gamma_min <- 1
  gamma_max <- 0
  
  R_max <- ((1-gamma_max)*Reward - gamma_max*Risk)
  R_min <- ((1-gamma_min)*Reward - gamma_min*Risk)
  
  # Estimate initial end point (max/min)
  res_max <- try(qlearning(K_list = K_list, A_list = A_list, P_list = P_list, Y = R_max), silent = TRUE)
  res_min <- try(qlearning(K_list = K_list, A_list = A_list, P_list = P_list, Y = R_min), silent = TRUE)
  
  if(!is.list(res_min)|!is.list(res_max)){
    stop('Initial estimation failed!')
  }
  
  opt_min <- est_opt(res_min, K_list = K_list)
  opt_max <- est_opt(res_max, K_list = K_list)
  
  Risk_min <- compute_risk(Risk, P, A_list, K_list, res_min, eta)
  Risk_max <- compute_risk(Risk, P, A_list, K_list, res_max, eta)
  
  # Break when risk constraint is too small or too large
  if(tau > Risk_max){
    warning('Risk constraint is too large, return unconstrained!')
    return(list(res = res_max, gamma = gamma_max, sigma_list = sigma_list))
  }
  
  # When risk constraint is too small
  if(tau < Risk_min){
    stop('Risk constraint is too small, please try larger risk constraint!')
  }
  
  # Bisection search 
  n_try <- 1
  while(gamma_min - gamma_max>=epsilon&n_try<=100){
    gamma_now <- (gamma_min + gamma_max)/2
    n_try <- n_try + 1
    
    R_now <- ((1-gamma_now)*Reward - gamma_now*Risk)
    res <- try(qlearning(K_list = K_list, A_list = A_list, P_list = P_list, Y = R_now), silent = TRUE)
    
    if(!is.list(res)){
      stop('Q-Learning failed!')
    }
    
    opt <- est_opt(res, K_list = K_list)
    Risk_now <- compute_risk(Risk, P, A_list, K_list, res, eta)
    
    if(Risk_now>tau){
      Risk_max <- Risk_now
      gamma_max <- gamma_now
      res_max <- res
    }else{
      Risk_min <- Risk_now
      gamma_min <- gamma_now
      res_min <- res
    }
    # cat("Risk_now: ", Risk_now, "  Gamma_min: ", gamma_min, " Risk_min: ", Risk_min,
    #     "  Gamma_max: ", gamma_max, " Risk_max: ", Risk_max, '\n')
  }
  
  # Output final estimated coefficients
  return(list(res = res_min, gamma = gamma_min, sigma_list = sigma_list))
}

