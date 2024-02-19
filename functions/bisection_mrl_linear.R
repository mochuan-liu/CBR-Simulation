# Bisection for mrl
bisection_mrl_linear <- function(H_list, A_list, P_list, P, Reward, Risk, tau, 
                           init_type = 'model', mu2, lambda_vec, epsilon, eta = eta){
  
  suppressMessages(require(glmnet))
  T <- length(A_list)
  N <- length(A_list[[1]])
  
  # Initialize gamma_max and gamma_min
  gamma_min <- 1
  gamma_max <- 0
  R_max <- ((1-gamma_max)*Reward - gamma_max*Risk)
  R_min <- ((1-gamma_min)*Reward - gamma_min*Risk)
  
  # Estimate initial end point (min)
  if(init_type == 'model'){
    init_min <- get_reg_init_linear(H_list, A_list, Y = R_min)
    init_max <- get_reg_init_linear(H_list, A_list, Y = R_max)
  }
  
  R_max <- subtract_mean(R_max, H_list[[1]])
  R_min <- subtract_mean(R_min, H_list[[1]])
  
  res_min <- try(single_dc_linear(H_list = H_list, A_list = A_list, R = R_min/P, 
                                  beta_start = init_min$beta, beta0_start = init_min$beta0, 
                                  mu2 = mu2, lambda_vec = lambda_vec), silent = TRUE)
  
  # Estimate initial end point (max)
  res_max <- try(single_dc_linear(H_list = H_list, A_list = A_list, R = R_max/P, 
                                  beta_start = init_max$beta, beta0_start = init_max$beta0, 
                                  mu2 = mu2, lambda_vec = lambda_vec), silent = TRUE)
  
  if(!is.list(res_min)|!is.list(res_max)){
    stop('Initial estimation failed!')
  }
  
  Risk_min <- compute_risk(Risk, P, A_list, H_list, res_min, eta = eta)
  Risk_max <- compute_risk(Risk, P, A_list, H_list, res_max, eta = eta)
  
  # Break when risk constraint is too small or too large
  if(tau > Risk_max){
    warning('Risk constraint is too large, return unconstrained!')
    return(list(res = res_max, gamma = gamma_max))
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
    
    if(init_type == 'model'){
      init <- get_reg_init_linear(H_list, A_list, Y = R_now)
    }
    
    R_now <- subtract_mean(R_now, H_list[[1]])
    
    res <- try(single_dc_linear(H_list = H_list, A_list = A_list, R = R_now/P, 
                                beta_start = init$beta, beta0_start = init$beta0, 
                                mu2 = mu2, lambda_vec = lambda_vec), silent = TRUE)
    
    if(is.list(res)){
      Risk_now <- compute_risk(Risk, P, A_list, H_list, res, eta = eta)
    }else{
      stop('DC algoritm failed to converge!')
    }
    
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
  return(list(res = res_min, gamma = gamma_min))
}


# Bisection search with CV
br2_mrl_linear <- function(H_list, A_list, P_list, P, Reward, Risk, tau, 
                           mu2_list = 10^(-8), lambda_list = list(c(0.01,0.1), c(0.01,0.1)),
                           init_type = 'model', cv.folds = 2, epsilon = 0.001, eta = 0.0001){

  T <- length(A_list)
  N <- length(A_list[[1]])


  # CV training/testing data
  flods.list <- caret::createFolds(1:N, cv.folds)
  dat.list <- list()

  for(k.cv in 1:cv.folds){
    temp.list <- list(H_list_training = list(),  A_list_training = list(), P_list_training = list(), P_training = c(),
                      H_list_testing  = list(),  A_list_testing  = list(), P_list_testing  = list(), P_testing  = c(),
                      Risk_training = c(), Reward_training = c(),
                      Risk_testing  = c(), Reward_testing  = c())

    index.training <- flods.list[[k.cv]]
    index.testing  <- setdiff(1:N, index.training)

    temp.list$Risk_training <- Risk[index.training]
    temp.list$Risk_testing  <- Risk[index.testing]
    temp.list$Reward_training  <- Reward[index.training]
    temp.list$Reward_testing   <- Reward[index.testing]
    
    temp.list$P_training <- P[index.training]
    temp.list$P_testing  <- P[index.testing]
  
    for(t in 1:T){
      temp.list$H_list_training[[t]] <- H_list[[t]][index.training, ]
      temp.list$A_list_training[[t]] <- A_list[[t]][index.training]
      temp.list$P_list_training[[t]] <- P_list[[t]][index.training]

      temp.list$H_list_testing[[t]] <- H_list[[t]][index.testing, ]
      temp.list$A_list_testing[[t]] <- A_list[[t]][index.testing]
      temp.list$P_list_testing[[t]] <- P_list[[t]][index.testing]
    }
    dat.list[[k.cv]] <- temp.list
  }

  # Get tuning parameter list
  para_list <- c(lambda_list, list(mu2_list))
  dat_para  <- do.call(expand.grid, para_list) 
  names(dat_para)[T+1] <- 'mu2'

  # CV
  value_training <- rep(0, dim(dat_para)[1])
  value_testing  <- rep(0, dim(dat_para)[1])
  risk_training  <- rep(0, dim(dat_para)[1])
  risk_testing   <- rep(0, dim(dat_para)[1])
  for(i in 1:dim(dat_para)[1]){

    mu2_now <- dat_para$mu2[i] %>% as.numeric()
    lambda_vec_now <- dat_para[i,1:T] %>% unlist()

    for(k.cv in 1:cv.folds){
      dat <- dat.list[[k.cv]]

      H_training_list <- dat$H_list_training
      A_training_list <- dat$A_list_training
      P_training_list <- dat$P_list_training
      P_training      <- dat$P_training
      Reward_training <- dat$Reward_training
      Risk_training   <- dat$Risk_training

      H_testing_list <- dat$H_list_testing
      A_testing_list <- dat$A_list_testing
      P_testing_list <- dat$P_list_testing
      P_testing      <- dat$P_testing
      Reward_testing <- dat$Reward_testing
      Risk_testing   <- dat$Risk_testing
      
      # Bisection search
      o_search <- try(bisection_mrl_linear(H_list = H_training_list, A_list = A_training_list, P_list = P_training_list, P = P_training, 
                                     Risk = Risk_training, Reward = Reward_training, tau = tau,
                                     mu2 = mu2_now, lambda_vec = lambda_vec_now, 
                                     init_type = init_type, epsilon = epsilon, eta = eta), silent = TRUE)
      if(!is.list(o_search)){
        warning('DC encounters error!')
        value_training[i] <- NA
        value_testing[i]  <- NA
        
        risk_training[i] <- NA
        risk_testing[i]  <- NA
        break
      }
      
      # Optimal estimated training/testing 
      opt_training <- est_opt(o_search$res, H_testing = H_training_list)
      opt_testing  <- est_opt(o_search$res, H_testing = H_training_list)
      
      if(k.cv == 1){
        value_training[i] <-  compute_value(Reward_training, P_training, A_training_list, opt_training)
        value_testing[i]  <-  compute_value(Reward_testing,  P_testing,  A_testing_list,  opt_testing)
        
        risk_training[i] <- compute_value(Risk_training, P_training, A_training_list, opt_training)
        risk_testing[i]  <- compute_value(Risk_testing,  P_testing,  A_testing_list,  opt_testing)
        
      }else{
        value_training[i] <-  value_training[i] + compute_value(Reward_training, P_training, A_training_list, opt_training)
        value_testing[i]  <-  value_testing[i]  + compute_value(Reward_testing,  P_testing,  A_testing_list,  opt_testing)
        
        risk_training[i] <-  risk_training[i] + compute_value(Risk_training, P_training, A_training_list, opt_training)
        risk_testing[i]  <-  risk_testing[i]  + compute_value(Risk_testing,  P_testing,  A_testing_list,  opt_testing)
      }
    }
   
    cat('Current lambda: ', paste(round(lambda_vec_now, 5), collapse = ' , '), ' Current mu2: ', round(mu2_now, 6), ' Completed!', '\n')
  }

  # Get optimal tuning pair
  if(all(is.na(value_testing))){
    stop('No valid tuning parameter found! All are NA!')
  }
  
  index_opt <- which.max(value_testing)
  mu2_max <- dat_para$mu2[index_opt]
  lambda_vec_max <- dat_para[index_opt, 1:T]
  
  cat('Max lambda: ', paste(round(lambda_vec_max, 5), collapse = ' , '), ' Max mu2: ', round(mu2_max, 6), '\n')

  # Final estimation
  o <- try(bisection_mrl_linear(H_list = H_list, A_list = A_list, P = P, Risk = Risk, Reward = Reward, tau = tau,
                          mu2 = mu2_max, lambda_vec = lambda_vec_max, 
                          init_type = init_type, epsilon = epsilon, eta = eta))
  if(!is.character(o)){
    o$lambda <- lambda_vec_max
    o$mu2 <- mu2_max
  }
  return(o)
}


# Compute cumulative value 
compute_value <- function(V, P, A, A_est){
  T <- length(A)
  for(t in 1:T){
    V <- V*(A[[t]]==A_est[[t]])
  }
  return(mean(V/P))
}

compute_risk <- function(V, P, A, H_testing, out, eta = 0.0001){
  
  T <- length(A)
  N <- length(A[[1]])
  beta <- out$beta
  beta0 <- out$beta0
  
  af_t_mat <- matrix(Inf, N, T)
  
  for(t in 1:T){
    f_t <- cbind(1, H_testing[[t]]) %*% (c(beta0[[t]], beta[[t]]))
    af_t <- f_t*A[[t]]
    af_t <- ifelse(af_t<=0, 0, af_t)
    af_t <- ifelse(af_t>=eta, 1, af_t/eta)
    af_t_mat[,t] <- af_t
  }
  
  V <- V*(apply(af_t_mat, 1, min))
  return(mean(V/P))
}

subtract_mean <- function(Y, H){
  cv.res <- cv.glmnet(x = H, y = Y, family = 'gaussian')
  cv.lambda <- cv.res$lambda.min
  y.model.lm <- glmnet(x = H, y = Y, lambda = cv.lambda, family = 'gaussian') 
  y.hat <- predict(y.model.lm, newx =  H) %>% as.vector()
  Y <- Y - y.hat
  return(Y)
}
