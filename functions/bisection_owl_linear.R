owl_single <- function(H, A, P, Y, lambda){

  # load solve.QP
  suppressMessages(require(quadprog))

  H <- as.matrix(H)
  N <- length(A)
  Y.origin <- Y
  A.origin <- A

  # standardize
  cv.res <- cv.glmnet(x = H, y = Y, family = 'gaussian')
  cv.lambda <- cv.res$lambda.min
  m.model <- glmnet(x = H, y = Y, lambda = cv.lambda, family = 'gaussian')
  m.hat <- predict(m.model, newx = H)

  A <- A*sign(Y-m.hat) %>% as.vector()
  Y <- abs(Y-m.hat)
  if(min(Y)<10^(-3))  Y <- Y + 10^(-3)

  # Solve QP
  K <- diag(A)%*%H%*%t(H)%*%diag(A) + diag(N)*10^(-7)
  d <- rep(1, N)
  W <- rbind(matrix(A, nrow = 1), diag(N), -diag(N))
  b <- c(0, rep(0, N), -(lambda)^(-1)*Y/P)

  s <- solve.QP(Dmat = K, dvec = d, Amat = t(W), bvec = b, meq = 1)

  alpha <- diag(A)%*%matrix(s$solution, ncol = 1)
  beta <- t(H)%*%alpha

  # search beta0
  beta0_grid <- seq(-20,20,0.01)
  value <- rep(NA, length(beta0_grid))
  K_alpha <- H%*%t(H)%*%alpha
  for(i in 1:length(beta0_grid)){
    l <- 1 - A*(K_alpha+beta0_grid[i])
    value[i] <- mean(ifelse(l>0, l, 0)*Y/P)
  }
  beta0 <- beta0_grid[which.min(value)]

  return(list(beta = beta %>% matrix(ncol = 1), beta0 = beta0))
}

weight_least_square <- function(H, Y, weight){
  cv.res <- cv.glmnet(x = H, y = Y, weights = weight, family = 'gaussian')
  cv.lambda <- cv.res$lambda.min
  m.model <- glmnet(x = H, y = Y, weights = weight, lambda = cv.lambda, family = 'gaussian') 
  return(m.model)
}

est_opt_owl_single <- function(out, H_testing = NULL){
  beta <- out$beta %>% matrix(ncol = 1)
  beta0 <- out$beta0
  H_testing <- as.matrix(H_testing)
  return(sign(H_testing%*%beta+beta0))
}

# AOWL
aowl <- function(H_list, A_list, P_list, Y, lambda_vec){
  
  T <- length(A_list)
  A_est_list <- list()
  N <- length(Y)
  
  model_list <- list()
  for(t in 1:T){
    model_list[[t]] <- list()
  }
  out <- list(beta = list(), beta0 = list())
  
  for(t in T:1){
    if(t == T){
      # Final stage estimation
      o <- owl_single(H = H_list[[T]], A = A_list[[T]], P = P_list[[T]], Y = Y, lambda = lambda_vec[T])
      A_est_list[[T]] <- est_opt_owl_single(out = o, H_testing = H_list[[T]])
      
      # Get imputation model
      weight.now <- ((A_est_list[[T]]==A_list[[T]])/P_list[[T]])*(1-P_list[[T]])/P_list[[T]]
      M.model <- weight_least_square(H_list[[T]], Y, weight.now)
      model_list[[T]][[T]] <- M.model
      
      Y_now <- Y
      weight.1 <- (A_est_list[[T]]==A_list[[T]])/P_list[[T]]
      out$beta[[T]]   <- o$beta
      out$beta0[[T]] <- o$beta0
    }else{
      Y_now <- (Y_now*(A_est_list[[t+1]]==A_list[[t+1]]))/P_list[[t+1]]
      for(s in (t+1):T){
        if(s==t+1){
          aug.weight <- rep(1, N)
          M.model <- model_list[[t+1]][[s]]
          M.est <- predict(M.model, newx = H_list[[s]])
          
          Y_aug <- ((A_est_list[[s]]==A_list[[s]])/P_list[[s]]-1)*M.est
        }else{
          aug.weight <- aug.weight*(A_est_list[[s-1]]==A_list[[s-1]])/P_list[[s-1]]
          M.model <- model_list[[t+1]][[s]]
          M.est <- predict(M.model, newx = H_list[[s]])
          
          Y_aug <- Y_aug + aug.weight*((A_est_list[[s]]==A_list[[1]])/P_list[[s]]-1)*M.est
        }
      }
      
      o <- owl_single(H = H_list[[t]], A = A_list[[t]], P = P_list[[t]], Y = Y_now - Y_aug, lambda = lambda_vec[t])
      A_est_list[[t]] <- est_opt_owl_single(out = o, H_testing = H_list[[t]])
      weight.1 <- weight.1*(A_est_list[[t]]==A_list[[t]])/P_list[[t]]
      
      weight.2 <- rep(1, N)
      for(j in t:T){
        weight.2 <- weight.2*P_list[[j]]
        weight.now <-  weight.1*(1-P_list[[j]])/weight.2
        M.model <- weight_least_square(H_list[[j]], Y, weight.now)
        model_list[[t]][[j]] <- M.model
      }
      out$beta[[t]]   <- o$beta
      out$beta0[[t]] <- o$beta0
    }
  }
  
  return(out)
}

# OWL
bowl <- function(H_list, A_list, P_list, Y, lambda_vec){
  
  T <- length(A_list)
  A_est_list <- list()
  
  out <- list(beta = list(), beta0 = list())
  
  for(t in T:1){
    if(t == T){
      # Final stage estimation
      o <- owl_single(H = H_list[[T]], A = A_list[[T]], P = P_list[[T]], Y = Y, lambda = lambda_vec[T])
      A_est_list[[T]] <- est_opt_owl_single(out = o, H_testing = H_list[[T]])
      
      Y_now <- Y
      out$beta[[T]]   <- o$beta
      out$beta0[[T]] <- o$beta0
    }else{
      Y_now <- Y_now*(A_est_list[[t+1]]==A_list[[t+1]])/P_list[[t+1]]
      o <- owl_single(H = H_list[[t]], A = A_list[[t]], P = P_list[[t]], Y = Y_now, lambda = lambda_vec[t])
      A_est_list[[t]] <- est_opt_owl_single(out = o, H_testing = H_list[[t]])
      
      out$beta[[t]]   <- o$beta
      out$beta0[[t]] <- o$beta0
    }
  }
  
  return(out)
}

owl2 <- function(H_list, A_list, P_list, Y, lambda_vec, method){
  suppressMessages(require(glmnet))
  
  if(method == 'aowl'){
    return(aowl(H_list, A_list, P_list, Y, lambda_vec))
  }else{
    return(bowl(H_list, A_list, P_list, Y, lambda_vec))
  }
}

# Bisection search
bisection_owl <- function(H_list, A_list, P_list, P, Reward, Risk, tau, 
                          lambda_vec, epsilon, eta, method){
  
  suppressMessages(require(glmnet))
  
  T <- length(A_list)
  N <- length(A_list[[1]])
  
  # Estimate initial end point risk
  gamma_max <- 0
  gamma_min <- 1
  R_max <- ((1-gamma_max)*Reward - gamma_max*Risk)
  R_min <- ((1-gamma_min)*Reward - gamma_min*Risk)
  
  res_min <- try(owl2(H_list = H_list, A_list = A_list, P_list = P_list, Y = R_min, lambda_vec = lambda_vec, method = method), silent = TRUE)
  res_max <- try(owl2(H_list = H_list, A_list = A_list, P_list = P_list, Y = R_max, lambda_vec = lambda_vec, method = method), silent = TRUE)
  
  if(!is.list(res_min)|!is.list(res_max)){
    stop('Initial estimation failed!')
  }
  
  Risk_min <- compute_risk(Risk, P, A_list, H_list, res_min, eta)
  Risk_max <- compute_risk(Risk, P, A_list, H_list, res_max, eta)
  
  # When risk constraint is too large
  if(tau > Risk_max){
    warning('Risk constraint is too large, return unconstrained!')
    return(list(res = res_max, gamma = gamma_max))
  }
  
  # When risk constraint is too small
  if(tau < Risk_min){
    stop('Risk constraint is too small, please try larger risk constraint!')
  }
  
  n_try <- 1
  while(gamma_min - gamma_max>=epsilon&n_try <=100){
    n_try <- n_try + 1
    gamma_now <- (gamma_min + gamma_max)/2
    
    R_now <- ((1-gamma_now)*Reward - gamma_now*Risk)
    
    res <- try(owl2(H_list = H_list, A_list = A_list, P_list = P_list, Y = R_now, lambda_vec = lambda_vec, method = method), silent = TRUE)
    
    if(is.list(res)){
      Risk_now <- compute_risk(Risk, P, A_list, H_list, res, eta)
    }else{
      stop('OWL2 algoritm failed to converge!')
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
  
  return(list(res = res_min, gamma = gamma_min))
}

# Bisection selection with CV
br2_owl_linear <- function(H_list, A_list, P_list, P, Reward, Risk, tau, 
                           lambda_list, 
                           cv.folds = 2, epsilon, eta, method){
  
  # Get number of stage and sample size
  T <- length(A_list)
  N <- length(A_list[[1]])
  
  # cv training/testing data
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
  dat_para  <- do.call(expand.grid, lambda_list) 
  
  # Cross validation
  value_training <- rep(NA, dim(dat_para)[1])
  value_testing  <- rep(NA, dim(dat_para)[1])
  risk_training  <- rep(NA, dim(dat_para)[1])
  risk_testing   <- rep(NA, dim(dat_para)[1])
  for(i in 1:dim(dat_para)[1]){
    lambda_vec_now <- dat_para[i,] %>% unlist()
    
    for(k.cv in 1:cv.folds){
      dat <- dat.list[[k.cv]]
      
      H_training_list <- dat$H_list_training
      A_training_list <- dat$A_list_training
      P_training_list <- dat$P_list_training
      P_training <- dat$P_training
      Reward_training <- dat$Reward_training
      Risk_training   <- dat$Risk_training
      
      H_testing_list <- dat$H_list_testing
      A_testing_list <- dat$A_list_testing
      P_testing_list <- dat$P_list_testing
      P_testing <- dat$P_testing
      Reward_testing <- dat$Reward_testing
      Risk_testing   <- dat$Risk_testing
      
      # Bisection search
      o_search <- try(bisection_owl(H_list = H_training_list, A_list = A_training_list, P_list = P_training_list, P = P_training, 
                                    Reward = Reward_training, Risk = Risk_training, tau = tau, 
                                    lambda_vec = lambda_vec_now, epsilon = epsilon, eta = eta, method = method), silent = TRUE)
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
      opt_testing  <- est_opt(o_search$res, H_testing = H_testing_list)
      
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
    
    cat('Current lambda: ', paste(round(lambda_vec_now, 5), collapse = ' , '), ' Completed!', '\n')
  }
  
  # Get optimal tuning pair
  if(all(is.na(value_testing))){
    stop('No valid tuning parameter found! All are NA!')
  }
  
  index_opt <- which.max(value_testing)
  lambda_vec_max <- dat_para[index_opt,] %>% unlist()
  
  cat('Max lambda: ', paste(round(lambda_vec_max, 5), collapse = ' , '), '\n')
  
  # initialize beta_start, beta0_start
  o <- try(bisection_owl(H_list = H_list, A_list = A_list, P_list = P_list, P = P,
                         Reward = Reward, Risk = Risk, tau = tau, 
                         lambda_vec = lambda_vec_max, epsilon = epsilon, eta = eta, method = method))
  if(!is.character(o)){
    o$lambda <- lambda_vec_max
  }
  return(o)
}

# Estimate optimal rule for owl
est_opt_owl_single <- function(out, H_testing = NULL){
  beta <- out$beta %>% matrix(ncol = 1)
  beta0 <- out$beta0
  H_testing <- as.matrix(H_testing)
  return(sign(H_testing%*%beta+beta0))
}


