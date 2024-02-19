#####################
# Base function for mrl (linear)
#####################

## Global setting
suppressMessages(require(kernlab))
suppressMessages(require(quadprog))
suppressMessages(require(osqp))
suppressMessages(require(caret))
suppressMessages(require(dplyr))
epsilon.max = 0.001
iteration.max = 100

# Compute quadratic matrix
B_matrix_list <- function(K_list, A_list, index_negative){
  
  T <- length(A_list)
  N <- length(A_list[[1]])
  N_negative <- length(index_negative)
  
  B_list <- list()
  one <- matrix(1, nrow = N, ncol = N)
  for(t in 1:T){
    K <- K_list[[t]] + one
    K_sub <- t(K[index_negative,])
    K_subsub <- K[index_negative, index_negative]
    B_t <- rbind(cbind(K, K_sub, -K_sub), cbind(t(K_sub), K_subsub, -K_subsub), cbind(-t(K_sub), -K_subsub, K_subsub))
    A_t <- diag(c(as.vector(A_list[[t]]), rep(1, N_negative), rep(1, N_negative)))
    B_list[[t]] <- A_t%*%B_t%*%A_t
  }
  return(B_list)
}

# Compute smoothed derivatives
c_vec_1_list <- function(H_list, A_list, d, beta_list, beta0_list, mu2){
  
  T <- length(A_list)
  N <- length(A_list[[1]])
  c_list <- list()
  
  c_mat <- matrix(0, N, T)
  f_mat <- matrix(0, N, T)
  for(t in 1:T){
    for(i in 1:N){
      f_mat[i,t] <- A_list[[t]][i]*(H_list[[t]][i,]%*%beta_list[[t]] + beta0_list[[t]])
    }
  }
  for(i in 1:N){
    for(t in 1:T){
      f_vec <- (c(-f_mat[i,], -(1-d[i])) + f_mat[i,t])/mu2 
      f_vec <- f_vec %>% exp()
      c_mat[i,t] <- 1/(sum(f_vec))
    }
  }
  
  for(t in 1:T){ c_list[[t]] <- c_mat[,t] }
  
  return(c_list)
}

c_vec_2_list <- function(H_list, A_list, beta_list, beta0_list, mu2){
  
  T <- length(A_list)
  N <- length(A_list[[1]])
  c_list <- list()
  
  c_mat <- matrix(0, N, T)
  f_mat <- matrix(0, N, T)
  for(t in 1:T){
    for(i in 1:N){
      f_mat[i,t] <- (H_list[[t]][i,]%*%beta_list[[t]] + beta0_list[[t]])
    }
  }
  for(i in 1:N){
    for(t in 1:T){
      c_mat[i,t] <- -1/(1+exp(2*f_mat[i,t]/mu2)) + 1/(1+exp(-2*f_mat[i,t]/mu2))
    }
  }
  
  for(t in 1:T){ c_list[[t]] <- -c_mat[,t] }
  
  return(c_list)
}

# Compute constraint matrix A
A_constraint <- function(N, N_negative, T){
  m <- (N+2*N_negative)
  
  # ui sum
  A_1 <- matrix(0, nrow = N, ncol = m*T)
  for(i in 1:N){
    for(t in 1:T){
      A_1[i, i+(m*(t-1))] <- 1
    }
  }
  
  # li sum
  A_2 <- matrix(0, nrow = N_negative*T, ncol = m*T)
  for(i in 1:N_negative){
    for(t in 1:T){
      A_2[(i-1)*(T)+t, i+(m*(t-1)+N)] <- 1
      A_2[(i-1)*(T)+t, i+(m*(t-1)+N+N_negative)] <- 1
    }
  }
  
  # li sum sum
  A_3 <- matrix(0, nrow = N_negative, ncol = m*T)
  for(i in 1:N_negative){
    for(t in 1:T){
      A_3[i, i+(m*(t-1)+N)] <- 1
      A_3[i, i+(m*(t-1)+N+N_negative)] <- 1
    }
  }
  
  # ui li
  A_4 <- diag(rep(1, m*T))
  
  return(rbind(A_1, A_2, A_3, A_4))
}

# Single DC iteration
single_dc_linear <- function(H_list, A_list, R, beta_start = NULL, beta0_start = NULL, 
                             mu2, lambda_vec){

  N <- length(A_list[[1]])
  T <- length(A_list)
  
  # rescale outcome 
  lambda <- 1
  ratio_vec <- lambda_vec/lambda
  
  # get negative index 
  index_negative <- which(R<0)
  O_negative <- -R[index_negative]/lambda
  O <- abs(R/lambda)
  N_negative <- length(O_negative)
  
  # get subdata 
  H_sub_list <- list()
  A_sub_list <- list()
  for(t in 1:T){
    H_sub_list[[t]] <- H_list[[t]][index_negative, ]
    A_sub_list[[t]] <- A_list[[t]][index_negative]
  }
  one <- matrix(1, nrow = 1, ncol = N)
  one_sub <- matrix(1, nrow = 1, ncol = N_negative)
  one_mat <- matrix(1, nrow = N, ncol = N)
  one_sub_mat <- matrix(1, nrow = N_negative, ncol = N)
  one_subsub_mat <- matrix(1, nrow = N_negative, ncol = N_negative)
  
  # Change O to positive vector
  d <- ifelse(R>0, 1, 0)
  d_vec <- c()
  for(t in 1:T){
    d_vec <- c(d_vec, -d, rep(1, N_negative), rep(1, N_negative))
  }
  
  # compute kernel matrix 
  K_list <- list()
  K_sub_list <- list()
  K_subsub_list <- list()
  for(t in 1:T){ 
    K_list[[t]] <- H_list[[t]]%*%t(H_list[[t]])
    K_sub_list[[t]] <- K_list[[t]][index_negative, ]
    K_subsub_list[[t]] <- K_list[[t]][index_negative, index_negative]
  }
  
  # compute A_mat
  A_mat_list <- list()
  A_submat_list <- list()
  for(t in 1:T){
    A_mat_list[[t]] <- diag(as.vector(A_list[[t]]))
    A_submat_list[[t]] <- A_mat_list[[t]][index_negative, index_negative]
  }
  
  # compute quadratic matrix list
  B_list <- B_matrix_list(K_list = K_list, A_list = A_list, index_negative = index_negative)
  H <- matrix(0, nrow = (N+2*N_negative)*T, ncol = (N+2*N_negative)*T)
  
  for(t in 1:T){
    H[((t-1)*(N+2*N_negative)+1):(t*(N+2*N_negative)),((t-1)*(N+2*N_negative)+1):(t*(N+2*N_negative))] <- B_list[[t]]/(as.numeric(ratio_vec[t]))
  }
  H <- H + diag(10^(-7), dim(H)[1])
  
  # compute W matrix 
  W <- diag(as.vector(O))
  W_negative <- diag(as.vector(O_negative))
  
  # compute constraint matrix 
  A_cons <- A_constraint(N, N_negative, T)
  l_cons <- c(rep(0, length(O)), rep(rep(0, length(O_negative)), each = T), (T-1)*O_negative, rep(0, (N+2*N_negative)*T))
  u_cons <- c(O, rep(O_negative, each = T), T*O_negative, rep(c(O, O_negative, O_negative), times = T))

  s <- 0
  beta_old <- beta_start
  beta0_old <- beta0_start
  while(TRUE){
    
    # compute c_vec
    c_1_list <- c_vec_1_list(H_list = H_list, A_list = A_list, d = d,  beta_list = beta_old, beta0_list = beta0_old, mu2 = mu2)
    c_2_list <- c_vec_2_list(H_list = H_sub_list, A_list = A_sub_list, beta_list = beta_old, beta0_list = beta0_old, mu2 = mu2)
    c_combine <- c()
    for(t in 1:T){
      c_t_1 <- matrix(c_1_list[[t]], nrow = 1)%*%W%*%A_mat_list[[t]]%*%(K_list[[t]] + one_mat)%*%A_mat_list[[t]] +
        matrix(c_2_list[[t]], nrow = 1)%*%W_negative%*%(K_sub_list[[t]] + one_sub_mat)%*%A_mat_list[[t]]
      
      c_t_2 <- matrix(c_1_list[[t]], nrow = 1)%*%W%*%A_mat_list[[t]]%*%t(K_sub_list[[t]] + one_sub_mat) + 
        matrix(c_2_list[[t]], nrow = 1)%*%W_negative%*%(K_subsub_list[[t]] + one_subsub_mat)
        
      c_t_1 <- c_t_1/(as.numeric(ratio_vec[t]))
      c_t_2 <- c_t_2/(as.numeric(ratio_vec[t]))
      c_combine <- c(c_combine, -as.vector(c_t_1), -as.vector(c_t_2), as.vector(c_t_2))
    }
    c_combine <- -(c_combine + d_vec)
    
    # solve quadratic problem
    sv <- solve_osqp(P = H, q = -c_combine, A = A_cons, l = l_cons, u = u_cons, 
                     pars = osqpSettings(verbose = FALSE, eps_abs = 1e-06, eps_rel = 1e-06, eps_prim_inf = 1e-06, eps_dual_inf = 1e-06))
    v <- sv$x
    
    if(length(v)!=T*(N+2*N_negative)){ stop('The dimension of dual problem solution is problematic')}
    
    u_list <- list()
    lp_list <- list()
    ln_list <- list()
    for(t in 1:T){
      v_all <- v[((t-1)*(N+2*N_negative)+1):(t*(N+2*N_negative))]
      
      u_list[[t]] <- v_all[1:N] %>% matrix(ncol = 1)
      lp_list[[t]] <- v_all[(N+1):(N+N_negative)] %>% matrix(ncol = 1)
      ln_list[[t]] <- v_all[(N+N_negative+1):(N+2*N_negative)] %>% matrix(ncol = 1)
    }
    
    beta_new <- list()
    beta0_new <- list()
    for(t in 1:T){
      beta_new[[t]] <- t(H_list[[t]])%*%A_mat_list[[t]]%*%(u_list[[t]] - W%*%c_1_list[[t]]) + 
        t(H_sub_list[[t]])%*%(lp_list[[t]] - ln_list[[t]] - W_negative%*%c_2_list[[t]])
      
      beta0_new[[t]] <- one%*%A_mat_list[[t]]%*%(u_list[[t]] - W%*%c_1_list[[t]]) + 
        one_sub%*%(lp_list[[t]] - ln_list[[t]] - W_negative%*%c_2_list[[t]])
      
      beta_new[[t]] <- beta_new[[t]]/(as.numeric(ratio_vec[t]))
      beta0_new[[t]] <- beta0_new[[t]]/(as.numeric(ratio_vec[t]))
    }
    
    diff <- c()
    for(t in 1:T){
      diff <- c(diff, beta_new[[t]]-beta_old[[t]], beta0_new[[t]]-beta0_old[[t]])
    }
    # print(max(abs(diff)))
    
    s <- s + 1
    if(max(abs(diff)) <= epsilon.max) break
    if(s>=iteration.max){
      warning(paste("Algorithm doesn't converge within", iteration.max, "iterations", sep = ' '))
      break
    }
    
    beta_old <- beta_new
    beta0_old <- beta0_new
  }
  
  return(list(beta = beta_new, beta0 = beta0_new))
}

## Estimate optimal treatment assignment under estimated rule
est_opt <- function(out, H_testing = NULL){
  T <- length(H_testing)
  beta <- out$beta
  beta0 <- out$beta0
  
  A_est <- list()
  for(t in 1:T){
    A_est[[t]] <- cbind(1, H_testing[[t]]) %*% (c(beta0[[t]], beta[[t]]) %>% matrix(ncol = 1)) %>% sign() %>% as.vector()
  }
  return(A_est)
}

## Estimate objective function under estimated rule
est_value <- function(o, R_testing, P_testing, A_testing, H_testing, mu2, type = 'smooth'){
  
  T <- length(A_testing)
  N_testing <- length(R_testing)
  
  f_1 <- matrix(0, N_testing, T)
  f_0 <- matrix(0, N_testing, T)
  
  for(t in 1:T){
    f_1[,t] <- (A_testing[[t]] %>% matrix(ncol = 1))*(H_testing[[t]] %*% (o$beta[[t]] %>% matrix(ncol = 1)) + o$beta0[[t]] %>% as.vector())
  }
  f_0 <- f_1
  f_1 <- cbind(f_1, 1)
  f_0 <- cbind(f_0, 0)
  
  l1 <- f_1 %>% apply(MARGIN = 1, FUN = min)
  l0 <- f_0 %>% apply(MARGIN = 1, FUN = min)
  
  l <- l1 - l0
  reward <- (R_testing*l/P_testing) %>% mean()
  
  return(reward)
}
