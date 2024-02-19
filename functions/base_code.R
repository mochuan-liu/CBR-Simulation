#####################
# Base function for mrl (nonlinear)
#####################

## Global setting
suppressMessages(require(kernlab))
suppressMessages(require(quadprog))
suppressMessages(require(osqp))
suppressMessages(require(caret))
suppressMessages(require(dplyr))
epsilon.max = 0.001
iteration.max = 100

# rbf bandwidth
rbf_bandwidth <- function(H, A, c.sigma = 1.25){
  
  H <- H %>% as.matrix()
  N <- dim(H)[1]
  
  d_m <- matrix(0, N, N)
  a_m <- matrix(TRUE, N, N)
  for(i in 1:N){
    for(j in i:N){
      d_m[i,j] <- norm(H[i,]-H[j,], '2')
      a_m[i,j] <- (A[i]==A[j])
      d_m[j,i] <- d_m[i,j]
      a_m[j,i] <- a_m[i,j]
    }
  }
  d_l <- as.vector(d_m)
  a_l <- as.vector(a_m)
  sigma <- median(d_l[a_l==FALSE])
  
  return(c.sigma*sigma)
}

# Single stage kernel matrix
kernel_matrix <- function(H, H_testing = NULL, sigma = NULL, kernel = 'gaussian'){
  
  H <- as.matrix(H)
  
  if(is.null(H_testing)){
    H_testing <- H
  }else{
    H_testing <- H_testing %>% as.matrix()
  }
  
  if(kernel == 'linear'){
    ker <- vanilladot()
    K_t <- kernelMatrix(ker, x = H_testing, y = H)
  }else if(kernel == 'gaussian'){
    if(is.null(sigma)) stop('Bandwidth sigma is missing')
    ker <- rbfdot(sigma = 1/sigma^2)
    K_t <- kernelMatrix(ker, x = H_testing, y = H)
  }else{
    stop('Kernel is not supported')
  }
  return(K_t)
}

# Compute kernel matrix list K_t
K_matrix_list <- function(H_list, H_testing_list = NULL, sigma_list = NULL, kernel = 'gaussian'){
  
  T <- length(H_list)
  
  K_list <- list()
  
  for(t in 1:T){
    if(kernel == 'gaussian'){
      K_t <- kernel_matrix(H = H_list[[t]], H_testing = H_testing_list[[t]], kernel = kernel, sigma = sigma_list[[t]])
    }else{
      K_t <- kernel_matrix(H = H_list[[t]], H_testing = H_testing_list[[t]], kernel = kernel)
    }
    K_list[[t]] <- K_t
  }
  return(K_list)
}

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
c_vec_1_list <- function(K_list, A_list, d, beta_list, beta0_list, mu2){
  
  T <- length(A_list)
  N <- length(A_list[[1]])
  c_list <- list()
  
  c_mat <- matrix(0, N, T)
  f_mat <- matrix(0, N, T)
  for(t in 1:T){
    for(i in 1:N){
      f_mat[i,t] <- A_list[[t]][i]*(K_list[[t]][i,]%*%beta_list[[t]] + beta0_list[[t]])
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

c_vec_2_list <- function(K_list, A_list, beta_list, beta0_list, mu2){
  
  T <- length(A_list)
  N <- length(A_list[[1]])
  c_list <- list()
  
  c_mat <- matrix(0, N, T)
  f_mat <- matrix(0, N, T)
  for(t in 1:T){
    for(i in 1:N){
      f_mat[i,t] <- (K_list[[t]][i,]%*%beta_list[[t]] + beta0_list[[t]])
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
single_dc <- function(K_list, A_list, R, beta_start, beta0_start, mu2, 
                      lambda_vec, kernel = 'gaussian', sigma_list = NULL, type = 'smooth'){

  N <- length(A_list[[1]])
  T <- length(A_list)
  
  # rescale outcome 
  lambda <- 1
  ratio_vec <- lambda_vec/lambda
  
  # get negative index 
  index_negative <- which(R<0) %>% sort()
  O_negative <- -R[index_negative]/lambda
  O <- abs(R/lambda)
  N_negative <- length(O_negative)
  
  # get subdata 
  A_sub_list <- list()
  for(t in 1:T){
    A_sub_list[[t]] <- A_list[[t]][index_negative]
  }
  one <- matrix(1, nrow = 1, ncol = N)
  one_sub <- matrix(1, nrow = 1, ncol = N_negative)
  one_mat <- matrix(1, nrow = N, ncol = N)
  one_sub_mat <- matrix(1, nrow = N_negative, ncol = N)
  one_subsub_mat <- matrix(1, nrow = N_negative, ncol = N_negative)
  I_sub <- diag(1, N)[, index_negative]
  
  # Change O to positive vector
  d <- ifelse(R>0, 1, 0)
  d_vec <- c()
  for(t in 1:T){
    d_vec <- c(d_vec, -d, rep(1, N_negative), rep(1, N_negative))
  }
  
  # compute kernel matrix 
  K_sub_list <- list()
  K_subsub_list <- list()
  for(t in 1:T){ 
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
    c_1_list <- c_vec_1_list(K_list = K_list, A_list = A_list, d = d,  beta_list = beta_old, beta0_list = beta0_old, mu2 = mu2)
    c_2_list <- c_vec_2_list(K_list = K_sub_list, A_list = A_sub_list, beta_list = beta_old, beta0_list = beta0_old, mu2 = mu2)
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
      beta_new[[t]] <- A_mat_list[[t]]%*%(u_list[[t]] - W%*%c_1_list[[t]]) + 
        I_sub%*%(lp_list[[t]] - ln_list[[t]] - W_negative%*%c_2_list[[t]])
      
      beta0_new[[t]] <- one%*%A_mat_list[[t]]%*%(u_list[[t]] - W%*%c_1_list[[t]]) + 
        one_sub%*%(lp_list[[t]] - ln_list[[t]] - W_negative%*%c_2_list[[t]])
      
      beta_new[[t]] <- beta_new[[t]]/(as.numeric(ratio_vec[t]))
      beta0_new[[t]] <- beta0_new[[t]]/(as.numeric(ratio_vec[t]))
    }
    
    diff <- c()
    for(t in 1:T){
      diff <- c(diff, (beta_new[[t]]-beta_old[[t]]), beta0_new[[t]]-beta0_old[[t]])
    }
    
    s <- s + 1
    # print(max(abs(diff)))
    
    if(max(abs(diff)) <= epsilon.max) break
    if(s>=iteration.max){
      warning(paste("Algorithm doesn't converge within", iteration.max, "iterations", sep = ' '))
      break
    }
    
    beta_old <- beta_new
    beta0_old <- beta0_new
  }
  
  return(list(beta = beta_new, beta0 = beta0_new, sigma_list = sigma_list))
}

## Estimate optimal treatment assignment under estimated rule
est_opt <- function(out, H_training, H_testing = NULL, K_list = NULL, kernel = 'gaussian', sigma_list = NULL){
  beta <- out$beta
  beta0 <- out$beta0
  T <- length(beta)
  
  A_est <- list()
  
  if(is.null(K_list)){
    if(is.null(H_testing)){
      H_testing <- H_training
    }
    for(t in 1:T){
      K <- kernel_matrix(H = H_training[[t]], H_testing = H_testing[[t]], kernel = kernel, sigma = sigma_list[[t]])
      A_est[[t]] <- cbind(1, K) %*% (c(beta0[[t]], beta[[t]]) %>% matrix(ncol = 1)) %>% sign() %>% as.vector()
    }
  }else{
    for(t in 1:T){
      A_est[[t]] <- cbind(1, K_list[[t]]) %*% (c(beta0[[t]], beta[[t]]) %>% matrix(ncol = 1)) %>% sign() %>% as.vector()
    }
  }
  
  return(A_est)
}
