suppressMessages(library(dplyr))
suppressMessages(library(truncnorm))

# reward/risk function
y <- function(H,A1,A2){
  H <- as.data.frame(H)
  colnames(H) <- sapply(1:dim(H)[2], FUN = function(x) paste('X', x, sep = ''))
  return(1 + 2*H$X2 + A1*(H$X8^2 + 1) + A2*(H$X9^2 + H$X1^2))
}

r <- function(H,A1,A2){
  H <- as.data.frame(H)
  colnames(H) <- sapply(1:dim(H)[2], FUN = function(x) paste('X', x, sep = ''))
  return(2 - H$X2 + A1*(H$X1 + 1) + A2*(A1*H$X9 + 1))
}

noisy <- function(n){
  noisy_base <- rtruncnorm(n, -0.25, 0.25, 0, 1)
  noisy_R <- noisy_base + rtruncnorm(n, -0.5, 0.5, 0, 1)
  noisy_Y <- noisy_base + rtruncnorm(n, -0.5, 0.5, 0, 1)
  return(list(noisy_R = noisy_R, noisy_Y = noisy_Y))
}

# get opt rule under constraint
searching_opt <- function(H, A1, fun_y, fun_r, tau_list, gamma_list = seq(0, 10, 0.01), epsilon = 0.01){ 
  
  K <- length(gamma_list)
  N <- dim(H)[1]
  S <- length(tau_list)
  
  out_list <- list()
  for(s in 1:S){
    out_list[[s]] <- list(tau = tau_list[s], value = -Inf, risk = NA, A2_opt = NULL, A1_opt = NULL, gamma2 = NA, gamma1 = NA)
  }
  
  # Stage 2 pseudo-outcome
  Y21 <- fun_y(H = H, A1 = A1, A2 = rep( 1, N))
  Y20 <- fun_y(H = H, A1 = A1, A2 = rep(-1, N))
  R21 <- fun_r(H = H, A1 = A1, A2 = rep( 1, N))
  R20 <- fun_r(H = H, A1 = A1, A2 = rep(-1, N))
  
  for(i in 1:K){
    gamma1 <- gamma_list[i]
    A2_opt <- sign(Y21-Y20 - gamma1*(R21-R20))
    
    # Stage 1 pseudo-outcome
    Y11 <- fun_y(H = H, A1 = rep( 1, N), A2 = A2_opt)
    Y10 <- fun_y(H = H, A1 = rep(-1, N), A2 = A2_opt)
    R11 <- fun_r(H = H, A1 = rep( 1, N), A2 = A2_opt)
    R10 <- fun_r(H = H, A1 = rep(-1, N), A2 = A2_opt)
    for(j in 1:K){
      gamma2 <- gamma_list[j]
      A1_opt <- sign(Y11-Y10 - gamma2*(R11-R10))
      Y_opt <- fun_y(H = H, A1 = A1_opt, A2 = A2_opt) %>% mean()
      R_opt <- fun_r(H = H, A1 = A1_opt, A2 = A2_opt) %>% mean()
      
      for(s in 1:S){
        if(abs(R_opt - tau_list[[s]])<=epsilon|tau_list[[s]]>=100){
          if(Y_opt >= out_list[[s]]$value){
            out_list[[s]]$value <- Y_opt
            out_list[[s]]$risk  <- R_opt
            out_list[[s]]$A2_opt <- A2_opt
            out_list[[s]]$A1_opt <- A1_opt
            out_list[[s]]$gamma2 <- gamma2
            out_list[[s]]$gamma1 <- gamma1
          }
        }
      }
    }
  }
  
  return(out_list)
}

