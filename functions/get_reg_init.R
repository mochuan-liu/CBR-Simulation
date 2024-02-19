get_reg_init <- function(K_list, H_list = NULL, A_list, Y, kernel = 'linear', sigma_list = NULL){
  
  suppressMessages(require(glmnet))
  
  T <- length(A_list)
  
  beta  <- list()
  beta0 <- list()
  
  if(is.null(K_list)){
    K_list <- K_matrix_list(H_list, kernel = kernel, sigma_list = sigma_list)
  }
  
  for(t in T:1){
    K <- K_list[[t]]
    A <- A_list[[t]]
    
    X <- cbind(K, A, sweep(K, 1, A, "*"))
    model.cv <- cv.glmnet(X, Y, family = 'gaussian', alpha = 1, intercept = TRUE, grouped = FALSE)
    model <- glmnet(X, Y,  family = 'gaussian', alpha = 1, lambda = model.cv$lambda.min, intercept = TRUE, standardize = FALSE)
    
    p <- dim(K_list[[t]])[1]
    
    coe <- model$beta
    beta0[[t]] <- coe[p+1,]
    beta[[t]]  <- coe[(p+2):(2*p+1)]
  }
  
  return(list(beta = beta, beta0 = beta0))
}

