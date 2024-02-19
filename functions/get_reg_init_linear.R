get_reg_init_linear <- function(H_list, A_list, Y){
  
  suppressMessages(require(glmnet))
  
  T <- length(A_list)
  
  beta  <- list()
  beta0 <- list()
  
  for(t in T:1){
    H <- as.matrix(H_list[[t]])
    A <- A_list[[t]]
    
    X <- cbind(H, A, sweep(H, 1, A, "*"))
    model.cv <- cv.glmnet(X, Y, family = 'gaussian', alpha = 1, intercept = TRUE, grouped = FALSE)
    model <- glmnet(X, Y,  family = 'gaussian', alpha = 1, lambda = model.cv$lambda.min, intercept = TRUE, standardize = FALSE)
    
    p <- dim(H_list[[t]])[2]
    
    coe <- model$beta
    beta0[[t]] <- coe[p+1,]
    beta[[t]]  <- coe[(p+2):(2*p+1)]
  }
  
  return(list(beta = beta, beta0 = beta0))
}