# ## Require R packages kernlab, quadprog, osqp, caret, dplyr,
# ## glmnet and DTRlearn2.
# install.packages(c('kernlab', 'quadprog', 'osqp', 'caret', 'dplyr', 'glmnet', 'DTRlearn2', 'optparse', 'doParallel'))
library(kernlab)
library(quadprog)
library(osqp)
library(caret)
library(dplyr)
library(glmnet)
library(DTRlearn2)

## Set working directory
working_path <- 'your_own_working_path'   # change to your own working directory 
saving_path  <- 'your_own_saving_path'    # change to your own saving directory 
setwd(working_path)

## Load CBR functions 
source('./functions/base_code.R')
source('./functions/get_reg_init.R')
source('./functions/bisection_mrl.R')
source('./functions/bisection_owl.R')
source('./functions/bisection_qlearning.R')

## Load simulation functions 
source('./model_2/complex_model_2.R')

## Load testing data
k <- 6   # Get true optimal treatment assignments under tau = 1.5
testing_list <- readRDS('./data/complex_testing_2.rds')
tau <- testing_list$O[[k]]$tau
A1_testing_opt <- testing_list$O[[k]]$A1_opt
A2_testing_opt <- testing_list$O[[k]]$A2_opt
H1 <- testing_list$H1
H2 <- testing_list$H2
H_testing_list <- list(H1 %>% scale(), H2 %>% scale())
A_testing_list <- list(testing_list$A1, testing_list$A2)
P_testing_list <- list(testing_list$P1, testing_list$P2)

# Argument setting
suppressPackageStartupMessages(library(optparse))
option_list <- list( 
  make_option(c("--seed"), type  = "double", default = 2233, help = "Set seed [default %default]"),
  make_option(c("--method"), type  = "character", default = 'mrl', help = "Method [default %default]"),
  make_option(c("--N"), type  = "double", default = 200, help = "Size of training data [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Set parallel setting
suppressPackageStartupMessages(library(doParallel))
numCores <- as.numeric(Sys.getenv("SLURM_NTASKS"))
if(Sys.getenv("SLURM_NTASKS") == ''){
  numCores = 1
}
registerDoParallel(cores=numCores)

## Simulation setting 
method <- opt$method
N.training <- opt$N
seed <- opt$seed
epsilon <- 10^(-3)
eta <- 10^(-4)
mu2_list <- 10^(-8)
lambda_list <- list(2^(seq(-8,8,2)), 2^(seq(-8,8,2)))
kernel <- 'gaussian'
sigma_list <- NULL

set.seed(seed, kind = "L'Ecuyer-CMRG")
Output = foreach(i = 1:numCores, .combine = 'c') %dopar% {
  
  ## Generate training data
  H <- runif(7*N.training, 0, 1) %>% matrix(ncol = 7) %>% as.data.frame()
  colnames(H) <- sapply(1:dim(H)[2], FUN = function(x) paste('X', x, sep = ''))
  X8_based <- runif(N.training, 0.5, 1)
  H$X8  <- X8_based + runif(N.training, -0.5, 0.5)
  H$X9  <- X8_based + runif(N.training, -0.5, 0.5)
  noisy_vec <- noisy(N.training)
  
  H1.training <- H[,c(1:8)]
  A1.training <- (rbinom(N.training, 1, 0.5) - 0.5)*2
  H2.training <- cbind(H, A1.training)
  A2.training <- (rbinom(N.training, 1, 0.5) - 0.5)*2
  P1.training <- rep(0.5, N.training)
  P2.training <- rep(0.5, N.training)
  
  R <- r(H = H, A1 = A1.training, A2 = A2.training) + noisy_vec$noisy_R
  Y <- y(H = H, A1 = A1.training, A2 = A2.training) + noisy_vec$noisy_Y
  
  H_list <- list(H1.training %>% scale() %>% as.matrix(), H2.training %>% scale() %>% as.matrix())
  A_list <- list(A1.training, A2.training)
  P_list <- list(P1.training, P2.training)
  P <- P1.training*P2.training
  Reward <- Y 
  Risk   <- R 
 
  if(method == 'mrl'){
    res_mrl <- try(br2_mrl(H_list = H_list, A_list = A_list, P_list = P_list, P = P,
                           Reward = Reward, Risk = Risk, tau = tau, 
                           kernel = kernel, sigma_list = sigma_list, 
                           mu2 = mu2_list, lambda_list = lambda_list,
                           epsilon = epsilon, eta = eta))
    res <- res_mrl
  }else if(method == 'owl'){
    res_owl <- try(br2_owl(H_list = H_list, A_list = A_list, P_list = P_list, P = P,
                           Reward = Reward, Risk = Risk, tau = tau, 
                           kernel =  kernel, sigma_list = sigma_list, 
                           lambda_list = lambda_list, 
                           epsilon = epsilon, eta = eta, method = 'owl')) 
    res <- res_owl
  }else if(method == 'aowl'){
    res_aowl <- try(br2_owl(H_list = H_list, A_list = A_list, P_list = P_list, P = P,
                            Reward = Reward, Risk = Risk, tau = tau, 
                            kernel =  kernel, sigma_list = sigma_list, 
                            lambda_list = lambda_list, 
                            epsilon = epsilon, eta = eta, method = 'aowl')) 
    res <- res_aowl
  }else if(method == 'Q-learning'){
    res_q <- try(br2_q(H_list = H_list, A_list = A_list, P_list = P_list, P = P,
                       Reward = Reward, Risk = Risk, tau = tau, sigma_list = sigma_list,
                       epsilon = epsilon, eta = eta))
    res <- res_q
  }else{
    sigma_list <- list()
    for(t in 1:length(A_list)){
      sigma_list[[t]] <- rbf_bandwidth(H_list[[t]], A_list[[t]])
    }
    K_list <- K_matrix_list(H_list = H_list, kernel = kernel, sigma_list = sigma_list)
    res_un <- try(bisection_mrl(H_list = H_list, K_list = K_list, A_list = A_list, P = P, Risk = Risk, Reward = Reward, tau = 1000,
                                kernel = kernel, mu2 = 10^(-8), lambda_vec = c(1,1), 
                                init_type = 'model', sigma_list = sigma_list, 
                                epsilon = epsilon, eta = eta))
    res <- res_un
  }
  
  if(!is.character(res)){
    ## Get estimated treatment recommendations and testing reward/risk
    est <- est_opt(res$res, H_training = H_list, H_testing = H_testing_list, kernel = kernel, sigma_list = res$sigma_list)
    
    Reward_testing <- mean((testing_list$Y)*(est[[1]]==testing_list$A1&est[[2]]==testing_list$A2)/(0.25))
    Risk_testing   <- mean((testing_list$R)*(est[[1]]==testing_list$A1&est[[2]]==testing_list$A2)/(0.25))
    Efficacy_ratio <- (Reward_testing+0.289)/(Risk_testing+0.687)
  }else{
    Reward_testing <- NA
    Risk_testing <- NA
    Efficacy_ratio <- NA
  }
  
  out <- list(H_list = H_list, A_list = A_list, P_list = P_list,
              Reward = Reward, Risk = Risk, tau = tau,
              res = res, method = method, 
              Reward_testing = Reward_testing,
              Risk_testing = Risk_testing, 
              Efficacy_ratio = Efficacy_ratio)
  
  res_file_name <- sprintf('Setting_II_gaussian_%s_%s_%s_%s_out.rds', method, N.training, seed, i)
  saveRDS(out, file.path(saving_path, res_file_name))
}
