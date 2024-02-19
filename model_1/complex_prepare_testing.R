work_path <- 'your_own_path'
setwd(work_path)

set.seed(1234)
source('model_1/complex_model_1.R')

N <- 5000

H <- runif(7*N, -1, 1) %>% matrix(ncol = 7) %>% as.data.frame()
colnames(H) <- sapply(1:dim(H)[2], FUN = function(x) paste('X', x, sep = ''))

X8_based <- runif(N, -0.5, 0.5)
H$X8  <- X8_based + runif(N, -0.5, 0.5)
H$X9  <- X8_based + runif(N, -0.5, 0.5)

A1 <- (rbinom(N, 1, 0.5) - 0.5)*2
A2 <- (rbinom(N, 1, 0.5) - 0.5)*2

noisy_vec <- noisy(N)

R <- r(H = H, A1 = A1, A2 = A2) + noisy_vec$noisy_R
Y <- y(H = H, A1 = A1, A2 = A2) + noisy_vec$noisy_Y

O <- searching_opt(H = H, A1 = A1, fun_y = y, fun_r = r, tau_list = c(seq(0.1, 0.9, 0.1), seq(1, 2, 0.1), 100), gamma_list = seq(0, 20, 0.02), epsilon = 0.01)

saveRDS(list(O = O, H = H, H1 = H[,c(1:8)], H2 = cbind(H, A1), A1 = A1, A2 = A2, 
             Y = Y, R = R), file.path(work_path, 'data/complex_testing_1.rds'))

cat('Done!\n')
