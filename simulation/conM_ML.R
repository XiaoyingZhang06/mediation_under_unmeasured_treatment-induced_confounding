library(doParallel)
library(foreach)
library(nnet)
library(caret)
library(rdetools)
library(expm)
library(MASS)
library(ks)
time_start <- Sys.time()

data_generate <- function(seed, n, return_true = FALSE){
  
  expit <- function(x){exp(x)/{1 + exp(x)}}
  # parameters
  beta_A0 <- 0 ; beta_AX <- 1
  beta_L0 <- 0.2 ; beta_LA <- -0.25 ; beta_LX <- 0.3
  beta_U0 <- 0.5 ; beta_UA <- -0.3 ; beta_UX <- 0.8 ; beta_UL <- 0.7
  beta_M0 <- -0.6 ; beta_MA <- 0.7 ; beta_MX <- -0.4 ; beta_ML <- 0.25 ; beta_MU <- 0.5
  beta_Z0 <- -0.3 ; beta_ZA <- 0.4 ; beta_ZX <- 0.5 ; beta_ZL <- -0.6 ; beta_ZU <- 0.5
  beta_W0 <- 0.5 ; beta_WX <- -0.15 ; beta_WL <- -0.6 ; beta_WU <- 0.5
  beta_Y0 <- 0.2 ; beta_YA <- 0.5 ; beta_YX <- -0.6 ; beta_YL <- -0.3 ; beta_YU <- 1 ; beta_YM <- 0.8
  sigma_U <- 0.5; sigma_L <- 0.25; sigma_W <- 0.25; sigma_Z <- 0.25; sigma_Y <- 0.5; sigma_M <- 1
  
  # data generate
  set.seed(seed)
  X_prime <- rnorm(n,0,1)
  prob_A <- expit(beta_A0 + beta_AX*X_prime)
  A <- rbinom(n,1,prob_A)
  L_prime <- rnorm(n, beta_L0 + beta_LA*A + beta_LX*X_prime, sigma_L)
  U <- rnorm(n, beta_U0 + beta_UA*A + beta_UL*L_prime + beta_UX*X_prime, sigma_U)
  M <- rnorm(n, beta_M0 + beta_MA*A + beta_MX*X_prime + beta_ML*L_prime + beta_MU*U + 0.5*A*U - 0.5*X_prime*L_prime, sigma_M)
  W <- rlnorm(n, beta_W0 + beta_WX*X_prime + beta_WL*L_prime + beta_WU*U, sigma_W)
  Z <- rlnorm(n, beta_Z0 + beta_ZA*A + beta_ZX*X_prime + beta_ZL*L_prime + beta_ZU*U, sigma_Z)
  Y <- rnorm(n, beta_Y0 + beta_YA*A + beta_YM*M + beta_YX*X_prime + beta_YL*L_prime + beta_YU*U + 0.5*A*U - 0.5*X_prime*L_prime, sigma_Y)
  obs_data <- data.frame(X=X_prime^3,A,U,L=L_prime^3,M,Z,W,Y)
  true_data <- data.frame(X,A,U,L,M,Z,W,Y)
  if (return_true) {
    return(true_data)
  } else {
    return(obs_data)
  }
  
}


# ### true value
# true_value <- function(a,a_star){
#   beta_L0 <- 0.2 ; beta_LA <- -0.25 ; beta_LX <- 0.3
#   beta_U0 <- 0.5 ; beta_UA <- -0.3 ; beta_UX <- 0.8 ; beta_UL <- 0.7
#   beta_M0 <- -0.6 ; beta_MA <- 0.7 ; beta_MX <- -0.4 ; beta_ML <- 0.25 ; beta_MU <- 0.5
#   beta_Y0 <- 0.2 ; beta_YA <- 0.5 ; beta_YX <- -0.6 ; beta_YL <- -0.3 ; beta_YU <- 1 ; beta_YM <- 0.8
#   data <- data_generate(sample(1:1000,1),10000000,return_true=TRUE)
#   with(data, {
#     L_a <- beta_L0 + beta_LA*a + beta_LX*X
#     U_a <- beta_U0 + beta_UA*a + beta_UX*X + beta_UL*L_a
#     L_astar <- beta_L0 + beta_LA*a_star + beta_LX*X
#     U_astar <- beta_U0 + beta_UA*a_star + beta_UX*X + beta_UL*L_astar
#     M_astar <- beta_M0 + beta_MA*a_star + beta_MX*X + beta_ML*L_astar + beta_MU*U_astar + 0.5*a_star*U_astar - 0.5*X*L_astar
#     Y_a_M_astar <- beta_Y0 + beta_YA*a + beta_YM*M_astar + beta_YX*X + beta_YL*L_a + beta_YU*U_a + 0.5*a*U_a - 0.5*X*L_a
#     mean(Y_a_M_astar)
#   })
# }
# Y10_true <- true_value(1,0) # 0.5084


point_est <- function(data){
  
  est_MR <- numeric(CF_K_fold)
  est_M1 <- numeric(CF_K_fold)
  est_M2 <- numeric(CF_K_fold)
  est_M3 <- numeric(CF_K_fold)
  est_M4 <- numeric(CF_K_fold)
  kernel_sigma <- 15
  appro_rate <- 0.05
  lm_Hh <- 0.01
  lm_Qh <- 1
  lm_Hq <- 1
  lm_Qq <- 0.01
  
  n <- nrow(data)
  pseudo_Y <- numeric(n)
  random_sample <- createFolds(1:n, k = CF_K_fold)
  for (cf_rep in 1:CF_K_fold) {
    te_idx <- random_sample[[cf_rep]]
    tr_idx  <- setdiff(1:n, te_idx)
    tr_data <- data[-random_sample[[cf_rep]],]
    te_data <- data[random_sample[[cf_rep]],]
    n_tr <- n - length(random_sample[[cf_rep]])
    n_te <- length(random_sample[[cf_rep]])
    
    # P(a|X) & P(a_star|X)
    model_A <- nnet(factor(A) ~ X, data = tr_data, size = 3, decay = 0.1, maxit = 500, trace = FALSE)
    prob_A_hat <- as.numeric(predict(model_A, newdata = data.frame(X = te_data$X), type = "raw"))
    pa_X <- prob_A_hat^a*(1-prob_A_hat)^(1-a)
    pastar_X <- prob_A_hat^a_star*(1-prob_A_hat)^(1-a_star)
    
    # P(W,L|a,X)
    W <- tr_data[tr_data$A == a,]$W
    L <- tr_data[tr_data$A == a,]$L
    X <- tr_data[tr_data$A == a,]$X
    model_WL <- nnet(x = X, y = matrix(c(W,L), ncol = 2), size = 3, decay = 0.1, maxit = 500, linout = TRUE, trace = FALSE)
    muWL_aX <- predict(model_WL, newdata = matrix(te_data$X, ncol = 1))
    resid_WL <- cbind(tr_data$W[tr_data$A == a],tr_data$L[tr_data$A == a]) - predict(model_WL, newdata = matrix(tr_data[tr_data$A == a,]$X, ncol = 1))
    kde_WL <- kde(resid_WL)
    
    # P(M|a_star,X)
    model_M <- nnet(M ~ X, data = tr_data[tr_data$A == a_star,], size = 3, decay = 0.1, maxit = 500, linout = TRUE, trace = FALSE)
    muM_astarX <- as.numeric(predict(model_M, newdata = data.frame(X = te_data$X)))
    muM_astarX_tr <- as.numeric(predict(model_M, newdata = data.frame(X = tr_data$X)))
    resid_M <- tr_data$M[tr_data$A == a_star] - predict(model_M, newdata = data.frame(X = tr_data[tr_data$A == a_star,]$X))
    kde_M <- kde(resid_M)
    pM_astarX <- rep(0,n_te)
    for (i in 1:n_te){
      resid_i <- te_data$M[i] - muM_astarX[i]
      pM_astarX[i] <- approx(kde_M$eval.points, kde_M$estimate, xout = resid_i, rule = 2)$y
    }
    pM_astarX_tr <- rep(0,n_tr)
    for (i in 1:n_tr){
      resid_i <- tr_data$M[i] - muM_astarX_tr[i]
      pM_astarX_tr[i] <- approx(kde_M$eval.points, kde_M$estimate, xout = resid_i, rule = 2)$y
    }
    
    # h_a(W,L,M,X) & q_a(Z,L,M,X)
    ## Gaussian RBF kernel ############################################################
    h_arg <- cbind(tr_data$W, tr_data$A, tr_data$L, tr_data$M, tr_data$X)
    K_H1 <- rbfkernel(h_arg, sigma = kernel_sigma)
    q_arg <- cbind(tr_data$Z, tr_data$A, tr_data$L, tr_data$M, tr_data$X)
    K_Q <- rbfkernel(q_arg, sigma = kernel_sigma)
    n_sample <- 100
    K_H2 <- matrix(0, n_tr, n_tr)
    features_i <- as.matrix(tr_data[, c("W", "A", "L", "M", "X")])
    for (j in 1:n_tr) {
      M_j_samples <- rkde(n_sample, kde_M) + muM_astarX_tr[j]
      M_mat <- cbind(tr_data$W[j], tr_data$A[j], tr_data$L[j], M_j_samples, tr_data$X[j])
      kernel_mat <- rbfkernel(X = features_i, Y = M_mat, sigma = kernel_sigma)
      K_H2[, j] <- rowMeans(kernel_mat)
    }
    ###################################################################################
    
    ## Optimization with approximation###################################################################
    r <- max(10, ceiling(n_tr * appro_rate))
    S <- diag(1, n_tr)[, sort(sample(1:n_tr, r))]
    ident_n <- diag(1, n_tr)
    ident_r <- diag(1, r)
    pM <- diag(pM_astarX_tr)
    ## Optimization ###################################################################
    sqrtmtSK_QS <- sqrtm(ginv(t(S) %*% K_Q %*% S))
    D <- K_Q %*% S %*% sqrtmtSK_QS
    sqrtmtSK_HS <- sqrtm(ginv(t(S) %*% K_H1 %*% S))
    V <- K_H1 %*% S %*% sqrtmtSK_HS
    ## Optimization h ##
    Gamm <- D %*% ginv(t(D) %*% D/n_tr + lm_Qh/n_tr^0.8 * ident_r) %*% t(D)
    alpha_h <- ginv(t(V) %*% Gamm %*% V + n_tr^2*lm_Hh/n_tr^0.8 * ident_r) %*%
      t(V) %*% Gamm %*% tr_data$Y
    ## Optimization q ##
    Gamm <- V %*% ginv(t(V) %*% V/n_tr + lm_Hq/n_tr^0.8 * ident_r) %*% t(V)
    approx_inv <- 1/(lm_Hq/n_tr^0.8) * ident_n - (1/n_tr)/(lm_Hq/n_tr^0.8)^2 * V %*% 
      ginv(ident_r + (1/n_tr)/(lm_Hq/n_tr^0.8) * t(V) %*% V) %*% t(V)
    alpha_q <- ginv(t(D) %*% t(pM) %*% Gamm %*% pM %*% D + n_tr^2*lm_Qq/n_tr^0.8 * ident_r) %*%
      t(D) %*% t(pM) %*% approx_inv %*% K_H2 %*% rep(1, n_tr)
    ###################################################################################
    
    ## Estimate function h_hat ########################################################
    h_hat <- function(w, a, l, m, x){
      W <- cbind(w, a, l, m, x)
      K <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(W, nrow(W), ncol(W)))
      return(t(K) %*% S %*% sqrtmtSK_HS %*% alpha_h)
    }
    ###################################################################################
    
    ## Estimate function q_hat ########################################################
    q_hat <- function(z, a, l, m, x){
      Z <- cbind(z, a, l, m, x)
      K <- rbfkernel(q_arg, sigma = kernel_sigma, matrix(Z, nrow(Z), ncol(Z)))
      return(t(K) %*% S %*% sqrtmtSK_QS %*% alpha_q)
    }
    ###################################################################################
    
    h_est <- h_hat(te_data$W,a,te_data$L,te_data$M,te_data$X)
    q_est <- q_hat(te_data$Z,a,te_data$L,te_data$M,te_data$X)
    
    # g0(W,X),g1(M,X),g(X)
    g0_est <- rep(0,n_te); g1_est <- rep(0,n_te); g_est <- rep(0,n_te)
    n_sample <- 100
    for (i in 1:n_te){
      sample_WL <- rkde(n_sample, kde_WL)
      sample_W <-  sample_WL[,1] + muWL_aX[i,1]
      sample_L <-  sample_WL[,2] + muWL_aX[i,2]
      sample_M <- rkde(n_sample, kde_M) + muM_astarX[i]
      g0_est[i] <- mean(h_hat(te_data$W[i],a,te_data$L[i],sample_M,te_data$X[i]), na.rm = TRUE)
      g1_est[i] <- mean(h_hat(sample_W,a,sample_L,te_data$M[i],te_data$X[i]), na.rm = TRUE)
      g_est[i] <- mean(h_hat(sample_W,a,sample_L,sample_M,te_data$X[i]), na.rm = TRUE)
    }
    ###################################################################################
    
    # psi(a,a_star)
    indicator_a <- as.numeric(te_data$A == a)
    indicator_astar <- as.numeric(te_data$A == a_star)
    pseudo_Y[te_idx] <- as.numeric(indicator_a/pa_X*pM_astarX*q_est*(te_data$Y-h_est) 
                                   + indicator_a/pa_X*(g0_est-g_est)
                                   + indicator_astar/pastar_X*(g1_est-g_est) 
                                   + g_est)
  }
  
  est <- mean(pseudo_Y)
  bias <- est - true_value
  se <- sd(pseudo_Y) / sqrt(n)
  lower <- est - 1.96 * se
  upper <- est + 1.96 * se
  cover <- true_value <= upper & true_value >= lower
  res <- c(bias = bias, se = se, cover = cover)
  return(res)
  
}

n_vec <- c(250,500,1000)
n_sim <- 1000
CF_K_fold <- 5
a <- 1; a_star <- 0
true_value <- 0.5084
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
n_cores <- detectCores()
cl <- makeCluster(getOption('cl.cores',n_cores))
print(paste("detectCores:",length(cl)))
clusterExport(cl, c("data_generate", "point_est", "n_vec", "true_value", "CF_K_fold", "a", "a_star"))
registerDoParallel(cl)
options(warn = -1)


results_list <- list()
for (n in n_vec) {
  cat("Starting simulations for n =", n, "\n")
  res <- foreach(i = 1:n_sim,
                          .combine = rbind,
                          .packages = c("nnet", "caret", "rdetools", "MASS", "expm", "ks")) %dopar% {
                            tryCatch({
                              data <- data_generate(i, n)
                              point_est(data)
                            }, error = function(e) {
                              cat("Iteration", i, "for n =", n, "failed:", e$message, "\n")
                              return(c(bias = NA, se = NA, cover = NA))
                            })
                          }
  res_df <- as.data.frame(res)
  res_df$n <- n
  results_list[[paste0('n',as.character(n))]] <- res_df
}
stopCluster(cl)


summary_by_n <- aggregate(cbind(bias, se, cover) ~ n, data = do.call(rbind, results_list), 
                          FUN = function(x) c(mean = mean(x, na.rm = TRUE)))
print(summary_by_n)
time_end <- Sys.time()
time <- difftime(time_end,time_start)
print(time)
save.image(paste("nonlinear_conM_ML",".RData",sep=''))
