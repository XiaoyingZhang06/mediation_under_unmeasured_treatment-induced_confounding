library(doParallel)
library(foreach)
library(nleqslv)
time_start <- Sys.time()
expit <- function(x){exp(x)/{1 + exp(x)}}

data_generate <- function(seed, n){
  
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
  X <- rnorm(n,0,1)
  prob_A <- expit(beta_A0 + beta_AX*X) 
  A <- rbinom(n,1,prob_A)
  L <- rnorm(n, beta_L0 + beta_LA*A + beta_LX*X, sigma_L)
  U <- rnorm(n, beta_U0 + beta_UA*A + beta_UL*L + beta_UX*X, sigma_U)
  M <- rnorm(n, beta_M0 + beta_MA*A + beta_MX*X + beta_ML*L + beta_MU*U, sigma_M)
  W <- rnorm(n, beta_W0 + beta_WX*X + beta_WL*L + beta_WU*U, sigma_W)
  Z <- rnorm(n, beta_Z0 + beta_ZA*A + beta_ZX*X + beta_ZL*L + beta_ZU*U, sigma_Z)
  Y <- rnorm(n, beta_Y0 + beta_YA*A + beta_YM*M + beta_YX*X + beta_YL*L + beta_YU*U, sigma_Y)
  data <- data.frame(X,A,U,L,M,Z,W,Y)
  return(data)
}

# ### true value
# true_value <- function(a,a_star){
#   beta_L0 <- 0.2 ; beta_LA <- -0.25 ; beta_LX <- 0.3
#   beta_U0 <- 0.5 ; beta_UA <- -0.3 ; beta_UX <- 0.8 ; beta_UL <- 0.7
#   beta_M0 <- -0.6 ; beta_MA <- 0.7 ; beta_MX <- -0.4 ; beta_ML <- 0.25 ; beta_MU <- 0.5
#   beta_Y0 <- 0.2 ; beta_YA <- 0.5 ; beta_YX <- -0.4 ; beta_YL <- -0.3 ; beta_YU <- 1 ; beta_YM <- 0.8
#   data <- data_generate(sample(1:1000,1),10000000)
#   with(data, {
#     L_a <- beta_L0 + beta_LA*a + beta_LX*X
#     U_a <- beta_U0 + beta_UA*a + beta_UX*X + beta_UL*L_a
#     L_astar <- beta_L0 + beta_LA*a_star + beta_LX*X
#     U_astar <- beta_U0 + beta_UA*a_star + beta_UX*X + beta_UL*L_astar
#     M_astar <- beta_M0 + beta_MA*a_star + beta_MX*X + beta_ML*L_astar + beta_MU*U_astar
#     Y_a_M_astar <- beta_Y0 + beta_YA*a + beta_YM*M_astar + beta_YX*X + beta_YL*L_a + beta_YU*U_a
#     mean(Y_a_M_astar)
#   })
# }
# Y10_true <- true_value(1,0) # 0.696

point_est <- function(data, mis_A = FALSE, mis_W = FALSE, mis_M = FALSE, mis_h = FALSE, mis_q = FALSE){
  
  X <- with(data, X)
  A <- with(data, A)
  M <- with(data, M)
  Y <- with(data, Y)
  Z <- with(data, Z)
  W <- with(data, W)
  L <- with(data, L)
  
  # P(A|X)
  if (mis_A) {
    model_A <- glm(A ~ abs(X-X^2+1), family = binomial)
    prob_A_hat <- as.numeric(expit(cbind(1,abs(X-X^2+1)) %*% coef(model_A)))
  } else {
    model_A <- glm(A ~ X, family = binomial)
    prob_A_hat <- as.numeric(expit(cbind(1,X) %*% coef(model_A)))
  }
  pa_X <- prob_A_hat^a*(1-prob_A_hat)^(1-a)
  pastar_X <- prob_A_hat^a_star*(1-prob_A_hat)^(1-a_star)
  
  # P(W,L|a,X)
  if (mis_W) {
    model_LaX <- lm(L[A==a] ~ abs(X-X^2+1)[A==a])
    sigma_LaX <- sigma(model_LaX)
    mu_LaX <- as.numeric(cbind(1,abs(X-X^2+1)) %*% coef(model_LaX))
    model_WaX <- lm(W[A==a] ~ abs(X-X^2+1)[A==a])
    sigma_WaX <- sigma(model_WaX)
    mu_WaX <- as.numeric(cbind(1,abs(X-X^2+1)) %*% coef(model_WaX))
  } else {
    model_LaX <- lm(L[A==a] ~ X[A==a])
    sigma_LaX <- sigma(model_LaX)
    mu_LaX <- as.numeric(cbind(1,X) %*% coef(model_LaX))
    model_WaX <- lm(W[A==a] ~ X[A==a])
    sigma_WaX <- sigma(model_WaX)
    mu_WaX <- as.numeric(cbind(1,X) %*% coef(model_WaX))
  }
  
  # f(M|a_star,X)
  if (mis_M) {
    model_MastarX <- lm(M[A==a_star] ~ abs(X-X^2+1)[A==a_star])
    sigma_MastarX <- sigma(model_MastarX)
    mu_MastarX <- as.numeric(cbind(1,abs(X-X^2+1)) %*% coef(model_MastarX))
  } else {
    model_MastarX <- lm(M[A==a_star] ~ X[A==a_star])
    sigma_MastarX <- sigma(model_MastarX)
    mu_MastarX <- as.numeric(cbind(1,X) %*% coef(model_MastarX))
  }
  fM_astarX <- dnorm(M,mu_MastarX,sigma_MastarX)
  
  # h_a(W,L,M,X)
  if (mis_h) {
    h_equation <- function(gamma) {
      h <- gamma[1] + gamma[2]*W + gamma[3]*L + gamma[4]*M + gamma[5]*abs(X-X^2+1)
      f1 <- mean((Y-h)*A^a*(1-A)^(1-a))
      f2 <- mean((Y-h)*A^a*(1-A)^(1-a)*Z)
      f3 <- mean((Y-h)*A^a*(1-A)^(1-a)*L)
      f4 <- mean((Y-h)*A^a*(1-A)^(1-a)*M)
      f5 <- mean((Y-h)*A^a*(1-A)^(1-a)*abs(X-X^2+1))
      c(f1,f2,f3,f4,f5)
    }
    gamma_start <- c(0,0,0,0,0)
    gamma_hat <- nleqslv(gamma_start, h_equation)$x
    h_hat <- as.numeric(cbind(1,W,L,M,abs(X-X^2+1)) %*% gamma_hat)
  } else {
    h_equation <- function(gamma) {
      h <- gamma[1] + gamma[2]*W + gamma[3]*L + gamma[4]*M + gamma[5]*X
      f1 <- mean((Y-h)*A^a*(1-A)^(1-a))
      f2 <- mean((Y-h)*A^a*(1-A)^(1-a)*Z)
      f3 <- mean((Y-h)*A^a*(1-A)^(1-a)*L)
      f4 <- mean((Y-h)*A^a*(1-A)^(1-a)*M)
      f5 <- mean((Y-h)*A^a*(1-A)^(1-a)*X)
      c(f1,f2,f3,f4,f5)
    }
    gamma_start <- c(0,0,0,0,0)
    gamma_hat <- nleqslv(gamma_start, h_equation)$x
    h_hat <- as.numeric(cbind(1,W,L,M,X) %*% gamma_hat)
  }
  
  # q_a(Z,L,M,X)
  if (mis_q) {
    q_equation <- function(theta) {
      q <- exp(theta[1]+theta[2]*Z^2+theta[3]*M*Z+theta[4]*abs(X-X^2+1)*Z+theta[5]*M^2+theta[6]*abs(X-X^2+1)*M
               +theta[7]*abs(X-X^2+1)^2+theta[8]*L^2+theta[9]*L*Z+theta[10]*L*M+theta[11]*abs(X-X^2+1)*L+theta[12]*M+theta[13]*L+theta[14]*abs(X-X^2+1)+theta[15]*Z)
      f1 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)))
      f2 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*W^2))
      f3 <- mean(((q*fM_astarX*M-mu_MastarX)*A^a*(1-A)^(1-a)*W))
      f4 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*abs(X-X^2+1)*W))
      f5 <- mean(((q*fM_astarX*M^2-(mu_MastarX^2+sigma_MastarX^2))*A^a*(1-A)^(1-a)))
      f6 <- mean(((q*fM_astarX*M-mu_MastarX)*A^a*(1-A)^(1-a)*abs(X-X^2+1)))
      f7 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*abs(X-X^2+1)^2))
      f8 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*L^2))
      f9 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*L*W))
      f10 <- mean(((q*fM_astarX*M-mu_MastarX)*A^a*(1-A)^(1-a)*L))
      f11 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*abs(X-X^2+1)*L))
      f12 <- mean(((q*fM_astarX*M-mu_MastarX)*A^a*(1-A)^(1-a)))
      f13 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*L))
      f14 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*abs(X-X^2+1)))
      f15 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*W))
      c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15)
    }
    theta_start <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    theta_hat <- nleqslv(theta_start, q_equation)$x
    q_hat <- as.numeric(exp(cbind(1,Z^2,M*Z,abs(X-X^2+1)*Z,M^2,abs(X-X^2+1)*M,abs(X-X^2+1)^2,L^2,L*Z,L*M,abs(X-X^2+1)*L,M,L,abs(X-X^2+1),Z) %*% theta_hat))
  } else {
    q_equation <- function(theta) {
      q <- exp(theta[1]+theta[2]*Z^2+theta[3]*M*Z+theta[4]*X*Z+theta[5]*M^2+theta[6]*X*M
               +theta[7]*X^2+theta[8]*L^2+theta[9]*L*Z+theta[10]*L*M+theta[11]*X*L+theta[12]*M+theta[13]*L+theta[14]*X+theta[15]*Z)
      f1 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)))
      f2 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*W^2))
      f3 <- mean(((q*fM_astarX*M-mu_MastarX)*A^a*(1-A)^(1-a)*W))
      f4 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*X*W))
      f5 <- mean(((q*fM_astarX*M^2-(mu_MastarX^2+sigma_MastarX^2))*A^a*(1-A)^(1-a)))
      f6 <- mean(((q*fM_astarX*M-mu_MastarX)*A^a*(1-A)^(1-a)*X))
      f7 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*X^2))
      f8 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*L^2))
      f9 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*L*W))
      f10 <- mean(((q*fM_astarX*M-mu_MastarX)*A^a*(1-A)^(1-a)*L))
      f11 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*X*L))
      f12 <- mean(((q*fM_astarX*M-mu_MastarX)*A^a*(1-A)^(1-a)))
      f13 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*L))
      f14 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*X))
      f15 <- mean(((q*fM_astarX-1)*A^a*(1-A)^(1-a)*W))
      c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15)
    }
    theta_start <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    theta_hat <- nleqslv(theta_start, q_equation)$x
    q_hat <- as.numeric(exp(cbind(1,Z^2,M*Z,X*Z,M^2,X*M,X^2,L^2,L*Z,L*M,X*L,M,L,X,Z) %*% theta_hat))
  }
  
  # g0(W,X),g1(M,X),g(X)
  if (mis_h) {
    g0_hat <- as.numeric(cbind(1,W,L,mu_MastarX,abs(X-X^2+1)) %*% gamma_hat)
    g1_hat <- as.numeric(cbind(1,mu_WaX,mu_LaX,M,abs(X-X^2+1)) %*% gamma_hat)
    g_hat <- as.numeric(cbind(1,mu_WaX,mu_LaX,mu_MastarX,abs(X-X^2+1)) %*% gamma_hat)
  } else {
    g0_hat <- as.numeric(cbind(1,W,L,mu_MastarX,X) %*% gamma_hat)
    g1_hat <- as.numeric(cbind(1,mu_WaX,mu_LaX,M,X) %*% gamma_hat)
    g_hat <- as.numeric(cbind(1,mu_WaX,mu_LaX,mu_MastarX,X) %*% gamma_hat)
  }
  
  # psi(a,a_star)
  indicator_a <- as.numeric(A == a)
  indicator_astar <- as.numeric(A == a_star)
  est_MR <- mean(indicator_a/pa_X*fM_astarX*q_hat*(Y-h_hat) 
                 + indicator_a/pa_X*(g0_hat-g_hat)
                 + indicator_astar/pastar_X*(g1_hat-g_hat) 
                 + g_hat)
  est_M1 <- mean(g_hat)
  est_M2 <- mean(indicator_astar/pastar_X*g1_hat)
  est_M3 <- mean(indicator_a/pa_X*g0_hat)
  est_M4 <- mean(indicator_a/pa_X*fM_astarX*q_hat*Y)
  est <- c(MR = est_MR, M1 = est_M1, M2 = est_M2, M3 = est_M3, M4 = est_M4)
  return(est)
}

boot <- function(data, mis_A = FALSE, mis_W = FALSE, mis_M = FALSE, mis_h = FALSE, mis_q = FALSE){
  
  result <- data.frame(matrix(NA,nrow=n_boot,ncol=5))
  for (i_boot in 1:n_boot){
    indices <- sample(1:n, size = n, replace = TRUE)
    boot_sample <- data[indices, ]
    result[i_boot,] <- point_est(boot_sample, mis_A, mis_W, mis_M, mis_h, mis_q)
  }
  est <- point_est(data, mis_A, mis_W, mis_M, mis_h, mis_q)
  std <- apply(result,2,sd)
  bias <- est - true_value
  MSE <- std^2 + bias^2
  cover <- true_value <= est + 1.96 * std & true_value >= est - 1.96 * std
  return(list(estimates=est, std_errors=std, bias=bias, MSE=MSE, cover=cover))
}

n <- 1000
n_sim <- 1000
n_boot <- 200
true_value <- 0.696
a <- 1; a_star <- 0
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
n_cores <- detectCores()
cl <- makeCluster(n_cores)
print(paste("detectCores:",length(cl)))
registerDoParallel(cl)
clusterExport(cl, c("data_generate", "expit", "point_est", "n", "true_value", "a", "a_star", "boot", "n_boot"))
settings <- list(
  all_correct = function(data) boot(data),
  M1 = function(data) boot(data, mis_A = 1, mis_q = 1),
  M2 = function(data) boot(data, mis_M = 1, mis_q = 1),
  M3 = function(data) boot(data, mis_W = 1, mis_q = 1),
  M4 = function(data) boot(data, mis_W = 1, mis_h = 1)
)
combine_results <- function(result1, result2){
  for (setting_name in names(settings)) {
    result1[[setting_name]]$estimates <- rbind(result1[[setting_name]]$estimates,result2[[setting_name]]$estimates)
    result1[[setting_name]]$std_errors <- rbind(result1[[setting_name]]$std_errors,result2[[setting_name]]$std_errors)
    result1[[setting_name]]$bias <- rbind(result1[[setting_name]]$bias,result2[[setting_name]]$bias)
    result1[[setting_name]]$MSE <- rbind(result1[[setting_name]]$MSE,result2[[setting_name]]$MSE)
    result1[[setting_name]]$cover <- rbind(result1[[setting_name]]$cover,result2[[setting_name]]$cover)
  }
  return(result1)
}

results_list <- foreach(i = 1:n_sim,
                        .combine = combine_results,
                        .packages = c("nleqslv")) %dopar% {
                          
                          tryCatch({
                            data <- data_generate(i, n)
                            results_list <- list()
                            for (setting_name in names(settings)) {
                              res <- settings[[setting_name]](data)
                              results_list[[setting_name]]$estimates <- res$estimates
                              results_list[[setting_name]]$std_errors <- res$std_errors
                              results_list[[setting_name]]$bias <- res$bias
                              results_list[[setting_name]]$MSE <- res$MSE
                              results_list[[setting_name]]$cover <- res$cover
                            }
                            return(results_list)
                          }, error = function(e) {
                            warning("iteration ", i, " error: ", e$message)
                            error_results <- list()
                            for (setting_name in names(settings)) {
                              error_results[[setting_name]] <- list(
                                estimates = rep(NA,5),
                                std_errors = rep(NA,5),
                                bias = rep(NA,5),
                                MSE = rep(NA,5),
                                cover = rep(NA,5)
                              )
                            }
                            return(error_results)
                          })
                        }

stopCluster(cl)

### summary
summary_stats <- list()
for (setting_name in names(settings)) {
  summary_stats[[setting_name]] <- list(
    mean_estimates = colMeans(results_list[[setting_name]]$estimates, na.rm=TRUE),
    mean_std_errors = colMeans(results_list[[setting_name]]$std_errors, na.rm=TRUE),
    mean_MSE = colMeans(results_list[[setting_name]]$MSE, na.rm=TRUE),
    mean_bias = colMeans(results_list[[setting_name]]$bias, na.rm=TRUE),
    coverage = colMeans(results_list[[setting_name]]$cover, na.rm=TRUE)
  )
}
time_end <- Sys.time()
time <- difftime(time_end,time_start)
print(time)
save.image(file = "simresults_continuousM_N1000_para.RData")
