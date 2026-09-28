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
  beta_M0 <- -0.6 ; beta_MA <- 0.7 ; beta_MX <- 1 ; beta_ML <- 0.25 ; beta_MU <- 0.5
  beta_Z0 <- -0.3 ; beta_ZA <- 0.4 ; beta_ZX <- 0.5 ; beta_ZL <- -0.6 ; beta_ZU <- 0.5
  beta_W0 <- 0.5 ; beta_WX <- -0.15 ; beta_WL <- -0.6 ; beta_WU <- 0.5
  beta_Y0 <- 0.2 ; beta_YA <- 0.5 ; beta_YX <- -0.6 ; beta_YL <- -0.3 ; beta_YU <- 1 ; beta_YM <- 0.8
  sigma_U <- 0.5; sigma_L <- 0.25; sigma_W <- 0.25; sigma_Z <- 0.25; sigma_Y <- 0.5
  
  # data generation
  set.seed(2025 + seed)
  X <- rnorm(n,0,1)
  prob_A <- expit(beta_A0 + beta_AX*X)
  A <- rbinom(n,1,prob_A)
  L <- rnorm(n, beta_L0 + beta_LA*A + beta_LX*X, sigma_L)
  U <- rnorm(n, beta_U0 + beta_UA*A + beta_UL*L + beta_UX*X, sigma_U)
  prob_M <- expit(beta_M0 + beta_MA*A + beta_MX*X + beta_ML*L + beta_MU*U)
  M <- rbinom(n, 1, prob_M)
  W <- rnorm(n, beta_W0 + beta_WX*X + beta_WL*L + beta_WU*U, sigma_W)
  Z <- rnorm(n, beta_Z0 + beta_ZA*A + beta_ZX*X + beta_ZL*L + beta_ZU*U, sigma_Z)
  Y <- rnorm(n, beta_Y0 + beta_YA*A + beta_YM*M + beta_YX*X + beta_YL*L + beta_YU*U, sigma_Y)
  data <- data.frame(X,A,U,L,M,Z,W,Y)
  return(data)
  
}

# # true values
# true_value <- function(a,a_star){
#   beta_L0 <- 0.2 ; beta_LA <- -0.25 ; beta_LX <- 0.3
#   beta_U0 <- 0.5 ; beta_UA <- -0.3 ; beta_UX <- 0.8 ; beta_UL <- 0.7
#   beta_M0 <- -0.6 ; beta_MA <- 0.7 ; beta_MX <- 1 ; beta_ML <- 0.25 ; beta_MU <- 0.5
#   beta_Y0 <- 0.2 ; beta_YA <- 0.5 ; beta_YX <- -0.6 ; beta_YL <- -0.3 ; beta_YU <- 1 ; beta_YM <- 0.8
#   data <- data_generate(sample(1:1000,1),10000000)
#   with(data, {
#     mu_La <- beta_L0 + beta_LA*a + beta_LX*X
#     mu_Ua <- beta_U0 + beta_UA*a + beta_UX*X + beta_UL*mu_La
#     model_M <- glm(M[A==a_star] ~ X[A==a_star], family = binomial)
#     mu_Mastar <- as.numeric(expit(cbind(1,X) %*% coef(model_M)))
#     mu_Ya_Mastar <- beta_Y0 + beta_YA*a + beta_YM*mu_Mastar + beta_YX*X + beta_YL*mu_La + beta_YU*mu_Ua
#     mean(mu_Ya_Mastar)
#   })
# }
# Y10_true <- true_value(1,0) # 1.2485

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
    model_A <- glm(A ~ sqrt(abs(X)), family = binomial)
    prob_A_hat <- as.numeric(expit(cbind(1,sqrt(abs(X))) %*% coef(model_A)))
  } else {
    model_A <- glm(A ~ X, family = binomial)
    prob_A_hat <- as.numeric(expit(cbind(1,X) %*% coef(model_A)))
  }
  pa_X <- prob_A_hat^a*(1-prob_A_hat)^(1-a)
  pastar_X <- prob_A_hat^a_star*(1-prob_A_hat)^(1-a_star)
  
  # P(W,L|a,X)
  if (mis_W) {
    model_WaX <- lm(W[A==a] ~ sqrt(abs(X))[A==a])
    muW_aX <- as.numeric(cbind(1,sqrt(abs(X))) %*% coef(model_WaX))
    model_LaX <- lm(L[A==a] ~ sqrt(abs(X))[A==a])
    muL_aX <- as.numeric(cbind(1,sqrt(abs(X))) %*% coef(model_LaX))
  } else {
    model_WaX <- lm(W[A==a] ~ X[A==a])
    muW_aX <- as.numeric(cbind(1,X) %*% coef(model_WaX))
    model_LaX <- lm(L[A==a] ~ X[A==a])
    muL_aX <- as.numeric(cbind(1,X) %*% coef(model_LaX))
  }
  
  # P(M|a_star,X)
  if (mis_M) {
    model_MastarX <- glm(M[A==a_star] ~ sqrt(abs(X))[A==a_star], family = binomial)
    probM_astarX <- as.numeric(expit(cbind(1,sqrt(abs(X))) %*% coef(model_MastarX)))
  } else {
    model_MastarX <- glm(M[A==a_star] ~ X[A==a_star], family = binomial)
    probM_astarX <- as.numeric(expit(cbind(1,X) %*% coef(model_MastarX)))
  }
  pM_astarX <- probM_astarX^M*(1-probM_astarX)^(1-M)
  
  # h_a(W,L,M,X)
  if (mis_h) {
    h_equation <- function(gamma) {
      h <- gamma[1] + gamma[2]*W + gamma[3]*L + gamma[4]*M + gamma[5]*sqrt(abs(X))
      f1 <- mean((Y-h)*A^a*(1-A)^(1-a))
      f2 <- mean((Y-h)*A^a*(1-A)^(1-a)*Z)
      f3 <- mean((Y-h)*A^a*(1-A)^(1-a)*L)
      f4 <- mean((Y-h)*A^a*(1-A)^(1-a)*M)
      f5 <- mean((Y-h)*A^a*(1-A)^(1-a)*sqrt(abs(X)))
      c(f1,f2,f3,f4,f5)
    }
    gamma_start <- c(0,0,0,0,0)
    gamma_hat <- nleqslv(gamma_start, h_equation)$x
    h_hat <- as.numeric(cbind(1,W,L,M,sqrt(abs(X))) %*% gamma_hat)
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
    q1_equation <- function(theta) {
      m <- 1
      indicator_M <- as.numeric(M == m)
      indicator_A <- as.numeric(A == a)
      q <- 1+exp((-1)^m*(theta[1]+theta[2]*Z+theta[3]*L+theta[4]*sqrt(abs(X))))
      f1 <- mean(indicator_A*(indicator_M*q-1)*1)
      f2 <- mean(indicator_A*(indicator_M*q-1)*W)
      f3 <- mean(indicator_A*(indicator_M*q-1)*L)
      f4 <- mean(indicator_A*(indicator_M*q-1)*sqrt(abs(X)))
      c(f1,f2,f3,f4)
    }
    theta1_start <- c(0,0,0,0)
    theta1_hat <- nleqslv(theta1_start, q1_equation)$x
    q0_equation <- function(theta) {
      m <- 0
      indicator_M <- as.numeric(M == m)
      indicator_A <- as.numeric(A == a)
      q <- 1+exp((-1)^m*(theta[1]+theta[2]*Z+theta[3]*L+theta[4]*sqrt(abs(X))))
      f1 <- mean(indicator_A*(indicator_M*q-1)*1)
      f2 <- mean(indicator_A*(indicator_M*q-1)*W)
      f3 <- mean(indicator_A*(indicator_M*q-1)*L)
      f4 <- mean(indicator_A*(indicator_M*q-1)*sqrt(abs(X)))
      c(f1,f2,f3,f4)
    }
    theta0_start <- c(0,0,0,0)
    theta0_hat <- nleqslv(theta0_start, q0_equation)$x
    q_hat <- ifelse(M == 1,
                    as.numeric(1 + exp((-1) * cbind(1,Z,L,sqrt(abs(X))) %*% theta1_hat)),
                    as.numeric(1 + exp(cbind(1,Z,L,sqrt(abs(X))) %*% theta0_hat)))
  } else {
    q1_equation <- function(theta) {
      m <- 1
      indicator_M <- as.numeric(M == m)
      indicator_A <- as.numeric(A == a)
      q <- 1+exp((-1)^m*(theta[1]+theta[2]*Z+theta[3]*L+theta[4]*X))
      f1 <- mean(indicator_A*(indicator_M*q-1)*1)
      f2 <- mean(indicator_A*(indicator_M*q-1)*W)
      f3 <- mean(indicator_A*(indicator_M*q-1)*L)
      f4 <- mean(indicator_A*(indicator_M*q-1)*X)
      c(f1,f2,f3,f4)
    }
    theta1_start <- c(0,0,0,0)
    theta1_hat <- nleqslv(theta1_start, q1_equation)$x
    q0_equation <- function(theta) {
      m <- 0
      indicator_M <- as.numeric(M == m)
      indicator_A <- as.numeric(A == a)
      q <- 1+exp((-1)^m*(theta[1]+theta[2]*Z+theta[3]*L+theta[4]*X))
      f1 <- mean(indicator_A*(indicator_M*q-1)*1)
      f2 <- mean(indicator_A*(indicator_M*q-1)*W)
      f3 <- mean(indicator_A*(indicator_M*q-1)*L)
      f4 <- mean(indicator_A*(indicator_M*q-1)*X)
      c(f1,f2,f3,f4)
    }
    theta0_start <- c(0,0,0,0)
    theta0_hat <- nleqslv(theta0_start, q0_equation)$x
    q_hat <- ifelse(M == 1,
                    as.numeric(1 + exp((-1) * cbind(1,Z,L,X) %*% theta1_hat)),
                    as.numeric(1 + exp(cbind(1,Z,L,X) %*% theta0_hat)))
  }
  # g0(W,X),g1(M,X),g(X)
  if (mis_h) {
    g0_hat <- as.numeric(cbind(1,W,L,probM_astarX,sqrt(abs(X))) %*% gamma_hat)
    g1_hat <- as.numeric(cbind(1,muW_aX,muL_aX,M,sqrt(abs(X))) %*% gamma_hat)
    g_hat <- as.numeric(cbind(1,muW_aX,muL_aX,probM_astarX,sqrt(abs(X))) %*% gamma_hat)
  } else {
    g0_hat <- as.numeric(cbind(1,W,L,probM_astarX,X) %*% gamma_hat)
    g1_hat <- as.numeric(cbind(1,muW_aX,muL_aX,M,X) %*% gamma_hat)
    g_hat <- as.numeric(cbind(1,muW_aX,muL_aX,probM_astarX,X) %*% gamma_hat)
  }
  # psi(a,a_star)
  indicator_a <- as.numeric(A == a)
  indicator_astar <- as.numeric(A == a_star)
  est_MR <- mean(indicator_a/pa_X*pM_astarX*q_hat*(Y-h_hat) 
                 + indicator_a/pa_X*(g0_hat-g_hat)
                 + indicator_astar/pastar_X*(g1_hat-g_hat) 
                 + g_hat)
  est_M1 <- mean(g_hat)
  est_M2 <- mean(indicator_astar/pastar_X*g1_hat)
  est_M3 <- mean(indicator_a/pa_X*g0_hat)
  est_M4 <- mean(indicator_a/pa_X*pM_astarX*q_hat*Y)
  est <- c(MR = est_MR, M1 = est_M1, M2 = est_M2, M3 = est_M3, M4 = est_M4)
  return(est)
}

bootstrap_parallel <- function(data, mis_A = FALSE, mis_W = FALSE, mis_M = FALSE, mis_h = FALSE, mis_q = FALSE) {
  
  results <- foreach(i = 1:n_boot, .combine = 'rbind', .packages = c("nleqslv")) %dopar% {
    indices <- sample(1:n, size = n, replace = TRUE)
    boot_sample <- data[indices, ]
    point_est(boot_sample, mis_A, mis_W, mis_M, mis_h, mis_q)
  }
  
  std <- apply(results,2,sd)
  est <- point_est(data, mis_A, mis_W, mis_M, mis_h, mis_q)
  bias <- est - true_value
  MSE <- std^2 + bias^2
  cover <- true_value <= est + 1.96 * std & true_value >= est - 1.96 * std
  
  list(estimates=est, std_errors=std, bias=bias, MSE=MSE, coverage=cover)
}

n <- 1000
n_sim <- 1000
n_boot <- 200
true_value <- 1.2485
a <- 1; a_star <- 0
n_cores <- detectCores()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c("data_generate", "expit", "point_est", "n", "true_value", "a", "a_star"))
settings <- list(
  all_correct = function() bootstrap_parallel(data),
  M1 = function() bootstrap_parallel(data, mis_A = 1, mis_q = 1),
  M2 = function() bootstrap_parallel(data, mis_M = 1, mis_q = 1),
  M3 = function() bootstrap_parallel(data, mis_W = 1, mis_q = 1),
  M4 = function() bootstrap_parallel(data, mis_W = 1, mis_h = 1)
)
results_list <- list()
for (setting_name in names(settings)) {
  results_list[[setting_name]] <- list(
    estimates = matrix(NA, nrow = n_sim, ncol = 5, 
                       dimnames = list(NULL, c("MR", "M1", "M2", "M3", "M4"))),
    std_errors = matrix(NA, nrow = n_sim, ncol = 5, 
                        dimnames = list(NULL, c("MR", "M1", "M2", "M3", "M4"))),
    bias = matrix(NA, nrow = n_sim, ncol = 5, 
                  dimnames = list(NULL, c("MR", "M1", "M2", "M3", "M4"))),
    MSE = matrix(NA, nrow = n_sim, ncol = 5, 
                 dimnames = list(NULL, c("MR", "M1", "M2", "M3", "M4"))),
    coverage = matrix(NA, nrow = n_sim, ncol = 5, 
                      dimnames = list(NULL, c("MR", "M1", "M2", "M3", "M4")))
  )
}
for (i in 1:n_sim) {
  data <- data_generate(i, n)
  for (setting_name in names(settings)) {
    res <- settings[[setting_name]]()
    results_list[[setting_name]]$estimates[i, ] <- res$estimates
    results_list[[setting_name]]$std_errors[i, ] <- res$std_errors
    results_list[[setting_name]]$bias[i, ] <- res$bias
    results_list[[setting_name]]$MSE[i, ] <- res$MSE
    results_list[[setting_name]]$coverage[i, ] <- res$coverage
  }
  if (i %% 5 == 0){
    print(paste("Completed simulation", i))
  }
}
stopCluster(cl)

summary_stats <- list()
for (setting_name in names(settings)) {
  summary_stats[[setting_name]] <- list(
    mean_estimates = colMeans(results_list[[setting_name]]$estimates, na.rm=TRUE),
    mean_std_errors = colMeans(results_list[[setting_name]]$std_errors, na.rm=TRUE),
    mean_MSE = colMeans(results_list[[setting_name]]$MSE, na.rm=TRUE),
    mean_bias = colMeans(results_list[[setting_name]]$bias, na.rm=TRUE),
    coverage = colMeans(results_list[[setting_name]]$coverage, na.rm=TRUE)
  )
}
time_end <- Sys.time()
time <- difftime(time_end,time_start)
print(time)
save.image(file = "simresults_binaryM_N1000_para.RData")
