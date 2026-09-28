library(haven)
library(dplyr)
library(nleqslv)
library(labelled)
library(nnet)
library(caret)
library(rdetools)
library(expm)
library(MASS)
library(ks)
library(mice)
library(doParallel)
library(foreach)
time_start <- Sys.time()

########## data process ##########
rand <- read_dta(
  "randhrs1992_2020v2.dta",
  col_select = c("hhid", "pn", "r15agey_e", "ragender", "raracem", "r15shlt", "raedegrm",
                 "rameduc", "rafeduc", "h15atotb", "h15atotn", "r15mstath", "r15lbsatwlf"))
rand2020 <- rand %>%
  filter(!is.na(r15agey_e)) %>%
  rename(HHID = hhid, PN = pn) # n = 15723
lb_data <- read_dta(
  "H20LB_R.dta",
  col_select = c("HHID", "PN",
                 "RLBELIG", "RLB029A", "RLB029B", "RLB029C", "RLB029D", "RLB029E", "RLB029F"))
lb_data <- lb_data %>%
  filter(RLBELIG==1) %>% # n = 7444
  rowwise() %>%
  mutate(discrim = mean(c(RLB029A,RLB029B,RLB029C,RLB029D,RLB029E,RLB029F), na.rm = TRUE)) %>%
  filter(!is.na(discrim)) %>% # n = 4626
  dplyr::select(HHID, PN, discrim)
merge_data <- left_join(lb_data,rand2020,by=c("HHID","PN"))
merge_data <- merge_data %>%
  dplyr::select(-c(1, 2)) %>%
  rename(
    age = r15agey_e,
    gender = ragender,
    race = raracem,
    marital = r15mstath,
    medu = rameduc,
    fedu = rafeduc,
    degree = raedegrm,
    health = r15shlt,
    networth = h15atotb,
    nonhousing = h15atotn,
    wellbeing = r15lbsatwlf
  )  %>%
  filter(race!=3) %>%
  filter(age>=50) %>%
  filter(!is.na(wellbeing)) # n = 4090
summary(merge_data)
merge_data <- zap_labels(merge_data)
merge_data$gender<-ifelse(merge_data$gender==2,0,1)
merge_data$marital<-ifelse(merge_data$marital<=2,0,1)
merge_data$race<-ifelse(merge_data$race==1,0,1)
imp_data <- mice(merge_data, method = "rf", m = 5, printFlag = FALSE)
merge_data<-complete(imp_data)
summary(merge_data)
saveRDS(merge_data, "wellbeing_data.rds")

appdata <- readRDS("wellbeing_data.rds")
appdata$discrim<-ifelse(appdata$discrim>5,0,1)
normalize <- function(x) {(x - min(x)) / (max(x) - min(x))}
appdata$networth <- normalize(appdata$networth)
appdata$nonhousing <- normalize(appdata$nonhousing)
appdata$medu <- normalize(appdata$medu)
appdata$fedu <- normalize(appdata$fedu)
appdata$age <- normalize(appdata$age)
appdata$degree <- normalize(appdata$degree)
appdata$health <- normalize(appdata$health)
summary(appdata)
choose_proxies <- function(data, proxy = 1:4) {
  proxies <- list(
    list(Z = "medu", W = "networth"),
    list(Z = "fedu", W = "networth"),
    list(Z = "medu", W = "nonhousing"),
    list(Z = "fedu", W = "nonhousing")
  )
  sel <- proxies[[proxy]]
  data %>%
    rename(
      X1 = age, X2 = gender, X3 = marital, L1 = degree, L2 = health,
      A = race, M = discrim, Y = wellbeing,
      Z = !!sym(sel$Z), W = !!sym(sel$W)
    )
}


########## parametric ##########
point_est_para <- function(data,a,a_star){

  expit <- function(x){exp(x)/{1 + exp(x)}}
  X <- with(data, cbind(X1,X2,X3))
  L <- with(data, cbind(L1,L2))
  A <- with(data, A)
  M <- with(data, M)
  Y <- with(data, Y)
  Z <- with(data, Z)
  W <- with(data, W)

  # P(A|X)
  model_A <- glm(A ~ X, family = binomial, data = data)
  prob_A_hat <- as.numeric(expit(predict(model_A)))
  pa_X <- prob_A_hat^a*(1-prob_A_hat)^(1-a)
  pastar_X <- prob_A_hat^a_star*(1-prob_A_hat)^(1-a_star)

  # P(W,L|a,X)
  model_WaX <- lm(W[A==a] ~ X[A==a,])
  mu_WaX <- as.numeric(cbind(1,X) %*% coef(model_WaX))
  model_LaX <- lm(L[A==a,] ~ X[A==a,])
  mu_L1aX <- as.numeric(cbind(1,X) %*% coef(model_LaX)[,1])
  mu_L2aX <- as.numeric(cbind(1,X) %*% coef(model_LaX)[,2])
  mu_LaX <- cbind(mu_L1aX, mu_L2aX)

  # P(M|a_star,X)
  model_Mastar <- glm(M[A==a_star] ~ X[A==a_star,], family = binomial)
  prob_MastarX <- as.numeric(expit(cbind(1,X) %*% coef(model_Mastar)))
  pM_astarX <- prob_MastarX^M*(1-prob_MastarX)^(1-M)

  # h_a(W,L,M,X)
  h_equation <- function(gamma) {
    h <- as.numeric(gamma[1] + gamma[2]*W + L %*% gamma[3:4] + gamma[5]*M + X %*% gamma[6:8])
    residual <- (Y-h)*A^a*(1-A)^(1-a)
    f1 <- mean(residual*1)
    f2 <- mean(residual*Z)
    f3 <- colMeans(residual*L)
    f4 <- mean(residual*M)
    f5 <- colMeans(residual*X)
    c(f1,f2,f3,f4,f5)
  }
  gamma_start <- rep(0,8)
  gamma_hat <- nleqslv(gamma_start, h_equation)$x
  h_hat <- as.numeric(cbind(1,W,L,M,X) %*% gamma_hat)

  # q_a(Z,L,M,X)
  q1_equation <- function(theta) {
    m <- 1
    indicator_M <- as.numeric(M == m)
    indicator_A <- as.numeric(A == a)
    q <- as.numeric(1 + exp((-1)^m*(theta[1] + theta[2]*Z + L %*% theta[3:4] + X %*% theta[5:7])))
    residual <- indicator_A*(indicator_M*q-1)
    f1 <- mean(residual*1)
    f2 <- mean(residual*W)
    f3 <- colMeans(residual*L)
    f4 <- colMeans(residual*X)
    c(f1,f2,f3,f4)
  }
  theta1_start <- rep(0,7)
  theta1_hat <- nleqslv(theta1_start, q1_equation)$x
  q0_equation <- function(theta) {
    m <- 0
    indicator_M <- as.numeric(M == m)
    indicator_A <- as.numeric(A == a)
    q <- as.numeric(1 + exp((-1)^m*(theta[1] + theta[2]*Z + L %*% theta[3:4] + X %*% theta[5:7])))
    residual <- indicator_A*(indicator_M*q-1)
    f1 <- mean(residual*1)
    f2 <- mean(residual*W)
    f3 <- colMeans(residual*L)
    f4 <- colMeans(residual*X)
    c(f1,f2,f3,f4)
  }
  theta0_start <- rep(0,7)
  theta0_hat <- nleqslv(theta0_start, q0_equation)$x
  q_hat <- ifelse(M == 1,
                  as.numeric(1 + exp((-1) * cbind(1,Z,L,X) %*% theta1_hat)),
                  as.numeric(1 + exp(cbind(1,Z,L,X) %*% theta0_hat)))

  # g0(W,X),g1(M,X),g(X)
  g0_hat <- as.numeric(cbind(1,W,L,prob_MastarX,X) %*% gamma_hat)
  g1_hat <- as.numeric(cbind(1,mu_WaX,mu_LaX,M,X) %*% gamma_hat)
  g_hat <- as.numeric(cbind(1,mu_WaX,mu_LaX,prob_MastarX,X) %*% gamma_hat)

  # psi(a,a_star)
  indicator_a <- as.numeric(A == a)
  indicator_astar <- as.numeric(A == a_star)
  est_MR <- mean(indicator_a/pa_X*pM_astarX*q_hat*(Y-h_hat)
                 + indicator_a/pa_X*(g0_hat-g_hat)
                 + indicator_astar/pastar_X*(g1_hat-g_hat)
                 + g_hat)
  return(est_MR)

}

########## bootstrap ##########
n <- nrow(appdata)
n_boot <- 200
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
n_cores <- detectCores()
cl <- makeCluster(getOption('cl.cores',n_cores))
print(paste("detectCores:",length(cl)))
clusterExport(cl, c("point_est_para", "n", "n_boot"))
registerDoParallel(cl)
options(warn = -1)

res_para <- data.frame()
for (proxy in 1:4) {

  cat("Starting caculation for proxy =", proxy, "\n")
  data <- choose_proxies(appdata, proxy)
  Y11_para_est <- point_est_para(data,1,1)
  Y10_para_est <-point_est_para(data,1,0)
  Y00_para_est <-point_est_para(data,0,0)
  IDE_para_est <- Y10_para_est - Y00_para_est
  IIE_para_est <- Y11_para_est - Y10_para_est

  res <- foreach(i = 1:n_boot, .combine = rbind, .packages = c("nleqslv")) %dopar% {
    tryCatch({
      set.seed(i)
      indices <- sample(1:n, size = n, replace = TRUE)
      boot_sample <- data[indices, ]
      Y11_para_i <- point_est_para(boot_sample,1,1)
      Y10_para_i <- point_est_para(boot_sample,1,0)
      Y00_para_i <- point_est_para(boot_sample,0,0)
      IDE_para_i <- Y10_para_i - Y00_para_i
      IIE_para_i <- Y11_para_i - Y10_para_i
      c(IDE_para_i=IDE_para_i, IIE_para_i=IIE_para_i)
    }, error = function(e) {
      c(IDE_para_i=NA, IIE_para_i=NA)
    })
  }
  IDE_para_se <- sd(res[,1])
  IIE_para_se <- sd(res[,2])
  IDE_para_lower <- IDE_para_est - 1.96 * IDE_para_se
  IDE_para_upper <- IDE_para_est + 1.96 * IDE_para_se
  IIE_para_lower <- IIE_para_est - 1.96 * IIE_para_se
  IIE_para_upper <- IIE_para_est + 1.96 * IIE_para_se
  
  res_para <- rbind(res_para,
                   data.frame(proxy = proxy,
                              effect = "IDE",
                              est = IDE_para_est,
                              se = IDE_para_se,
                              lower = IDE_para_lower,
                              upper = IDE_para_upper),
                   data.frame(proxy = proxy,
                              effect = "IIE",
                              est = IIE_para_est,
                              se = IIE_para_se,
                              lower = IIE_para_lower,
                              upper = IIE_para_upper))

}
stopCluster(cl)
print(res_para)


########## machine learning ##########
point_est_dml <- function(data,a,a_star){
  
  CF_K_fold <- 5
  n <- nrow(data)
  pseudo_Y <- numeric(n)
  kernel_sigma <- 15
  appro_rate <- 0.01
  lm_Hh <- 0.01
  lm_Qh <- 1
  lm_Hq <- 1
  lm_Qq <- 0.01
  
  random_sample <- createFolds(1:n, k = CF_K_fold)
  for (cf_rep in 1:CF_K_fold) {
    te_idx <- random_sample[[cf_rep]]
    tr_idx  <- setdiff(1:n, te_idx)
    tr_data <- data[-random_sample[[cf_rep]],]
    te_data <- data[random_sample[[cf_rep]],]
    n_tr <- n - length(random_sample[[cf_rep]])
    n_te <- length(random_sample[[cf_rep]])
    
    # P(a|X) & P(a_star|X)
    model_A <- nnet(factor(A) ~ X1+X2+X3, data = tr_data, size = 3, decay = 0.1, maxit = 500, trace = FALSE)
    prob_A_hat <- as.numeric(predict(model_A, newdata = data.frame(X1=te_data$X1, X2=te_data$X2, X3=te_data$X3), type = "raw"))
    pa_X <- prob_A_hat^a*(1-prob_A_hat)^(1-a)
    pastar_X <- prob_A_hat^a_star*(1-prob_A_hat)^(1-a_star)
    
    # P(W,L|a,X)
    W <- tr_data[tr_data$A == a,]$W
    L1 <- tr_data[tr_data$A == a,]$L1
    L2 <- tr_data[tr_data$A == a,]$L2
    X1 <- tr_data[tr_data$A == a,]$X1
    X2 <- tr_data[tr_data$A == a,]$X2
    X3 <- tr_data[tr_data$A == a,]$X3
    modell_WL <- nnet(x = matrix(c(X1,X2,X3), ncol=3), y = matrix(c(W,L1,L2), ncol=3), size = 3, decay = 0.1, maxit = 500, linout = TRUE, trace = FALSE)
    muWL_aX <- predict(modell_WL, newdata = matrix(c(te_data$X1,te_data$X2,te_data$X3), ncol = 3))
    resid_WL <- cbind(W,L1,L2) - predict(modell_WL, newdata = matrix(c(X1,X2,X3), ncol = 3))
    kde_WL <- kde(resid_WL)
    
    # P(M|a_star,X)
    model_M <- nnet(factor(M) ~ X1+X2+X3, data = tr_data[tr_data$A == a_star,], size = 3, decay = 0.1, maxit = 500, trace = FALSE)
    probM_astarX <- as.numeric(predict(model_M, newdata = data.frame(X1 = te_data$X1, X2 = te_data$X2, X3 = te_data$X3), type = "raw"))
    pM_astarX <- probM_astarX^te_data$M*(1-probM_astarX)^(1-te_data$M)
    
    # h_a(W,L,M,X) & q_a(Z,L,M,X)
    ## Gaussian RBF kernel ############################################################
    h_arg <- cbind(tr_data$W, tr_data$A, tr_data$L1, tr_data$L2, tr_data$X1, tr_data$X2, tr_data$X3)
    K_H <- rbfkernel(h_arg, sigma = kernel_sigma)
    q_arg <- cbind(tr_data$Z, tr_data$A, tr_data$L1, tr_data$L2, tr_data$X1, tr_data$X2, tr_data$X3)
    K_Q <- rbfkernel(q_arg, sigma = kernel_sigma)
    ###################################################################################
    
    ## Optimization with approximation#################################################
    r <- max(10, ceiling(n_tr * appro_rate))
    S <- diag(1, n_tr)[, sort(sample(1:n_tr, r))]
    ident_n <- diag(1, n_tr)
    ident_r <- diag(1, r)
    pM1 <- diag(as.numeric(tr_data$M),n_tr,n_tr)
    pM0 <- diag(as.numeric(1-tr_data$M),n_tr,n_tr)
    ## Optimization ###################################################################
    sqrtmtSK_QS <- sqrtm(ginv(t(S) %*% K_Q %*% S))
    D <- K_Q %*% S %*% sqrtmtSK_QS
    sqrtmtSK_HS <- sqrtm(ginv(t(S) %*% K_H %*% S))
    V <- K_H %*% S %*% sqrtmtSK_HS
    ## Optimization h ##
    Gamm <- D %*% ginv(t(D) %*% D/n_tr + lm_Qh/n_tr^0.8 * ident_r) %*% t(D)
    alpha_h1 <- ginv(t(V) %*% t(pM1) %*% Gamm %*% pM1 %*% V + n_tr^2*lm_Hh/n_tr^0.8 * ident_r) %*%
      t(V) %*% t(pM1) %*% Gamm %*% (tr_data$Y * tr_data$M)
    alpha_h0 <- ginv(t(V) %*% t(pM0) %*% Gamm %*% pM0 %*% V + n_tr^2*lm_Hh/n_tr^0.8 * ident_r) %*%
      t(V) %*% t(pM0) %*% Gamm %*% (tr_data$Y * (1-tr_data$M))
    ## Optimization q ##
    Gamm <- V %*% ginv(t(V) %*% V/n_tr + lm_Hq/n_tr^0.8 * ident_r) %*% t(V)
    alpha_q1 <- ginv(t(D) %*% t(pM1) %*% Gamm %*% pM1 %*% D + n_tr^2*lm_Qq/n_tr^0.8 * ident_r) %*%
      t(D) %*% t(pM1) %*% Gamm %*% rep(1, n_tr)
    alpha_q0 <- ginv(t(D) %*% t(pM0) %*% Gamm %*% pM0 %*% D + n_tr^2*lm_Qq/n_tr^0.8 * ident_r) %*%
      t(D) %*% t(pM0) %*% Gamm %*% rep(1, n_tr)
    ###################################################################################
    
    ## Estimate function h_hat ########################################################
    h_hat <- function(w,a,l1,l2,m,x1,x2,x3){
      W <- cbind(w,a,l1,l2,x1,x2,x3)
      K <- rbfkernel(h_arg, sigma = kernel_sigma, matrix(W, nrow(W), ncol(W)))
      return(m * t(K) %*% S %*% sqrtmtSK_HS %*% alpha_h1 + (1-m) * t(K) %*% S %*% sqrtmtSK_HS %*% alpha_h0)
    }
    ###################################################################################
    
    ## Estimate function q_hat ########################################################
    q_hat <- function(z,a,l1,l2,m,x1,x2,x3){
      Z <- cbind(z,a,l1,l2,x1,x2,x3)
      K <- rbfkernel(q_arg, sigma = kernel_sigma, matrix(Z, nrow(Z), ncol(Z)))
      return(m * t(K) %*% S %*% sqrtmtSK_QS %*% alpha_q1 + (1-m) * t(K) %*% S %*% sqrtmtSK_QS %*% alpha_q0)
    }
    ###################################################################################
    
    h_est <- h_hat(te_data$W,a,te_data$L1,te_data$L2,te_data$M,te_data$X1,te_data$X2,te_data$X3)
    q_est <- q_hat(te_data$Z,a,te_data$L1,te_data$L2,te_data$M,te_data$X1,te_data$X2,te_data$X3)
    
    # g0(W,X) & g1(M,X) & g(X)
    g0_est <- rep(0,n_te); g1_est <- rep(0,n_te); g_est <- rep(0,n_te)
    n_sample <- 100
    for (i in 1:n_te){
      sample_WL <- rkde(n_sample, kde_WL)
      sample_W <-  sample_WL[,1] + muWL_aX[i,1]
      sample_L1 <-  sample_WL[,2] + muWL_aX[i,2]
      sample_L2 <-  sample_WL[,3] + muWL_aX[i,3]
      sample_M <- rbinom(n_sample,1,probM_astarX[i])
      g0_est[i] <- mean(h_hat(te_data$W[i],a,te_data$L1[i],te_data$L2[i],sample_M,te_data$X1[i],te_data$X2[i],te_data$X3[i]), na.rm = TRUE)
      g1_est[i] <- mean(h_hat(sample_W,a,sample_L1,sample_L2,te_data$M[i],te_data$X1[i],te_data$X2[i],te_data$X3[i]), na.rm = TRUE)
      g_est[i] <- mean(h_hat(sample_W,a,sample_L1,sample_L2,sample_M,te_data$X1[i],te_data$X2[i],te_data$X3[i]), na.rm = TRUE)
    }
    
    # psi(a,a_star)
    indicator_a <- as.numeric(te_data$A == a)
    indicator_astar <- as.numeric(te_data$A == a_star)
    pseudo_Y[te_idx] <- as.numeric(indicator_a/pa_X*pM_astarX*q_est*(te_data$Y-h_est) 
                                   + indicator_a/pa_X*(g0_est-g_est)
                                   + indicator_astar/pastar_X*(g1_est-g_est) 
                                   + g_est)
  }
  
  return(pseudo_Y)
  
}

# seed <- sample(999999,1)
seed <- 498205
print(paste0("seed=",seed))
set.seed(seed)

res_dml <- data.frame()
for (proxy in 1:4) {
  cat("Starting caculation for proxy ", proxy, "\n")
  data <- choose_proxies(appdata, proxy)
  n <- nrow(data)
  
  pseudo_Y11_dml <- point_est_dml(data, 1, 1)
  pseudo_Y10_dml <- point_est_dml(data, 1, 0)
  pseudo_Y00_dml <- point_est_dml(data, 0, 0)
  
  IDE_dml_est <- mean(pseudo_Y10_dml - pseudo_Y00_dml)
  IIE_dml_est <- mean(pseudo_Y11_dml - pseudo_Y10_dml)
  
  IDE_dml_se <- sd(pseudo_Y10_dml - pseudo_Y00_dml) / sqrt(n)
  IIE_dml_se <- sd(pseudo_Y11_dml - pseudo_Y10_dml) / sqrt(n)
  
  IDE_dml_lower <- IDE_dml_est - 1.96 * IDE_dml_se
  IDE_dml_upper <- IDE_dml_est + 1.96 * IDE_dml_se
  IIE_dml_lower <- IIE_dml_est - 1.96 * IIE_dml_se
  IIE_dml_upper <- IIE_dml_est + 1.96 * IIE_dml_se
  
  res_dml <- rbind(res_dml,
                   data.frame(proxy = proxy,
                              effect = "IDE",
                              est = IDE_dml_est,
                              se = IDE_dml_se,
                              lower = IDE_dml_lower,
                              upper = IDE_dml_upper),
                   data.frame(proxy = proxy,
                              effect = "IIE",
                              est = IIE_dml_est,
                              se = IIE_dml_se,
                              lower = IIE_dml_lower,
                              upper = IIE_dml_upper))
}

print(res_dml)
time_end <- Sys.time()
time <- difftime(time_end,time_start)
print(time)
save.image(paste("app_results.RData",sep=''))
