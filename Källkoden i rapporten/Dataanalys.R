library(ismev)
library(extRemes)
library(gnFit)
source("projectedLevel.R")
library(NHPoisson)
library(lubridate)

# BM0 - Block maxima without trends
BM_function <- function(k, waterflow_data_transformad) {
  # Max per year
  year_max <- apply.yearly(waterflow_data_transformad, max)
  
  # The data needs to be in this form
  x <- matrix(ncol = 1, nrow = length(year_max$Value))
  x[,1] <- seq(1,length(year_max$Value),1)
  
  # ML estimates
  fit.gev_start_value <- gev.fit(year_max$Value, ydat= x, show = FALSE, method="BFGS")
  fit.gev <- gev.fit(year_max$Value, ydat= x, siglink = exp, show = FALSE, method="BFGS", muinit = fit.gev_start_value$mle[1], siginit = log(fit.gev_start_value$mle[2]), shinit = 0.1)
  
  estimates <- fit.gev$mle
  SError <- fit.gev$se

  # Covaraince matrix
  mat <- (fit.gev$cov)  
  cov <- matrix(0,6,6)
  cov[1,1] = mat[1,1]
  cov[2,1] = cov[1,2] = mat[2,1]
  cov[3,1] = cov[1,3] = mat[3,1]
  cov[2,2] = mat[2,2]
  cov[3,2] = cov[2,3] = mat[3,2]
  cov[3,3] = mat[3,3]
  
  # Tests
  ad <- gnfit(year_max$Value, "gev", pr = c(estimates[1], exp(estimates[2]), estimates[3]))$Apval
  xi_p <- 2*(1-pnorm(abs(estimates[3]/SError[3]), 0, 1)) 
  test <- ifelse(0 > estimates[3]- 1.96*SError[3] & 0 < estimates[3]+ 1.96*SError[3], 0, estimates[3])
  
  # Weather predictions 
  yellow_now = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[3], cov=cov, samples=10000)
  orange_now = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[3], cov=cov, samples=10000)
  red_now = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[3], cov=cov, samples=10000)
  
  return(list("name"= k,
              "loc_bm0"=  estimates[1], 
              "low_loc_bm0"=  estimates[1]- 1.96*SError[1],
              "high_loc_bm0"=  estimates[1]+ 1.96*SError[1],
              "loc_se_bm0" = SError[1],
              "phi_bm0" = estimates[2],
              "low_phi_bm0"=  estimates[2]- 1.96*SError[2],
              "high_phi_bm0"=  estimates[2]+ 1.96*SError[2],
              "phi_se_bm0" = SError[2],
              "xi_bm0" = estimates[3],
              "low_xi_bm0"=  estimates[3]- 1.96*SError[3],
              "high_xi_bm0"=  estimates[3]+ 1.96*SError[3],
              "xi_tested_bm0" = test,
              "xi_se_bm0" = SError[3],
              "xi_p_bm0" = xi_p,
              "ad_bm0" = ad,
              "yellow_now_bm0" = yellow_now,
              "orange_now_bm0" = orange_now,
              "red_now_bm0" = red_now,
              "cov_11_bm0" = mat[1,1],
              "cov_21_bm0" = mat[2,1],
              "cov_31_bm0" = mat[3,1],
              "cov_41_bm0" = 0,
              "cov_51_bm0" = 0,
              "cov_61_bm0" = 0,
              "cov_22_bm0" = mat[2,2],
              "cov_32_bm0" = mat[3,2],
              "cov_42_bm0" = 0,
              "cov_52_bm0" = 0,
              "cov_62_bm0" = 0,
              "cov_33_bm0" = mat[3,3],
              "cov_43_bm0" = 0,
              "cov_53_bm0" = 0,
              "cov_63_bm0" = 0,
              "cov_44_bm0" = 0,
              "cov_54_bm0" = 0,
              "cov_64_bm0" = 0,
              "cov_55_bm0" = 0,
              "cov_65_bm0" = 0,
              "cov_66_bm0" = 0,
              "yellow_level" = max(tail(year_max$Value, n=5)),
              "orange_level" = max(tail(year_max$Value, n=25)),
              "red_level"=max(tail(year_max$Value, n=50))
  ))
}

# BM1 - Block maxima with trends in my
BM_my_function <- function(k,waterflow_data_transformad) { 
  # Max per year
  year_max <- apply.yearly(waterflow_data_transformad, max)
  
  # The data needs to be in this form
  x <- matrix(ncol = 1, nrow = length(year_max$Value))
  x[,1] <- seq(1,length(year_max$Value),1)
  
  # ML estimate
  fit.gev <- gev.fit(year_max$Value, ydat= x, siglink= exp, mul = TRUE, show = FALSE, method="SANN")
  mat <- (fit.gev$cov)
  if (min(eigen(mat)$value) < 0){
    fit.gev <- gev.fit(year_max$Value, ydat= x, siglink= exp, mul = TRUE, show = FALSE, method="SANN", maxit = 1000000)    
  }
  
  estimates <- fit.gev$mle
  SError <- fit.gev$se
  
  # Covariance matrix
  mat <- (fit.gev$cov)
  cov <- matrix(0,6,6)
  cov[1,1] = mat[1,1]
  cov[1,2] = cov[2,1] = mat[3,1]
  cov[3,1] = cov[1,3] = mat[4,1]
  cov[1,4] = cov[4,1] = mat[2,1]
  cov[2,2] = mat[3,3]
  cov[2,3] = cov[3,2] = mat[4,3]
  cov[2,4] = cov[4,2] = mat[3,2]
  cov[3,3] = mat[4,4]
  cov[3,4] = cov[4,3] = mat[4,2]
  cov[4,4] = mat[2,2]
  
  # Tests
  xi_test <- ifelse(0 > estimates[4]- 1.96*SError[4] & 0 < estimates[4]+ 1.96*SError[4], 0, estimates[4])
  loc1_p <- 2*(1-pnorm(abs(estimates[2]/SError[2]), 0, 1)) 
  xi_p <- 2*(1-pnorm(abs(estimates[4]/SError[4]), 0, 1)) 
  
  # Weather predictions
  yellow_30 = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2024+29, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  orange_30 = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2024+29, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  red_30 = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2024+29, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  
  yellow_now = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2023, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  orange_now = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2023, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  red_now = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2023, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  
  yellow_ratio = yellow_30/yellow_now - 1
  orange_ratio = orange_30/orange_now - 1
  red_ratio = red_30/red_now - 1
  
  return(list("name"= k, 
              "loc0_bm1"=  estimates[1], 
              "low_loc0_bm1"=  estimates[1]- 1.96*SError[1],
              "high_loc0_bm1"=  estimates[1]+ 1.96*SError[1],
              "loc0_se_bm1" = SError[1],
              "loc1_bm1" = estimates[2],
              "low_loc1_bm1"=  estimates[2]- 1.96*SError[2],
              "high_loc1_bm1"=  estimates[2]+ 1.96*SError[2],
              "loc1_se_bm1" = SError[2],
              "phi_bm1" = estimates[3],
              "low_phi_bm1"=  estimates[3]- 1.96*SError[3],
              "high_phi_bm1"=  estimates[3]+ 1.96*SError[3],
              "phi_se_bm1" = SError[3],
              "xi_bm1" = estimates[4],
              "low_xi_bm1"=  estimates[4]- 1.96*SError[4],
              "high_xi_bm1"=  estimates[4]+ 1.96*SError[4],
              "xi_tested_bm1" = xi_test,
              "xi_se_bm1" = SError[4],
              #"ad_bm1" = ad, 
              "loc1_p_bm1" = loc1_p,
              "xi_p_bm1" = xi_p,
              "yellow_30_bm1" = yellow_30,
              "orange_30_bm1" = orange_30,
              "red_30_bm1" = red_30,
              "yellow_now_bm1" = yellow_now,
              "orange_now_bm1" = orange_now,
              "red_now_bm1" = red_now,
              "yellow_ratio_bm1" = yellow_ratio,
              "orange_ratio_bm1" = orange_ratio,
              "red_ratio_bm1" = red_ratio,
              "cov_11_bm1" = mat[1,1],
              "cov_21_bm1" = mat[3,1],
              "cov_31_bm1" = mat[4,1],
              "cov_41_bm1" = mat[2,1],
              "cov_51_bm1" = 0,
              "cov_61_bm1" = 0,
              "cov_22_bm1" = mat[3,3],
              "cov_32_bm1" = mat[4,3],
              "cov_42_bm1" = mat[3,2],
              "cov_52_bm1" = 0,
              "cov_62_bm1" = 0,
              "cov_33_bm1" = mat[4,4],
              "cov_43_bm1" = mat[4,2],
              "cov_53_bm1" = 0,
              "cov_63_bm1" = 0,
              "cov_44_bm1" = mat[2,2],
              "cov_54_bm1" = 0,
              "cov_64_bm1" = 0,
              "cov_55_bm1" = 0,
              "cov_65_bm1" = 0,
              "cov_66_bm1" = 0
  ))
}

# BM2 - Block maxima with trends in phi
BM_phi_function <- function(k,waterflow_data_transformad) { 
  # Max per year
  year_max <- apply.yearly(waterflow_data_transformad, max)
  
  # The data needs to be in this form
  x <- matrix(ncol = 1, nrow = length(year_max$Value))
  x[,1] <- seq(1,length(year_max$Value),1)
  
  # ML estimate
  fit.gev_start_value <- gev.fit(year_max$Value, ydat= x, sigl = TRUE, show=FALSE, method="BFGS", maxit=100000)
  fit.gev <- gev.fit(year_max$Value, ydat= x, sigl = TRUE, siglink = exp, show=FALSE, method="BFGS", maxit=100000, muinit = fit.gev_start_value$mle[1], siginit = c(log(fit.gev_start_value$mle[2]), fit.gev_start_value$mle[3]/fit.gev_start_value$mle[2]))#, shinit = fit.gev_start_value$mle[4])
  
  if (min(diag(fit.gev$cov)) < 0){
    fit.gev <- gev.fit(year_max$Value, ydat= x, sigl = TRUE, siglink = exp, show=FALSE, maxit=100000, muinit = fit.gev_start_value$mle[1], siginit = c(log(fit.gev_start_value$mle[2]), fit.gev_start_value$mle[3]/fit.gev_start_value$mle[2]))#, shinit = fit.gev_start_value$mle[4])
  }
  
  if (min(diag(fit.gev$cov)) < 0){
    fit.gev <- gev.fit(year_max$Value, ydat= x, sigl = TRUE, siglink = exp, show=FALSE, method="SANN", maxit=100000, muinit = fit.gev_start_value$mle[1], siginit = c(log(fit.gev_start_value$mle[2]), fit.gev_start_value$mle[3]/fit.gev_start_value$mle[2]))#, shinit = fit.gev_start_value$mle[4])
  }
  
  estimates <- fit.gev$mle
  SError <- fit.gev$se
  
  # Covariance matrix
  mat <- (fit.gev$cov)
  cov <- matrix(0,6,6)
  cov[1,1] = mat[1,1]
  cov[1,2] = cov[2,1] = mat[2,1]
  cov[3,1] = cov[1,3] = mat[4,1]
  cov[1,5] = cov[5,1] = mat[3,1]
  cov[2,2] = mat[2,2]
  cov[2,3] = cov[3,2] = mat[4,2]
  cov[2,5] = cov[5,2] = mat[3,2]
  cov[3,3] = mat[4,4]
  cov[3,5] = cov[5,3] = mat[4,3]
  cov[5,5] = mat[3,3]
  
  # Tests
  phi1_p <- 2*(1-pnorm(abs(estimates[3]/SError[3]), 0, 1)) 
  xi_p <- 2*(1-pnorm(abs(estimates[4]/SError[4]), 0, 1)) 
  
  # Weather predictions 
  yellow_30 = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2024+29, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  orange_30 = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2024+29, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  red_30 = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2024+29, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  
  yellow_now = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  orange_now = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  red_now = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  
  yellow_ratio = yellow_30/yellow_now - 1
  orange_ratio = orange_30/orange_now - 1
  red_ratio = red_30/red_now - 1
  
  return(list("name"= k, 
              "loc_bm2"=  estimates[1], 
              "low_loc_bm2"=  estimates[1]- 1.96*SError[1],
              "high_loc_bm2"=  estimates[1]+ 1.96*SError[1],
              "phi0_bm2" = (estimates[2]),
              "low_phi0_bm2"=  (estimates[2])- 1.96*(SError[2]),
              "high_phi0_bm2"=  (estimates[2])+ 1.96*(SError[2]),
              "phi1_bm2" = (estimates[3]),
              "low_phi1_bm2"=  (estimates[3])- 1.96*(SError[3]),
              "high_phi1_bm2"=  (estimates[3])+ 1.96*(SError[3]),
              "xi_bm2" = estimates[4],
              "low_xi_bm2"=  estimates[4]- 1.96*SError[4],
              "high_xi_bm2"=  estimates[4]+ 1.96*SError[4],
              #"xi_tested_bm2" = xi_test,
              #"ad_bm2" = ad,
              "phi1_p_bm2" = phi1_p,
              "xi_p_bm2" = xi_p,
              "loc_se_bm2" = SError[1],
              "phi0_se_bm2" = SError[2], 
              "phi1_se_bm2" = SError[3],
              "xi_se_bm2" = SError[4],
              "yellow_30_bm2" = yellow_30,
              "orange_30_bm2" = orange_30,
              "red_30_bm2" = red_30,
              "yellow_now_bm2" = yellow_now,
              "orange_now_bm2" = orange_now,
              "red_now_bm2" = red_now,
              "yellow_ratio_bm2" = yellow_ratio,
              "orange_ratio_bm2" = orange_ratio,
              "red_ratio_bm2" = red_ratio,
              "cov_11_bm2" = mat[1,1],
              "cov_21_bm2" = mat[2,1],
              "cov_31_bm2" = mat[4,1],
              "cov_41_bm2" = 0,
              "cov_51_bm2" = mat[3,1],
              "cov_61_bm2" = 0,
              "cov_22_bm2" = mat[2,2],
              "cov_32_bm2" = mat[4,2],
              "cov_42_bm2" = 0,
              "cov_52_bm2" = mat[3,2],
              "cov_62_bm2" = 0,
              "cov_33_bm2" = mat[4,4],
              "cov_43_bm2" = 0,
              "cov_53_bm2" = mat[4,3],
              "cov_63_bm2" = 0,
              "cov_44_bm2" = 0,
              "cov_54_bm2" = 0,
              "cov_64_bm2" = 0,
              "cov_55_bm2" = mat[3,3],
              "cov_65_bm2" = 0,
              "cov_66_bm2" = 0
              
  ))
}

# PoT0 - Peaks over thereshold without trends
PoT_function <- function(k,waterflow_data_transformad) {
  # Max per year
  year_max <- apply.yearly(waterflow_data_transformad, max)
  
  # The data needs to be in this form
  x <- matrix(ncol = 1, nrow = length(waterflow_data_transformad$Value))
  x[,1] <- seq(1,length(waterflow_data_transformad$Value),1)
  
  # Decide peak
  u <- sort(waterflow_data_transformad$Value)[length(waterflow_data_transformad$Value)*0.98]
  declust <- c(decluster(waterflow_data_transformad$Value, u,replace.with=u-1))
  
  # Count clusters
  v=0
  for (i in declust){
    if (i > u){
      v <- v+1
    }
  }
  clusters <- v
  
  pre_lambda <- 365.25*(gpd.fit(waterflow_data_transformad$Value, threshold = u, ydat= x, siglink = exp, show = FALSE, method="BFGS")$rate)
  
  # ML estimate
  fit.gpd_start_value <- gpd.fit(declust, threshold = u, ydat= x, show = FALSE, method="BFGS", maxit = 100000)
  fit.gpd <- gpd.fit(declust, threshold = u, ydat= x, siglink = exp, show = FALSE, method="BFGS", maxit = 100000, siginit = log(fit.gpd_start_value$mle[1]))
  
  estimates <- fit.gpd$mle
  SError <- fit.gpd$se
  
  # Covariance matrix
  mat <- fit.gpd$cov  
  cov = matrix(0,4,4)
  cov[1,1] = mat[1,1]
  cov[2,2] = mat[2,2]
  cov[1,2] = cov[2,1] = mat[2,1]
  
  # Make a vector with time between clusters
  time <- 0
  time_between_clusters <- c()
  for (i in declust){
    if (i > u){
      time_between_clusters <- c(time_between_clusters, time)
    }
    time <- time +1
  }
  
  # Estimate lambdas
  lambda <- 365.25*fit.gpd_start_value$rate  
  all.years <-seq(1,length(declust),length.out=length(declust))
  out<-fitPP.fun(covariates=cbind(all.years), posE=time_between_clusters, start=list(b0=0,b1=0), modSim = TRUE)
  omegaSE <- (sqrt(diag(out@vcov)))
  
  # Tests
  ad <- gnfit(waterflow_data_transformad$Value, "gpd", pr = c(exp(estimates[1]), estimates[2]), threshold = u)$Apval
  xi_p <- 2*(1-pnorm(abs(estimates[2]/SError[2]), 0, 1))
  omega_error <- 2*(1-pnorm(abs((out@detailsb$par[2]/omegaSE[2])), 0, 1))

  # Weather predictions 
  yellow2=probAboveGP(level=max(tail(year_max$Value, n=5)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  orange2=probAboveGP(level=max(tail(year_max$Value, n=25)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  red2=probAboveGP(level=max(tail(year_max$Value, n=50)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  
  yellow = probAboveDayGP(level=max(tail(year_max$Value, n=5)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  orange=probAboveDayGP(level=max(tail(year_max$Value, n=25)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  red=probAboveDayGP(level=max(tail(year_max$Value, n=50)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)

  return(list("name"= k,
              "pre_lambda" = pre_lambda,
              "post_lambda" = lambda,
              "threshold_u" = u,
              "phi_pot0"=  estimates[1], 
              "low_phi_pot0"=  estimates[1]- 1.96*SError[1],
              "high_phi_pot0"=  estimates[1]+ 1.96*SError[1],
              "phi_se_pot0" = SError[1],
              "xi_pot0" = estimates[2],
              "low_xi_pot0"=  estimates[2]- 1.96*SError[2],
              "high_xi_pot0"=  estimates[2]+ 1.96*SError[2],
              "xi_se_pot0" = SError[2],
              "xi_p_pot0" = xi_p,
              "ad_pot0" = ad,
              "lambda0" = out@detailsb$par[1] + log(365.25),
              "lambda1" = out@detailsb$par[2]*365.25,
              "lambda_cov_11" = out@vcov[1,1],
              "lambda_cov_22" = out@vcov[2,2]*365.25^2,
              "lambda_cov_21" = out@vcov[1,2]*365.25,
              "lambda0_se" = omegaSE[1],
              "lambda1_se" = omegaSE[2]*365.25,
              "lambda1_p" = omega_error,
              "post_datapts" = clusters,
              "yellow_tomorrow_pot0" = yellow,
              "orange_tomorrow_pot0" = orange,
              "red_tomorrow_pot0" = red,
              "yellow_next_pot0" = yellow2,
              "orange_next_pot0" = orange2,
              "red_next_pot0" = red2,
              "cov_11_pot0" = mat[1,1],
              "cov_21_pot0" = mat[2,1],
              "cov_31_pot0" = 0,
              "cov_41_pot0" = 0,
              "cov_22_pot0" = mat[2,2],
              "cov_32_pot0" = 0,
              "cov_42_pot0" = 0,
              "cov_33_pot0" = 0,
              "cov_43_pot0" = 0,
              "cov_44_pot0" = 0
              
  ))
}

# PoT2 - Peaks over thereshold with trends in phi
PoT_phi_function <- function(k,waterflow_data_transformad) {
  # max per year
  year_max <- apply.yearly(waterflow_data_transformad, max)
  
  # decide threshold
  u <- sort(waterflow_data_transformad$Value)[length(waterflow_data_transformad$Value)*0.98]
  declust <- c(decluster(waterflow_data_transformad$Value, u,replace.with=u-1))
  
  # count clusters
  v=0
  for (i in declust){
    if (i > u){
      v <- v+1
    }
  }
  clusters <- v
  
  # The data needs to be in this form
  x <- matrix(ncol = 1, nrow = length(declust))
  x[,1] <- seq(1,length(declust),1)
  
  # ML estimate, uses a modded version of gpd.fit that uses numderive instead of optim for hessian, the modded gpd.fit_mod is not included and should be changed to gpd.fit
  fit.gpd_start_value <- gpd.fit_mod(waterflow_data_transformad$Value, threshold = u, ydat= x, sigl=TRUE, show = FALSE, method="BFGS", maxit = 30000)
  pre_lambda <- 365.25* (gpd.fit_mod(waterflow_data_transformad$Value, threshold = u, ydat= x, siglink = exp, sigl=TRUE, show = FALSE, method="BFGS", maxit = 30000, siginit = c(log(fit.gpd_start_value$mle[1]), fit.gpd_start_value$mle[2]/fit.gpd_start_value$mle[1]), shinit = fit.gpd_start_value$mle[3])$rate)
  
  fit.gpd_start_value <- gpd.fit_mod(declust, threshold = u, ydat= x, sigl=TRUE, show = FALSE, method="BFGS", maxit = 100000)
  fit.gpd <- (gpd.fit_mod(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, show = FALSE, maxit = 1000000, siginit = c(log(fit.gpd_start_value$mle[1]), fit.gpd_start_value$mle[2]/fit.gpd_start_value$mle[1])))
  
  if (min(diag(fit.gpd$cov)) < 0){
    fit.gpd <- (gpd.fit_mod(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, show = FALSE, maxit = 100000, siginit = c(log(fit.gpd_start_value$mle[1]), fit.gpd_start_value$mle[2]/fit.gpd_start_value$mle[1]) ))
  }
  
  if (min(diag(fit.gpd$cov)) < 0){
    fit.gpd <- (gpd.fit_mod(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, method="SANN", show = FALSE, maxit = 100000, siginit = c(log(fit.gpd_start_value$mle[1]), fit.gpd_start_value$mle[2]/fit.gpd_start_value$mle[1])))
  }
  
  estimates <- fit.gpd$mle
  SError <- fit.gpd$se

  # covariance matrix
  mat <- fit.gpd$cov  
  cov = matrix(0,4,4)
  cov[1,1] = mat[1,1]
  cov[2,2] = mat[3,3]
  cov[1,2] = cov[2,1] = mat[3,1]
  cov[1,3] = cov[3,1] = mat[2,1]*365.25
  cov[3,3] = mat[2,2]*365.25^2
  cov[2,3] = cov[3,2] = mat[3,2]*365.25
  
  # tests
  phi1_p <- 2*(1-pnorm(abs(estimates[2]/SError[2]), 0, 1)) 
  xi_p <- 2*(1-pnorm(abs(estimates[3]/SError[3]), 0, 1))
  
  # weather predictions 
  yellow2=probAboveGP(level=max(tail(year_max$Value, n=5)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  orange2=probAboveGP(level=max(tail(year_max$Value, n=25)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  red2=probAboveGP(level=max(tail(year_max$Value, n=50)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  
  yellow=probAboveDayGP(level=max(tail(year_max$Value, n=5)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  orange=probAboveDayGP(level=max(tail(year_max$Value, n=25)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  red=probAboveDayGP(level=max(tail(year_max$Value, n=50)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  
  return(list("name"= k, 
              "phi0_pot2"=  estimates[1], 
              "low_phi0_pot2"=  estimates[1]- 1.96*SError[1],
              "high_phi0_pot2"=  estimates[1]+ 1.96*SError[1],
              "phi0_se_pot2" = SError[1],
              "phi1_pot2"=  estimates[2]*365.25, 
              "low_phi1_pot2"=  365.25*(estimates[2]- 1.96*SError[2]),
              "high_phi1_pot2"=  365.25*(estimates[2]+ 1.96*SError[2]),
              "phi1_se_pot2" = SError[2]*365.25,
              "xi_pot2" = estimates[3],
              "low_xi_pot2"=  estimates[3]- 1.96*SError[3],
              "high_xi_pot2"=  estimates[3]+ 1.96*SError[3],
              "xi_se_pot2" = SError[3],
              "xi_p_pot2" = xi_p,
              "yellow_tomorrow_pot2" =  yellow,
              "orange_tomorrow_pot2" = orange,
              "red_tomorrow_pot2" = red,
              "yellow_next_pot2" =  yellow2,
              "orange_next_pot2" = orange2,
              "red_next_pot2" = red2,
              "phi1_p_pot2" = phi1_p,
              "cov_11_pot2" = cov[1,1],
              "cov_21_pot2" = cov[2,1],
              "cov_31_pot2" = cov[3,1],
              "cov_41_pot2" = 0,
              "cov_22_pot2" = cov[2,2],
              "cov_32_pot2" = cov[3,2],
              "cov_42_pot2" = 0,
              "cov_33_pot2" = cov[2,2],
              "cov_43_pot2" = 0,
              "cov_44_pot2" = 0
  ))
}

# Function for Benjamini-Hochberg
BH_function <- function(data, column, k){
  BH <- list((data["name"]), p.adjust(t(data[column]), method = "BH"))
  names(BH) <- c("name", k)
  data = merge(data, data.frame(BH), by="name")
  return(data)
} 

# Tests if the p-values are significant
tester <- function(data, column, replacer, k){
  c <-c()
  for (i in 1:length(data[,column])){
    if (is.nan(data[column][i,1])){
      data[column][i,1] <- 1 
    }
    if (data[column][i,1] < 0.05){
      c <-c(c, data[replacer][i,1])
    }
    else{
      c <- c(c, 0)
    }
  }
  ans <- list(data["name"], c)
  names(ans) <- c("name", k)
  return(merge(data, data.frame(ans), by = "name"))
}

list <- list.files("data")

bm0 = NULL
bm1 = NULL
bm2 = NULL
pot0 = NULL
pot2 = NULL

# Loops through all files
for (i in list) {
  # Read data from CSV
  path <- c("data/","")
  data <- read.csv(paste(path, collapse=i), sep = ";")
  
  # Rewrites the data into date and value form
  waterflow_data_transformad <- data.frame(
    Date = as.Date(data$Datum),
    Value = data$Vatten
  )

  waterflow_data_transformad$Date <- waterflow_data_transformad$Date %m+% months(9)
  
  # Run the models
  bm0 = rbind(bm0, data.frame(BM_function(i, waterflow_data_transformad)))
  bm1 = rbind(bm1, data.frame(BM_my_function(i, waterflow_data_transformad)))
  bm2 = rbind(bm2, data.frame(BM_phi_function(i, waterflow_data_transformad)))
  pot0 = rbind(pot0, data.frame(PoT_function(i, waterflow_data_transformad)))
  pot2 = rbind(pot2, data.frame(PoT_phi_function(i, waterflow_data_transformad)))
}

# Benjamini-Hochberg for Anderson-Darling and null tests
bm0 = BH_function(bm0, "ad_bm0", "adbh_bm0")
bm0 = BH_function(bm0, "xi_p_bm0", "xi_bh_bm0")

bm1 = BH_function(bm1, "loc1_p_bm1", "loc1_bh_bm1")
bm1 = BH_function(bm1, "xi_p_bm1", "xi_bh_bm1")

bm2 = BH_function(bm2, "phi1_p_bm2", "phi1_bh_bm2")
bm2 = BH_function(bm2, "xi_p_bm2", "xi_bh_bm2")

pot0 = BH_function(pot0, "ad_pot0", "adbh_pot0")
pot0 = BH_function(pot0, "xi_p_pot0", "xi_bh_pot0")
pot0 = BH_function(pot0, "lambda1_p", "lambda1_bh")

pot2 = BH_function(pot2, "phi1_p_pot2", "phi1_bh_pot2")
pot2 = BH_function(pot2, "xi_p_pot2", "xi_bh_pot2")

# Checks if the trend parameters are zero
bm1 = tester(bm1, "loc1_bh_bm1","loc1_bm1","loc1_tested_bm1")
bm2 = tester(bm2, "phi1_bh_bm2","phi1_bm2","phi1_tested_bm2")

# Merge all models
Merged_data <- merge(bm0,bm1,by="name")
Merged_data <- merge(Merged_data,bm2,by="name")
Merged_data <- merge(Merged_data,pot0,by="name")
Merged_data <- merge(Merged_data,pot2,by="name")

# Write the csv
write.csv2(Merged_data, "All data.csv", row.names=FALSE)
