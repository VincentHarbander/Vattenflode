library(ismev)
#source("gev.R")
library(xts)
library(extRemes)
library(lubridate)
library(nortest)
#library(gnFit)
source("gnfit.R")
library(eva)
library(dgof)
source("ProjectedLevel.R")

BM_phi_function <- function(k) { # create a function with the name my_function
  
  #name
  result <- c()
  
  result <- c(result, k)
  
  # Läser datan från CSV
  path <- c("data/","")
  data <- read.csv(paste(path, collapse=k), sep = ";")
  
  # Skriver om datan till datum och värde form
  data2 <- data.frame(
    Date = as.Date(data$Datum),
    Value = data$Vatten
  )
  
  # Årligt max
  data2$Date <- data2$Date %m+% months(9)
  year_max <- apply.yearly(data2, max)
  
  
  # För att göra datan non-stationary behövs x
  x <- matrix(ncol = 1, nrow = length(year_max$Value))
  x[,1] <- seq(1,length(year_max$Value),1)
  
  # Gör en GEV
  #fit.gev <- gev.fit(year_max$Value, ydat= x, sigl = TRUE, siglink = exp, method="BFGS", show=FALSE, maxit=100000)
  #print(fit.gev$cov)
  fit.gev2 <- gev.fit(year_max$Value, ydat= x, sigl = TRUE, show=FALSE, method="BFGS", maxit=100000)
  fit.gev <- gev.fit(year_max$Value, ydat= x, sigl = TRUE, siglink = exp, show=FALSE, method="BFGS", maxit=100000, muinit = fit.gev2$mle[1], siginit = c(log(fit.gev2$mle[2]), fit.gev2$mle[3]/fit.gev2$mle[2]))#, shinit = fit.gev2$mle[4])
  
  if (min(diag(fit.gev$cov)) < 0){
    fit.gev <- gev.fit(year_max$Value, ydat= x, sigl = TRUE, siglink = exp, show=FALSE, maxit=100000, muinit = fit.gev2$mle[1], siginit = c(log(fit.gev2$mle[2]), fit.gev2$mle[3]/fit.gev2$mle[2]))#, shinit = fit.gev2$mle[4])
  }
  
  if (min(diag(fit.gev$cov)) < 0){
    fit.gev <- gev.fit(year_max$Value, ydat= x, sigl = TRUE, siglink = exp, show=FALSE, method="SANN", maxit=100000, muinit = fit.gev2$mle[1], siginit = c(log(fit.gev2$mle[2]), fit.gev2$mle[3]/fit.gev2$mle[2]))#, shinit = fit.gev2$mle[4])
  }
  #print(fit.gev$cov)
  #fit_mle <- fevd(Value, year_max, method = "MLE", type="GEV")
  #returnlevel <- return.level(fit_mle, return.period = c(5,25,50))
  
  mat <- (fit.gev$cov)
  #print(min(eigen(mat)$value))
  
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
  
  
  
  # ML skattning
  estimates <- fit.gev$mle
  SError <- fit.gev$se
  #print(estimates[1]/SError[1])
  #print(estimates[2]/SError[2]) 
  #print(estimates[3]/SError[3]) 
  #print(estimates[4]/SError[4])
  #print(estimates)
  #print(SError)
  #xi_test <- ifelse(0 > estimates[4]- 1.96*SError[4] & 0 < estimates[4]+ 1.96*SError[4], 0, estimates[4])
  #phi1_test <- ifelse(0 > (estimates[3])- 1.96*(SError[3]) & 0 < (estimates[3])+ 1.96*(SError[3]), 0, (estimates[3]))

  # p-värden
  #print("hej")
  
  phi1_p <- 2*(1-pnorm(abs(estimates[3]/SError[3]), 0, 1)) 
  #("då")
  xi_p <- 2*(1-pnorm(abs(estimates[4]/SError[4]), 0, 1)) 
  
  trend <- c()
  
  for (x in 0:(length(year_max$Value))-1) {
    trend <- c(trend, year_max$Value[x+1]/exp(estimates[3]*x))
  }
  
  ad <- gnfit(trend, "gev", pr = c(estimates[1], exp(estimates[2]), estimates[4]))$Apval
  
  # gul or
  
  
  yellow_30 = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2024+29, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  orange_30 = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2024+29, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  red_30 = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2024+29, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  
  yellow_now = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  orange_now = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  red_now = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[4], loc1=0, phi1=estimates[3], shape1=0, cov=cov, samples=10000)
  
  yellow_ratio = yellow_30/yellow_now - 1
  orange_ratio = orange_30/orange_now - 1
  red_ratio = red_30/red_now - 1
  
  
 # yellow_30 = orange_30 = red_30 = yellow_now = orange_now = red_now = yellow_ratio = orange_ratio =red_ratio = 0
  
  # fullösning
  #for (i in 1:length(SError)){
  #  if (is.nan(SError[i])){
  #    SError[i] <- 0
  #  }
  #}
  
  #levels = (levelResultsGEV(50, 2024, 1960, estimates[1], (estimates[2]), estimates[4], phi1=estimates[3]))

  #yellow=(probAboveGEV(levels[14], 30, 2024, 1960, estimates[1], (estimates[2]), estimates[4], phi1=estimates[3], loc0SE=SError[1], phi0SE=SError[2], shape0SE=SError[4], phi1SE=SError[3]))
  #orange= (probAboveGEV(levels[25*3-1], 30, 2024, 1960, estimates[1], (estimates[2]), estimates[4], phi1=estimates[3], loc0SE=SError[1], phi0SE=SError[2], shape0SE=SError[4], phi1SE=SError[3]))
  #red=(probAboveGEV(levels[50*3-1], 30, 2024, 1960, estimates[1], (estimates[2]), estimates[4], phi1=estimates[3], loc0SE=SError[1], phi0SE=SError[2], shape0SE=SError[4], phi1SE=SError[3]))
  

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
              #"phi1_tested_bm2" = phi1_test,
              
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


# Tar namnet på alla filer som finns i data mappen
list <- list.files("data")

# Loopar igenom alla filer
#for (i in list) {
#  BM_phi_function(i)
#}

(BM_phi_function("FLÖTEMARKEN_1375.csv"))

"


yellow_bm2
orange_bm2
red_bm2

"

