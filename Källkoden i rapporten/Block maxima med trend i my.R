library(ismev)
library(xts)
library(extRemes)
library(lubridate)
library(nortest)
#library(gnFit)
source("gnfit.R")
library(eva)
library(dgof)
source("ProjectedLevel.R")

BM_my_function <- function(k) { # create a function with the name my_function
  
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
  fit.gev <- gev.fit(year_max$Value, ydat= x, siglink= exp, mul = TRUE, show = FALSE, method="SANN")
  mat <- (fit.gev$cov)
  if (min(eigen(mat)$value) < 0){
    fit.gev <- gev.fit(year_max$Value, ydat= x, siglink= exp, mul = TRUE, show = FALSE, method="SANN", maxit = 1000000)    
  }
  #fit_mle <- fevd(Value, year_max, method = "MLE", type="GEV")
  #returnlevel <- return.level(fit_mle, return.period = c(5,25,50))
  
  mat <- (fit.gev$cov)
  #print(min(eigen(mat)$value))
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
  
  
  

  # ML skattning
  estimates <- fit.gev$mle
  #print(estimates)
  SError <- fit.gev$se
  xi_test <- ifelse(0 > estimates[4]- 1.96*SError[4] & 0 < estimates[4]+ 1.96*SError[4], 0, estimates[4])

  loc1_p <- 2*(1-pnorm(abs(estimates[2]/SError[2]), 0, 1)) 
  xi_p <- 2*(1-pnorm(abs(estimates[4]/SError[4]), 0, 1)) 

  
  trend <- c()
  
  for (x in 0:(length(year_max$Value)-1)) {
    trend <- c(trend, year_max$Value[x+1]-estimates[2]*(x))
  }
  
  ad <- gnfit(trend, "gev", pr = c(estimates[1], exp(estimates[3]), estimates[4]))$Apval
  
  
  yellow_30 = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2024+29, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  orange_30 = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2024+29, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  red_30 = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2024+29, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)

  yellow_now = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2023, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  orange_now = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2023, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  red_now = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2023, 1960, estimates[1], (estimates[3]), estimates[4], loc1=estimates[2], phi1=0, shape1=0, cov=cov, samples=10000)
  
  yellow_ratio = yellow_30/yellow_now - 1
  orange_ratio = orange_30/orange_now - 1
  red_ratio = red_30/red_now - 1
  
  #yellow_30 = orange_30 = red_30 = yellow_now = orange_now = red_now = yellow_ratio = orange_ratio =red_ratio = 0
  
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


# Tar namnet på alla filer som finns i data mappen
list <- list.files("data")

# Loopar igenom alla filer
#for (i in list) {
#  BM_my_function(i)
#}

print(BM_my_function(list[1]))
#print(BM_my_function("KÅTASELET_1480.csv"))

