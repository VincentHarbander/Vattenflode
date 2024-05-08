library(ismev)
library(xts)
library(extRemes)
library(lubridate)
library(nortest)
#library(gnFit)
source("gnfit.R")
library(eva)
library(dgof)
source("projectedLevel.R")


BM_function <- function(k) { # create a function with the name my_function
  
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
  fit.gev2 <- gev.fit(year_max$Value, ydat= x, show = FALSE, method="BFGS")
  fit.gev <- gev.fit(year_max$Value, ydat= x, siglink = exp, show = FALSE, method="BFGS", muinit = fit.gev2$mle[1], siginit = log(fit.gev2$mle[2]), shinit = 0.1)
  #fit_mle <- fevd(Value, year_max, method = "MLE", type="GEV")
  #returnlevel <- return.level(fit_mle, return.period = c(5,25,50))
  
  mat <- (fit.gev$cov)
  #print(min(eigen(mat)$value))
  
  # ML skattning
  estimates <- fit.gev$mle
  SError <- fit.gev$se
  test <- ifelse(0 > estimates[3]- 1.96*SError[3] & 0 < estimates[3]+ 1.96*SError[3], 0, estimates[3])
  ad <- gnfit(year_max$Value, "gev", pr = c(estimates[1], exp(estimates[2]), estimates[3]))$Apval
  xi_p <- 2*(1-pnorm(abs(estimates[3]/SError[3]), 0, 1)) 
  
  cov <- matrix(0,6,6)
  cov[1,1] = mat[1,1]
  cov[2,1] = cov[1,2] = mat[2,1]
  cov[3,1] = cov[1,3] = mat[3,1]
  cov[2,2] = mat[2,2]
  cov[3,2] = cov[2,3] = mat[3,2]
  cov[3,3] = mat[3,3]
  
  yellow_now = probAboveGEV(max(tail(year_max$Value, n=5)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[3], cov=cov, samples=10000)
  orange_now = probAboveGEV(max(tail(year_max$Value, n=25)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[3], cov=cov, samples=10000)
  red_now = probAboveGEV(max(tail(year_max$Value, n=50)), 1, 2023, 1960, estimates[1], (estimates[2]), estimates[3], cov=cov, samples=10000)
  

  
  return(list("name"= k,
              #"loc_ratio_bm0" = estimates[1]/SError[1],
              #"phi_ratio_bm0" = estimates[2]/SError[2],
              #"xi_ratio_bm0" = estimates[3]/SError[3],
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


# Tar namnet på alla filer som finns i data mappen
list <- list.files("data")

# Loopar igenom alla filer
#for (i in list) {
#  my_function(i)
#}

#df5 = NULL

#for (i in list) {
#  df5 = rbind(df5, data.frame(BM_function(i)))
#}

#write.csv2(df5, "All data.csv", row.names=FALSE)

print(BM_function("RÖÅN_1378.csv"))

"
cov_11_bm0
cov_21_bm0
cov_31_bm0
cov_41_bm0
cov_51_bm0
cov_61_bm0
cov_22_bm0
cov_32_bm0
cov_42_bm0
cov_52_bm0
cov_62_bm0
cov_33_bm0
cov_43_bm0
cov_53_bm0
cov_63_bm0
cov_44_bm0
cov_54_bm0
cov_64_bm0
cov_55_bm0
cov_65_bm0
cov_66_bm0
"
