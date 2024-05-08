#library(ismev)
library(xts)
library(extRemes)
library(lubridate)
library(nortest)
source("gpd.R")
#library(gnFit)
source("gnfit.R")
library(eva)
library(dgof)
source("projectedLevel.R")
library(numDeriv)

PoT_phi_function <- function(k) { # create a function with the name my_function
  
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
  

  
  # Bestämmer Peak
  u <- sort(data2$Value)[length(data2$Value)*0.98]
  declust <- c(decluster(data2$Value, u,replace.with=u-1))
  v=0
  for (i in declust){
    if (i > u){
      v <- v+1
    }
  }
  clusters <- v
  # För att kunna göra non-stationary behövs x värden
  x <- matrix(ncol = 1, nrow = length(declust))
  x[,1] <- seq(1,length(declust),1)
  
  fit.gpd2 <- gpd.fit_mod(data2$Value, threshold = u, ydat= x, sigl=TRUE, show = FALSE, method="BFGS", maxit = 30000)
  pre_lambda <- 365.25* (gpd.fit_mod(data2$Value, threshold = u, ydat= x, siglink = exp, sigl=TRUE, show = FALSE, method="BFGS", maxit = 30000, siginit = c(log(fit.gpd2$mle[1]), fit.gpd2$mle[2]/fit.gpd2$mle[1]), shinit = fit.gpd2$mle[3])$rate)
  
  #print(declust)
  # Gör en gpd
  #fit.gpd <- gpd.fit(declust, threshold = u, ydat= x, sigl=TRUE, show = FALSE, method="BFGS", maxit = 30000)
  #fit_mle <- fevd(Value, data2, method = "MLE", type="GP", threshold = u)
  fit.gpd2 <- gpd.fit_mod(declust, threshold = u, ydat= x, sigl=TRUE, show = FALSE, method="BFGS", maxit = 100000)
  fit.gpd <- (gpd.fit_mod(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, show = FALSE, maxit = 1000000, siginit = c(log(fit.gpd2$mle[1]), fit.gpd2$mle[2]/fit.gpd2$mle[1])))
  
  if (min(diag(fit.gpd$cov)) < 0){
 #   fit.gpd2 <- gpd.fit(declust, threshold = u, ydat= x, sigl=TRUE, show = FALSE, maxit = 100000)
    fit.gpd <- (gpd.fit_mod(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, show = FALSE, maxit = 100000, siginit = c(log(fit.gpd2$mle[1]), fit.gpd2$mle[2]/fit.gpd2$mle[1]) ))
  }
  
  if (min(diag(fit.gpd$cov)) < 0){
    #   fit.gpd2 <- gpd.fit(declust, threshold = u, ydat= x, sigl=TRUE, show = FALSE, maxit = 100000)
    fit.gpd <- (gpd.fit_mod(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, method="SANN", show = FALSE, maxit = 100000, siginit = c(log(fit.gpd2$mle[1]), fit.gpd2$mle[2]/fit.gpd2$mle[1])))
  }
  
  #print((fit.gpd$cov))
  #if (min(diag(fit.gpd$cov)) < 0){
    #   fit.gpd2 <- gpd.fit(declust, threshold = u, ydat= x, sigl=TRUE, show = FALSE, maxit = 100000)
  #  fit.gpd <- (gpd.fit(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, show = FALSE, maxit = 100000, siginit = c(log(fit.gpd2$mle[1]), 0), shinit = fit.gpd2$mle[3]))
  #}
  
  #if (min(diag(fit.gpd$cov)) < 0){
    #   fit.gpd2 <- gpd.fit(declust, threshold = u, ydat= x, sigl=TRUE, show = FALSE, maxit = 100000)
#    fit.gpd <- (gpd.fit(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, method="CG", show = FALSE, maxit = 100000, siginit = c(log(fit.gpd2$mle[1]), 0), shinit = fit.gpd2$mle[3]))
 # }
  
  
  #print(fit.gpd$mle)
  #print((fit.gpd$cov))
  mat <- fit.gpd$cov
  
  #jpeg(file=paste(i, ".jpeg"))
  #plot(fit_mle, "primary")
  #dev.off()
  
  # Gör en ML-skattning
  estimates <- fit.gpd$mle
  SError <- fit.gpd$se
  
  phi1_p <- 2*(1-pnorm(abs(estimates[2]/SError[2]), 0, 1)) 
  xi_p <- 2*(1-pnorm(abs(estimates[3]/SError[3]), 0, 1))
  
  #trend <- c()
  
  #for (x in 0:(length(data2$Value)-1)) {
  #  trend <- c(trend, data2$Value[x+1]/exp(estimates[2]*x))
  #}

  #print(fit.gpd$npy)
  
  #ad <- gnfit(trend, "gpd", pr = c(exp(estimates[1]), estimates[3]), threshold = u)$Apval
  #ad <- gpdAd(trend, bootstrap = TRUE, bootnum = 100)$p.value
  
  # Årligt max
  data2$Date <- data2$Date %m+% months(9)
  year_max <- apply.yearly(data2, max)
  
  #lambda <- fit.gpd$rate
  
  cov = matrix(0,4,4)
  
  cov[1,1] = mat[1,1]
  cov[2,2] = mat[3,3]
  cov[1,2] = cov[2,1] = mat[3,1]
  cov[1,3] = cov[3,1] = mat[2,1]*365.25
  cov[3,3] = mat[2,2]*365.25^2
  cov[2,3] = cov[3,2] = mat[3,2]*365.25
  
  #print(max(tail(year_max$Value, n=5)))
  #print(u)
  #print(lambda)
  #print(estimates[1]/SError[1])
  #print(estimates[2]/SError[2])
  #print(estimates[3]/SError[3])
  #print(cov)
  #print(year_max$Value)
  
  
  #(level, u, year, zeroYear, lambda, phi0, shape0, phi1=0, shape1=0, cov=matrix(0, 4, 4), samples=10000)
  
  yellow2=probAboveGP(level=max(tail(year_max$Value, n=5)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  orange2=probAboveGP(level=max(tail(year_max$Value, n=25)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  red2=probAboveGP(level=max(tail(year_max$Value, n=50)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  
  yellow=probAboveDayGP(level=max(tail(year_max$Value, n=5)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  orange=probAboveDayGP(level=max(tail(year_max$Value, n=25)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  red=probAboveDayGP(level=max(tail(year_max$Value, n=50)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[3], phi1=estimates[2]*365.25, cov=cov)
  
  #yellow=yellow2=orange=orange2=red=red2=0
  
  # För att göra datan non-stationary behövs x
  #x <- matrix(ncol = 1, nrow = length(year_max$Value))
  #x[,1] <- seq(1,length(year_max$Value),1)
  
  # Gör en GEV
  #fit.gev <- gev.fit(year_max$Value, ydat= x, siglink= exp, mul = TRUE, show = FALSE)
  #fit_mle <- fevd(Value, year_max, method = "MLE", type="GEV")
  #returnlevel <- return.level(fit_mle, return.period = c(5,25,50))
  
  # ML skattning
  #estimates2 <- fit.gev$mle
  
  #levels = (levelResultsGEV(50, 2024, 1960, estimates2[1], (estimates2[3]), estimates2[4], loc1=estimates2[2]))
  
  # fullösning
  #for (i in 1:length(SError)){
  #  if (is.nan(SError[i])){
  #    SError[i] <- 0
  #  }
  #}
  #print(estimates)
  #print(estimates/SError)
  #yellow=(probAboveGP(levels[14], u, 2024, 1960, lambda, estimates[1], estimates[3], phi1 = estimates[2], phi0SE=SError[1], shape0SE=SError[3], phi1SE = SError[2]))
  #orange=(probAboveGP(levels[25*3-2], u, 2024, 1960, lambda, estimates[1], estimates[3], phi1 = estimates[2], phi0SE=SError[1], shape0SE=SError[3], phi1SE = SError[2]))
  #red=(probAboveGP(levels[50*3-1], u, 2024, 1960, lambda, estimates[1], estimates[3], phi1 = estimates[2], phi0SE=SError[1], shape0SE=SError[3], phi1SE = SError[2]))
  
  
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
              #"post_datapts_pot2" = clusters,
              #"ad_pot2" = ad,
              #"pre_lambda_pot2" = pre_lambda,
              #"post_lambda_pot2" = lambda,
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


# Tar namnet på alla filer som finns i data mappen
list <- list.files("data")

# Loopar igenom alla filer
#for (i in list) {
#  PoT_phi_function(i)
#}

print(PoT_phi_function(list[1]))
#print(PoT_phi_function("RÄKTFORS_17.csv"))


"
lambda_pot2
low_lambda_pot2
high_lambda_pot2

yellow_tomorrow_pot2
orange_tomorrow_pot2
red_tomorrow_pot2

"