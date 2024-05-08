library(ismev)
#source("gpd.R")
library(xts)
library(extRemes)
library(lubridate)
library(nortest)
#library(gnFit)
source("gnfit.R")
library(eva)
library(dgof)
source("projectedLevel.R")
library(NHPoisson)
library(dplyr)

PoT_function <- function(k) { # create a function with the name my_function
  
  #name
  result <- c()
  
  result <- c(result, k)
  
  # Läser datan från CSV
  path <- c("Data/","")
  data <- read.csv(paste(path, collapse=k), sep = ";")
  
  # Skriver om datan till datum och värde form
  data2 <- data.frame(
    Date = as.Date(data$Datum),
    Value = data$Vatten
  )
  data2$Date <- data2$Date %m+% months(9)
  
  # För att kunna göra non-stationary behövs x värden
  x <- matrix(ncol = 1, nrow = length(data2$Value))
  x[,1] <- seq(1,length(data2$Value),1)
  
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
  
  
  pre_lambda <- 365.25*(gpd.fit(data2$Value, threshold = u, ydat= x, siglink = exp, show = FALSE, method="BFGS")$rate)
  
  # Gör en gpd
  fit.gpd2 <- gpd.fit(declust, threshold = u, ydat= x, show = FALSE, method="BFGS", maxit = 100000)
  fit.gpd <- gpd.fit(declust, threshold = u, ydat= x, siglink = exp, show = FALSE, method="BFGS", maxit = 100000, siginit = log(fit.gpd2$mle[1]))
  
  #if (min(diag(fit.gpd$cov))<0){
  #  fit.gpd <- (gpd.fit(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, show = FALSE, maxit = 100000, siginit = c(log(fit.gpd2$mle[1]), fit.gpd2$mle[2]/fit.gpd2$mle[1]), shinit = fit.gpd2$mle[3]))
  #}
  
  #if (min(diag(fit.gpd$cov))<0){
  #  fit.gpd <- (gpd.fit(declust, threshold = u, ydat= x, siglink = exp, sigl=TRUE, show = FALSE, method="SANN", maxit = 100000, siginit = c(log(fit.gpd2$mle[1]), fit.gpd2$mle[2]/fit.gpd2$mle[1]), shinit = fit.gpd2$mle[3]))
  #}
  #fit_mle <- fevd(Value, data2, method = "MLE", type="GP", threshold = u)
  mat <- fit.gpd$cov
  # ^ changed this to mat, used to be called "matrix"
  #print(fit.gpd$cov)
  #jpeg(file=paste(i, ".jpeg"))
  #plot(fit_mle, "primary")
  #dev.off()
  
  # Gör en ML-skattning
  estimates <- fit.gpd$mle
  SError <- fit.gpd$se
  
  
  #print(data2$Value)
  # gnfit skiter sig!!
  ad <- gnfit(data2$Value, "gpd", pr = c(exp(estimates[1]), estimates[2]), threshold = u)$Apval
  xi_p <- 2*(1-pnorm(abs(estimates[2]/SError[2]), 0, 1))


  
  
  
  
  cov = matrix(0,4,4)
  
  cov[1,1] = mat[1,1]
  cov[2,2] = mat[2,2]
  cov[1,2] = cov[2,1] = mat[2,1]

  # Årligt max
  year_max <- apply.yearly(data2, max)
  lambda <- 365.25*fit.gpd2$rate
  
 # data3 %>%
#    group_by(unique(year(data2$Date))) %>%
 #   summarize(n_break_threshold=sum(data2$Value > u))
  
  
  dat = data.frame(date = data2$Date,
                   year = year(data2$Date),
                   value = declust)
  #daa <- dat %>%
  #  group_by(year) %>%
  #  filter(value > u) %>%
  #  summarize(exceedances = n())
  #print((daa$exceedances))
  s <- 0
  ex <- c()
  for (i in declust){
    if (i > u){
      ex <- c(ex, s)
    #  k <- k -1
    }
    s <- s +1
  }
  
  #ex <- (POTevents.fun(declust, u, date = data2$Value)$L)
  
  all.seconds<-seq(1,length(declust),length.out=length(declust))
  out<-fitPP.fun(covariates=cbind(all.seconds), posE=ex, start=list(b0=0,b1=0), modSim = TRUE) # b0=0 is our guess at initial value for optimization, which is internally made 
  #out<-fitPP.fun(n =length(declust), posE=ex, start=list(b0=0), modSim = TRUE) # b0=0 is our guess at initial value for optimization, which is internally made 
  s <- (sqrt(diag(out@vcov)))
  
  omega_error <- 2*(1-pnorm(abs((out@detailsb$par[2]/s[2])), 0, 1))
  #print(out@detailsb$par[1]/s[1])
  #print(out@detailsb$par[2]/s[2])
  #print(out)

  #print(max(tail(year_max$Value, n=5)))
  #print(u)
  #print(lambda)
  #print(estimates[1])
  #print(estimates[2])
  #print(cov)
  #print(year_max$Value)

  
  #(level, u, year, zeroYear, lambda, phi0, shape0, phi1=0, shape1=0, cov=matrix(0, 4, 4), samples=10000)
  
  yellow2=probAboveGP(level=max(tail(year_max$Value, n=5)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  orange2=probAboveGP(level=max(tail(year_max$Value, n=25)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  red2=probAboveGP(level=max(tail(year_max$Value, n=50)), u=u, year=2025, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  
  yellow = probAboveDayGP(level=max(tail(year_max$Value, n=5)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  orange=probAboveDayGP(level=max(tail(year_max$Value, n=25)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  red=probAboveDayGP(level=max(tail(year_max$Value, n=50)), u=u, year=2024, zeroYear=1960, lambda=pre_lambda, phi0=estimates[1], shape0=estimates[2], cov=cov)
  "
  yellow = 0 
  orange = 0 
  red = 0

  yellow2 = 0 
  orange2 = 0 
  red2 = 0
  "
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
              "lambda0_se" = s[1],
              "lambda1_se" = s[2]*365.25,
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


# Tar namnet på alla filer som finns i data mappen
list <- list.files("data")

# Loopar igenom alla filer
#for (i in list) {
#  PoT_function(i)
#}

print(PoT_function(list[1]))

"
lambda_pot0
low_lambda_pot0
high_lambda_pot0


yellow_tomorrow_pot0 - Sannolikheten att en gul eller värre kommer IMORGON. Poissonprocessen måste alltså ha minst en händelse imorgon och överstigningen måste vara minst gul.
orange_tomorrow_pot0
red_tomorrow_pot0

"
