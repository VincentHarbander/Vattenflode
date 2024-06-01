# Vincent Harbander
# 2024

library(eva)
library(mvtnorm)

sampleParams <- function(paramMean, paramCov, number=1) {
  return(mvtnorm::rmvnorm(number, paramMean, paramCov))
}

sampleGEV <- function(year, zeroYear, loc0, phi0, shape0, loc1=0, phi1=0, shape1=0) {
  # This function returns a sample from the maximum flow of a specific year, simulated
  
  # year is the current year you want to sample from
  # check levelResultsGEV for the meaning of the rest of the arguments
  
  mu <- loc0+loc1*(year-zeroYear)
  sigma <- exp(phi0+phi1*(year-zeroYear))
  xi <- shape0+shape1*(year-zeroYear)
  
  # eva::rgevr can take lists of parameters, BUT returns dependent samples! A bug in the library! Always gives the same sample for the same parameters!
  return(eva::rgevr(1, 1, mu, sigma, xi))
}

sampleVectorGEV <- function(year, zeroYear, loc0, phi0, shape0, loc1=c(), phi1=c(), shape1=c()) {
  # This function returns a VECTOR OF INDEPENDENT of samples from the maximum flow of a specific year, simulated
  
  # NOTE: The parameter arguments here should be vectors, unlike in sampleGEV!
  
  # year is the current year you want to sample from
  # check levelResultsGEV for the meaning of the rest of the arguments
  
  # Assumes all of phi0 and shape0 have the same length as loc0
  n <- length(loc0)
  res <- rep(0, n)
  
  # If no trend was specified, have it be 0
  if(length(loc1) == 0) {
    loc1 <- rep(0, n)
  }
  if(length(phi1) == 0) {
    phi1 <- rep(0, n)
  }
  if(length(shape1) == 0) {
    shape1 <- rep(0, n)
  }
  
  res <- rep(0, n)
  
  for(i in 1:n) {
    res[i] <- sampleGEV(year, zeroYear, loc0[i], phi0[i], shape0[i], loc1[i], phi1[i], shape1[i])
  }
  
  return(res)
}


dataGEV <- function(T, year, zeroYear, loc0, phi0, shape0, loc1=0, phi1=0, shape1=0, cov=matrix(0, 6, 6), samples=10000) {
  # This function returns results of sample futures over T years
  
  # see the rest of the arguments from levelDataGEV
  
  # Sample parameters for each sample
  # paramVectors is a samples x dim(paramMeans) matrix
  paramVectors <- sampleParams(c(loc0, phi0, shape0, loc1, phi1, shape1), cov, samples)
  l0 <- paramVectors[,1]
  p0 <- paramVectors[,2]
  s0 <- paramVectors[,3]
  l1 <- paramVectors[,4]
  p1 <- paramVectors[,5]
  s1 <- paramVectors[,6]
  
  res <- matrix(0, samples, T)
  
  # The first one is necessarily max, after that we have to pick the max so far
  # t=1 is handled here
  res[,1] <- sampleVectorGEV(year+1, zeroYear, l0, p0, s0, l1, p1, s1)
  
  # Has +t to NOT include the current year. NOTE that is starts on 2 since 1 is handled above
  if(T > 1) {
    for(t in 2:T) {
      # Pass vectors of parameters
      res[,t] <- sampleVectorGEV(year+t, zeroYear, l0, p0, s0, l1, p1, s1)
    }
  }
  
  return(res)
}

# Need to not just take a sample, but take t *timed* samples and take the max of those
# This takes a new future sample and returns the biggest one seen so far
levelDataGEV <- function(T, year, zeroYear, loc0, phi0, shape0, loc1=0, phi1=0, shape1=0, cov=matrix(0, 6, 6), samples=10000) {
  # This returns the samples from max{X_1,...,X_T} over times t=1:T in a matrix
  
  # phi=ln(sigma)
  # T is the time period, the number of years.
  # see the rest of the arguments from levelResultsGEV
  
  # Sample parameters for each sample
  # paramVectors is a samples x dim(paramMeans) matrix
  paramVectors <- sampleParams(c(loc0, phi0, shape0, loc1, phi1, shape1), cov, samples)
  l0 <- paramVectors[,1]
  p0 <- paramVectors[,2]
  s0 <- paramVectors[,3]
  l1 <- paramVectors[,4]
  p1 <- paramVectors[,5]
  s1 <- paramVectors[,6]
  
  res <- matrix(0, samples, T)
  
  # The first one is necessarily max, after that we have to pick the max so far
  # t=1 is handled here
  res[,1] <- sampleVectorGEV(year+1, zeroYear, l0, p0, s0, l1, p1, s1)
  
  # Has +t to NOT include the current year. NOTE that is starts on 2 since 1 is handled above
  if(T > 1) {
    for(t in 2:T) {
      # Pass vectors of parameters
      sample <- sampleVectorGEV(year+t, zeroYear, l0, p0, s0, l1, p1, s1)
      # If we found something even bigger, update the new max
      # Perform an element-wise max
      res[,t] <- pmax(res[,t-1], sample)
    }
  }
  
  return(res)
}

levelResultsGEV <- function(timePeriod, startYear, zeroYear, loc0, phi0, shape0, loc1=0, phi1=0, shape1=0, cov=matrix(0, 6, 6), samples=10000, confidence=c(0.95)) {
  # This function returns a vector of average simulated return levels for every year in a given time period and a [confidence]-prediction interval for the level, all simulated
  
  # timePeriod (T) is the time period you want to check the return level for
  # startYear is the year you want to check from (not including startYear)
  # zeroYear is the year where location=loc0, at zeroYear+1 location=loc0+loc1, etc
  # loc refers to mu, the location
  # phi = ln(sigma), where sigma is the scale
  # shape refers to xi
  # confidence gives the function the confidence for the inside of the prediction interval
  # confidence is given as a list of decreasing confidences
  # samples gives how many samples are sampled PER YEAR to determine mean and prediction interval, more is better
  # The cov is the covariance matrix of the parameters listed as mu0, phi0, xi0, mu1, phi1, xi1
  # The parameters are sampled normally around the parameter with the given standard deviation
  
  resData <- levelDataGEV(timePeriod, startYear, zeroYear, loc0, phi0, shape0, loc1, phi1, shape1, cov, samples)
  
  means <- colMeans(resData)
  
  alpha <- 1-confidence
  
  pred <- matrix(0, 2*length(confidence), timePeriod)
  
  for(t in 1:timePeriod) {
    pred[,t] <- quantile(resData[,t], c(alpha/2, rev(1-alpha/2)))
  }
  
  # deparse.level just makes sure it doesn't create any labels for the dataframe
  return(rbind(means, pred, deparse.level=0))
  #return(c(means, pred))
}

levelPlotGEV <- function(timePeriod, startYear, zeroYear, loc0, phi0, shape0, loc1=0, phi1=0, shape1=0, cov=matrix(0, 6, 6), samples=10000, confidence=(0.95)) {
  # This function simply calls levelResultsGEV and handles the plotting, returning nothing
  # See levelResultsGEV for the meaning of the arguments
  
  years <- matrix(rep((startYear+1):(startYear+timePeriod), 3), nrow=3, byrow=TRUE)
  result <- levelResultsGEV(timePeriod, startYear, zeroYear, loc0, phi0, shape0, loc1, phi1, shape1, cov, samples, confidence)
  
  n <- length(confidence)
  ramp <- colorRampPalette(c("orange", "red"), alpha=TRUE)(n)
  colors <- c(ramp, rev(ramp))
  
  conf <- as.character(100*confidence)

  for(i in 1:length(conf)) {
    conf[i] <- paste(conf[i], "%", sep="")
  }
  
  # Plots COLUMNS
  matplot(t(years), t(result), type="l", lty=1, 
          col=c("blue", colors), xlab="Year",
          ylab="Max encountered level", main="Prediction level plot")
  legend("bottomright", legend=c("Mean", conf[1:n]),
         col=c("blue", ramp),
         lty=1)
}

probAboveGEV <- function(level, timePeriod, startYear, zeroYear, loc0, phi0, shape0, loc1=0, phi1=0, shape1=0, cov=matrix(0, 6, 6), samples=1000000) {
  # Returns the probability of exceeding level within the specified time period
  # startYear is NOT included in the time period analyzed
  # See levelResultsGEV for the rest of the arguments
  
  # resData contains samples of the biggest occurrance within the time period
  tmp <- levelDataGEV(timePeriod, startYear, zeroYear, loc0, phi0, shape0, loc1, phi1, shape1, cov, samples)
  resData <- tmp[,timePeriod] # Pick out the last time, it's necessarily the max
  
  # Now we count exceedances!
  exceedances <- 0
  
  for(elem in resData) {
    if(elem > level) {
      exceedances <- exceedances + 1
    }
  }

  return(exceedances/samples)
}

sampleGP <- function(year, zeroYear, phi0, shape0, phi1=0, shape1=0) {
  # This function returns a sample from the maximum flow of a specific year, simulated
  
  # year is the current year you want to sample from
  # check levelResultsGEV for the meaning of the rest of the arguments
  
  sigma <- exp(phi0+phi1*(year-zeroYear))
  xi <- shape0+shape1*(year-zeroYear)
  
  # eva::rgpd is also bugged with parameters as vectors!!! Need to write my own vectorization!
  return(eva::rgpd(1, 0, sigma, xi))
}

sampleVectorGP <- function(year, zeroYear, phi0, shape0, phi1=c(), shape1=c()) {
  # This function returns a VECTOR of INDEPENDENT samples of GP given a vector of parameters and years, simulated
  
  # See sampleGP for the arguments
  
  # Assumes shape0 has the same length as phi0
  n <- length(phi0)
  res <- rep(0, n)
  
  # If no trend was specified, have it be 0
  if(length(phi1) == 0) {
    phi1 <- rep(0, n)
  }
  if(length(shape1) == 0) {
    shape1 <- rep(0, n)
  }
  
  res <- rep(0, n)
  
  for(i in 1:n) {
    res[i] <- sampleGP(year, zeroYear, phi0[i], shape0[i], phi1[i], shape1[i])
  }
  
  return(res)
}

probAboveGP <- function(level, u, year, zeroYear, lambda, phi0, shape0, phi1=0, shape1=0, cov=matrix(0, 4, 4), samples=10000000) {
  # This function returns the probability that an occurrance above level happens during year using PoT (Peaks over threshold)
  
  # year is the EXACT year you want to check
  # ASSUMES level > u!
  # u is the threshold for PoT
  # cov is the covariance matrix of the parameters listed as phi0, xi0, phi1, xi1
  
  excessLevel <- level - u
  samplePoi <- rpois(samples, lambda) # Sample the number of threshold excesses for each sample future
  
  # Sample parameters as vectors!
  sampleVectors <- sampleParams(c(phi0, shape0, phi1, shape1), cov, samples)
  
  p0 <- sampleVectors[,1]
  s0 <- sampleVectors[,2]
  p1 <- sampleVectors[,3]
  s1 <- sampleVectors[,4]
  
  exceedances <- 0
  
  for(i in 1:samples) {
    n <- samplePoi[i]
    if(n == 0) { next }
    for(j in 1:n) { # For each occurance of this sample year
      # Note! Not sampleVectorGP!
      res <- sampleGP(year, zeroYear, p0[i], s0[i], p1[i], s1[i])
      
      # Is this high enough to be relevant?
      if(res > excessLevel) {
        exceedances <- exceedances + 1
        # Need to break here! We only care if an exceedance occurs or not.
        break
      }
    }
  }
  
  return(exceedances/samples)
}

probAboveDayGP <- function(level, u, year, zeroYear, lambda, phi0, shape0, phi1=0, shape1=0, cov=matrix(0, 4, 4), samples=365*1000000) {
  # This function returns the probability that an occurrance above level happens during A DAY using PoT (Peaks over threshold)
  
  # year is the EXACT year you want to check
  # ASSUMES level > u!
  # u is the threshold for PoT
  # cov is the covariance matrix of the parameters listed as phi0, xi0, phi1, xi1
  
  # Use Poisson thinning and increased default samples to compensate
  return(probAboveGP(level, u, year, zeroYear, lambda/365, phi0, shape0, phi1, shape1, cov, samples))
}

# Ankarvattnet_1537 enligt BM
#mu0 <- 106.005731007067
#mu1 <- 0.286752012672177
#sigma <- 33.8274083582117
#xi <- -0.335594561819737

# Made up, purely for testing
#cov <- matrix(0, 6, 6)
#cov[1, 1] <- 1
#cov[4, 4] <- 0.01
#cov[1, 4] <- 0.01
#cov[4, 1] <- cov[1, 4]

# Ankarvattnet_1537 enligt PoT
#lambda <- 0.01
#gpphi <- 3.4
#gpxi <- -0.23

# Made up, purely for testing
#gpcov <- matrix(0, 4, 4)
#gpcov[3,3] <- 0.3

#levelPlotGEV(timePeriod=30, startYear=2024, zeroYear=1960, loc0=mu0, phi0=log(sigma), shape0=xi, loc1=mu1, cov=cov, confidence=(19:1)/20)
