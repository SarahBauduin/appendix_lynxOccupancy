# Temporal effect on the distance to highways

dataFull <- list(nSites = nSites, 
                nVisits = nVisits, 
                nYears = nYears, 
                #forest = forestS, 
                connectForest = connectForestS,
                shrub = shrubS,
                openLand = openLandS,
                agri21 = agri21S,
                agri22 = agri22S,
                agri23 = agri23S,
                agri24 = agri24S,
                distHighway = covScale[, "distHgwCov[, -1]"],
                distHighway2 = covScale[, "distHgwCov2[, -1]"],
                roadLength = roadLengthS,
                pathLength = pathLengthS,
                hDens = covScale[, "hDensCov[, -1]"],
                #elevMean = covScale[, "meanElev"],
                #elevSD = covScale[, "sdElev"],
                tri = covScale[, "triCov[, -1]"],
                hfi = hfiS,
                binaryEffort = binaryEffortSampl,
                nSitesYear = apply(binaryEffortSampl, 2, sum),
                y = apply(obsLynxSampl, c(1, 3), sum))

# Model
model <- 
  paste("
model
{
  # State model priors
  b[1] ~ dnorm(mu.b, 1) # Random walk prior on year effect
  for(j in 2:nYears){
    b[j] ~ dnorm(b[j - 1], tau.b)
  }
  tau.b <- 1 / (sd.b * sd.b)
  mu.b ~ dnorm(0, 1)
  sd.b ~ dunif(0, 5) # Half-uniform hyperpriors
  
  # Random site effect 
  for (i in 1:nSites){
    u[i] ~ dnorm(0, tau.u) # On occupancy    
  } 
  tau.u <- 1 / (sd.u * sd.u)
  sd.u ~ dunif(0, 5) # Half-uniform hyperpriors
  
  for (i in 1:nSites){
    for(j in 1:nYears){
      v[i, j] ~ dnorm(0, tau.v) # On detection 
    }
  } 
  tau.v <- 1 / (sd.v * sd.v)
  sd.v ~ dunif(0, 5) # Half-uniform hyperpriors
  
  # Priors for alpha and beta coeffs
  for(i in 1:2) {
    beta[i]  ~ dnorm(0, 1)
  }
  for(i in 1:12) {
    alpha[i] ~ dnorm(0, 1)
  }
  for(j in 1:nYears){
    alphaH[j] ~ dnorm(0, 1)
    alphaH2[j] ~ dnorm(0, 1)
  }
  
  for(i in 1:nSites) {
    for(j in 1:nYears) {
      logit(phi[i, j]) <- alpha[1] * connectForest[i, j] + alpha[2] * shrub[i, j] + 
        alpha[3] * openLand[i, j] + alpha[4] * agri21[i, j] + alpha[5] * agri22[i, j] + alpha[6] * agri23[i, j] + 
        alpha[7] * agri24[i, j] +  alphaH[j] * distHighway[i] +  alphaH2[j] * distHighway2[i] +
        alpha[8] * roadLength[i, j] + alpha[9] * pathLength[i, j] + alpha[10] * hDens[i] + 
        alpha[11] * tri[i] + alpha[12] * hfi[i, j] + b[j] + u[i]
      z[i, j] ~ dbern(phi[i, j]) # True presence/absences states
      lp[i, j] <- beta[1] + beta[2] * roadLength[i, j] + v[i, j]
      p[i, j] <- (1 / (1 + exp(- lp[i, j]))) * (1 - step(-binaryEffort[i, j]))
      y[i, j] ~ dbin(p[i, j] * z[i, j], nVisits) # Likelihood
    }
  }
  
  # Finite sample occupancy - proportion of occupied sites
  for (j in 1:nYears) {
    psi.fs[j] <- sum(z[1:nSites, j]) / nSitesYear[j]
  }
}
")
writeLines(model,"modelOcc.txt")

zst <- dataFull$y
zst[zst > 0] <- 1
init1 <- list(alpha = runif(12, -2, 2),
              alphaH = runif(27, -1, 1),
              alphaH2 = runif(27, -1, 1),
              beta = runif(2, -2, 2),
              z = zst,
              b = rep(0, nYears),
              mu.b = 0,
              u = rep(0, nSites),
              v = matrix(ncol = nYears, nrow = nSites, data = rep(0, nSites * nYears)),
              sd.b = 1,
              sd.u = 1,
              sd.v = 1)
init2 <- list(alpha = runif(12, -2, 2),
              alphaH = runif(27, -1, 1),
              alphaH2 = runif(27, -1, 1),
              beta = runif(2, -2, 2),
              z = zst,
              b = rep(0, nYears),
              mu.b = 0,
              u = rep(0, nSites),
              v = matrix(ncol = nYears, nrow = nSites, data = rep(0, nSites * nYears)),
              sd.b = 1,
              sd.u = 1,
              sd.v = 1)
inits <- list(init1, init2)

# Parameters to be monitored
parameters <- c("psi.fs",
                "alpha",
                "alphaH",
                "alphaH2",
                "beta",
                "b",
                "mu.b",
                "v",
                "sd.b",
                "sd.u",
                "sd.v",
                "z")

# Run the model
lynxSim <- jags(data = dataFull,
                inits = inits,
                parameters = parameters,
                n.iter = 20000,
                model.file = "modelOcc.txt",
                n.chains = 2,
                n.burnin = 10000)

save(lynxSim, file="outputs/lynxMCP100_distHw_fTime.RData")
load("outputs/lynxMCP100_distHw_fTime.RData")

# Function to recreate occupancy probabilities
occEstimate <- function(matrixAlpa, vectorAlphaH, vectorAlphaH2, vectorBeta, lengthCov, 
                        covConnF, covShrub, covOpenL, covA21, covA22, covA23, covA24, 
                        covDistH, covDistH2, covRoadLength, covPathLength, covHDens, 
                        covTri, covHfi){
  
  mapsValues <- matrix(nrow = lengthCov, ncol = nrow(matrixAlpa))
  
  for(i in 1:nrow(matrixAlpa)){
    
    logitVal <- matrixAlpa[i, 1] * covConnF + matrixAlpa[i, 2] * covShrub + 
      matrixAlpa[i, 3] * covOpenL + matrixAlpa[i, 4] * covA21 + matrixAlpa[i, 5] * covA22 + 
      matrixAlpa[i, 6] * covA23 + matrixAlpa[i, 7] * covA24 + vectorAlphaH[i] * covDistH + 
      vectorAlphaH2[i] * covDistH2 + matrixAlpa[i, 8] * covRoadLength +
      matrixAlpa[i, 9] * covPathLength + matrixAlpa[i, 10] * covHDens +
      matrixAlpa[i, 11] * covTri + matrixAlpa[i, 12] * covHfi + vectorBeta[i]
    
    mapsValues[, i] <- exp(logitVal) / (1 + exp(logitVal))
    #print(i)
  }
  
  return(mapsValues)
}

alphaEstim <- lynxSim$sims.list$alpha
alphaHEstim <- lynxSim$sims.list$alphaH
alphaH2Estim <- lynxSim$sims.list$alphaH2
bEstim <- lynxSim$sims.list$b[, 27]


## Distance to highways
covVal <- as.numeric(allDistH)
covValForScaling <- as.numeric(covSampl[, "distHgwCov[, -1]"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
covValForScaling2 <- as.numeric(covSampl[, "distHgwCov2[, -1]"])
covValS2 <- ((covValTest * covValTest) - mean(covValForScaling2)) / sd(covValForScaling2)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value

# Estimate values for each year, changing the coefficient for the distance to highways covariates
valForYearMean <- matrix(nrow = 100, ncol = nYears)
valForYearQmin <- matrix(nrow = 100, ncol = nYears)
valForYearQmax <- matrix(nrow = 100, ncol = nYears)

for(y in 1:nYears){
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, 
                        vectorAlphaH = alphaHEstim[,y], vectorAlphaH2 = alphaH2Estim[,y], 
                        vectorBeta = bEstim, lengthCov = 100,
                        covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                        covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                        covDistH = covValS, covDistH2 = covValS2, covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                        covHDens = mean(allDensHS), covTri = mean(allTriS), 
                        covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  
  valForYearMean[,y] <- meanOcc
  valForYearQmin[,y] <- quantileOcc[1,]
  valForYearQmax[,y] <- quantileOcc[2,]
}

# Plot the mean and quantiles
# One color per year
library(viridis)
colPal <- inferno(n = 27)
plot(x = covValTest, y = valForYearMean[,1], type = "l", ylim = c(min(valForYearMean), max(valForYearMean)), 
     lwd = 2, col = colPal[1],
     xlab = "Distance to highways", ylab = "Occupancy probability")
lines(x = covValTest, y = valForYearQmin[,1], lty = 2, col = colPal[1])
lines(x = covValTest, y = valForYearQmax[,1], lty = 2, col = colPal[1])
for(y in 2:nYears){
  lines(x = covValTest, y = valForYearMean[,y], lwd = 2, col = colPal[y])
  lines(x = covValTest, y = valForYearQmin[,y], lty = 2, col = colPal[y])
  lines(x = covValTest, y = valForYearQmax[,y], lty = 2, col = colPal[y])
  Sys.sleep(2)
}
