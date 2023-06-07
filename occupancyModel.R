setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_lynxOccupancy/")
# Appendix Bauduin, ... and Gimenez

library(raster)
library(lubridate)
library(sf)
library(jagsUI)
library(MCMCvis)
library(ggplot2)
library(forcats)

# Occupancy model

###################
# Load covariates #
###################
load("outputs/covariatesCover.RData") # forestCov, connectForCov, shrubCov, openLandCov, agri21Cov, agri22Cov, agri23Cov, agri24Cov
load("outputs/covariatesRoads.RData") # distHgwCov, roadLengthCov
load("outputs/covariatesOthers.RData") # hDensCov, triCov, preyCov, hfiCov
load("data/effort_NovDecJanFebMarApr_1993_2020.RData") # effort

cov <- cbind.data.frame(forestCov, connectForCov[,-1], shrubCov[,-1], 
                        openLandCov[,-1], agri21Cov[,-1], agri22Cov[,-1], 
                        agri23Cov[,-1], agri24Cov[,-1], distHgwCov[,-1], 
                        roadLengthCov[,-1], hDensCov[,-1], triCov[,-1], 
                        preyCov[,-1], hfiCov[,-1])


# Grid cell on which lynx presence data are identified
load("data/gridFrCompleteCrop.RData")

# # Look at the covariates
gridFrCompleteCov <- gridFrComplete
gridFrCompleteCov@data <- cbind.data.frame(gridFrCompleteCov@data, cov, effort)
gridFrCompleteCovSf <- st_as_sf(gridFrCompleteCov)
# for(i in 1:ncol(gridFrCompleteCovSf)){
#   plot(gridFrCompleteCovSf[,i])
# }

###################### 
# Lynx presence data #
######################
lynxData1 <- read.csv("data/lynx/lynx_INDallR.csv", header = TRUE, sep = ";")
lynxData2 <- read.csv("data/lynx/Export_Requete_pour_export.csv", header = TRUE, sep = ";")
lynxData2 <- lynxData2[lynxData2$Fiabilité == "Retenu", ] # confirmed
lynxData3 <- read.csv("data/lynx/Lynx_DOMLNE_1989-2018.csv", header = TRUE, sep = ";")
lynxData4 <- read.csv("data/lynx/Lynx_DOMLNE_2019-2020.csv", header = TRUE, sep = ";")
load("data/lynx/agentsIndicesOCS.RData")
lynxData5 <- agentsIndicesOCS

allLynxData <- cbind.data.frame(date = c(as.Date(lynxData1$Date.observation), 
                                         format(as.Date(lynxData2$Date, format = "%d/%m/%y"), "%Y-%m-%d"), 
                                         as.Date(lynxData3$DateAttaque), 
                                         as.Date(lynxData4$Date.attaque),
                                         as.Date(lynxData5$date)),
                                X = c(lynxData1$X, lynxData2$x_l93, lynxData3$X, lynxData4$X, lynxData5$x),
                                Y = c(lynxData1$Y, lynxData2$y_l93, lynxData3$Y, lynxData4$Y, lynxData5$y))

allLynxData <- allLynxData[!is.na(allLynxData$date),]
allLynxData <- allLynxData[!is.na(allLynxData$X),]
allLynxData <- allLynxData[!is.na(allLynxData$Y),]

# Add the coordinates system
lynxData_spdf <- SpatialPointsDataFrame(coords = allLynxData[,c("X", "Y")], data = allLynxData, 
                                         proj4string = CRS(as.character("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")),
                                         match.ID = TRUE, bbox = NULL)
lynxData_spdf$date <- as.Date(allLynxData$date, format = "%Y-%m-%d")
lynxData_spdfTr <- raster::intersect(spTransform(lynxData_spdf, gridFrComplete@proj4string), gridFrComplete)

# Which cells are occupied during which year
yearStart <- as.numeric(format(min(lynxData_spdfTr$date), "%Y"))
yearEnd <- as.numeric(format(max(lynxData_spdfTr$date), "%Y")) - 1
obsLynx <- array(data = 0, dim = c(length(gridFrComplete$ID), 6, yearEnd - yearStart + 1))

for(year in yearStart:yearEnd){
  
  # Extract the lynx observations for the different time period
  novYear <- lynxData_spdfTr[lynxData_spdfTr$date >= as.Date(paste(year, "11", "01", sep = "-")) &
                               lynxData_spdfTr$date < (as.Date(paste(year, "11", "01", sep = "-")))%m+% months(1),]
  decYear <- lynxData_spdfTr[lynxData_spdfTr$date >= as.Date(paste(year, "12", "01", sep = "-")) &
                               lynxData_spdfTr$date < (as.Date(paste(year, "12", "01", sep = "-")))%m+% months(1),]
  janYearPlus1 <- lynxData_spdfTr[lynxData_spdfTr$date >= as.Date(paste(year + 1, "01", "01", sep = "-")) &
                                    lynxData_spdfTr$date < (as.Date(paste(year + 1, "01", "01", sep = "-")))%m+% months(1),]
  febYearPlus1 <- lynxData_spdfTr[lynxData_spdfTr$date >= as.Date(paste(year + 1, "02", "01", sep = "-")) &
                                    lynxData_spdfTr$date < (as.Date(paste(year + 1, "02", "01", sep = "-")))%m+% months(1),]
  marYearPlus1 <- lynxData_spdfTr[lynxData_spdfTr$date >= as.Date(paste(year + 1, "03", "01", sep = "-")) &
                                    lynxData_spdfTr$date < (as.Date(paste(year + 1, "03", "01", sep = "-")))%m+% months(1),]
  aprYearPlus1 <- lynxData_spdfTr[lynxData_spdfTr$date >= as.Date(paste(year + 1, "04", "01", sep = "-")) &
                                    lynxData_spdfTr$date < (as.Date(paste(year + 1, "04", "01", sep = "-")))%m+% months(1),]
  
  # plot(gridFrComplete, main = paste0("Winter ", year, "-", year + 1))
  # plot(novYear, add = TRUE, col = "red", pch = 16)
  # plot(decYear, add = TRUE, col = "blue", pch = 16)
  # plot(janYearPlus1, add = TRUE, col = "green", pch = 16)
  # plot(febYearPlus1, add = TRUE, col = "yellow", pch = 16)
  # plot(marYearPlus1, add = TRUE, col = "orange", pch = 16)
  # plot(aprYearPlus1, add = TRUE, col = "purple", pch = 16)
  
  # Which cells were occupied
  if(length(novYear) != 0){
    novCells <- over(novYear, gridFrComplete)
    novCellsID <- unique(novCells)[, "ID"]
    novCellsID <- novCellsID[!is.na(novCellsID)]
    obsLynx[which(gridFrComplete$ID %in% novCellsID), 1, year - yearStart + 1] <- 1
  }
  if(length(decYear) != 0){
    decCells <- over(decYear, gridFrComplete)
    decCellsID <- unique(decCells)[, "ID"]
    decCellsID <- decCellsID[!is.na(decCellsID)]
    obsLynx[which(gridFrComplete$ID %in% decCellsID), 2, year - yearStart + 1] <- 1
  }
  if(length(janYearPlus1) != 0){
    janCells <- over(janYearPlus1, gridFrComplete)
    janCellsID <- unique(janCells)[, "ID"]
    janCellsID <- janCellsID[!is.na(janCellsID)]
    obsLynx[which(gridFrComplete$ID %in% janCellsID), 3, year - yearStart + 1] <- 1
  }
  if(length(febYearPlus1) != 0){
    febCells <- over(febYearPlus1, gridFrComplete)
    febCellsID <- unique(febCells)[, "ID"]
    febCellsID <- febCellsID[!is.na(febCellsID)]
    obsLynx[which(gridFrComplete$ID %in% febCellsID), 4, year - yearStart + 1] <- 1
  }
  if(length(marYearPlus1) != 0){
    marCells <- over(marYearPlus1, gridFrComplete)
    marCellsID <- unique(marCells)[, "ID"]
    marCellsID <- marCellsID[!is.na(marCellsID)]
    obsLynx[which(gridFrComplete$ID %in% marCellsID), 5, year - yearStart + 1] <- 1
  }
  if(length(aprYearPlus1) != 0){
    aprCells <- over(aprYearPlus1, gridFrComplete)
    aprCellsID <- unique(aprCells)[, "ID"]
    aprCellsID <- aprCellsID[!is.na(aprCellsID)]
    obsLynx[which(gridFrComplete$ID %in% aprCellsID), 6, year - yearStart + 1] <- 1
  }
}

# We only have the effort from winter 1993-1994 to winter 2020-2021
# Lynx data are only valid until winter 2019-2020
# We will restrict the effort data and lynx data from 1993-1994 to 2019-2020
obsLynx <- obsLynx[,,(1993 - yearStart + 1):(1993 - yearStart + 1 + 2019 - 1993)]
dim(obsLynx) # 1856 (n cells in the grid) x 6 months x 1993 to 2019
effort <- effort[,(1993 - 1993 + 1):(1993 - 1993 + 1 + 2019 - 1993)]

# Identify cells which were never sampled
# Transform the effort in a binary variable
binaryEffort <- ifelse(effort > 0, 1, 0)
# Place NA in effort and in obsLynx for cells not sampled
for(i in 1:nrow(binaryEffort)){
  for(j in 1:ncol(binaryEffort)){
    if(binaryEffort[i, j] == 0){
      obsLynx[i, 1:6, j] <- NA 
    } 
  }
}

# Which cell were never sampled across all years
noSampl <- apply(obsLynx, 1, function(x) all(is.na(x)))
sum(noSampl)
# Remove cells never samples
covSampl <- cov[!noSampl,]
binaryEffortSampl <- binaryEffort[!noSampl,]
obsLynxSampl <- obsLynx[!noSampl,,]

# Dataset characteristics
nSites <- dim(obsLynxSampl)[1]
nVisits <- dim(obsLynxSampl)[2]
nYears <- dim(obsLynxSampl)[3]

# Scaling the covariates
covScale <- scale(covSampl) # scaling by column
# Duplicate values in-between years for the missing years
rep.col <- function(x, n){
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}
forestS <- cbind(rep.col(covScale[, "forest1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "forest2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "forest2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "forest2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "forest2018"], 2)) # from winter 18/19 to winter 19/20
connectForestS <- cbind(rep.col(covScale[, "connectFor90"], 7), # from winter 93/94 to winter 99/00
                          rep.col(covScale[, "connectFor00"], 6), # from winter 00/01 to winter 05/06
                          rep.col(covScale[, "connectFor06"], 6), # from winter 06/07 to winter 11/12
                          rep.col(covScale[, "connectFor12"], 6), # from winter 12/13 to winter 17/18
                          rep.col(covScale[, "connectFor18"], 2)) # from winter 18/19 to winter 19/20
shrubS <- cbind(rep.col(covScale[, "shrub1990"], 7), # from winter 93/94 to winter 99/00
                  rep.col(covScale[, "shrub2000"], 6), # from winter 00/01 to winter 05/06
                  rep.col(covScale[, "shrub2006"], 6), # from winter 06/07 to winter 11/12
                  rep.col(covScale[, "shrub2012"], 6), # from winter 12/13 to winter 17/18
                  rep.col(covScale[, "shrub2018"], 2)) # from winter 18/19 to winter 19/20
openLandS <- cbind(rep.col(covScale[, "openLand1990"], 7), # from winter 93/94 to winter 99/00
                     rep.col(covScale[, "openLand2000"], 6), # from winter 00/01 to winter 05/06
                     rep.col(covScale[, "openLand2006"], 6), # from winter 06/07 to winter 11/12
                     rep.col(covScale[, "openLand2012"], 6), # from winter 12/13 to winter 17/18
                     rep.col(covScale[, "openLand2018"], 2)) # from winter 18/19 to winter 19/20
agri21S <- cbind(rep.col(covScale[, "agri21_1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "agri21_2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "agri21_2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "agri21_2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "agri21_2018"], 2)) # from winter 18/19 to winter 19/20
agri22S <- cbind(rep.col(covScale[, "agri22_1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "agri22_2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "agri22_2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "agri22_2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "agri22_2018"], 2)) # from winter 18/19 to winter 19/20
agri23S <- cbind(rep.col(covScale[, "agri23_1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "agri23_2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "agri23_2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "agri23_2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "agri23_2018"], 2)) # from winter 18/19 to winter 19/20
agri24S <- cbind(rep.col(covScale[, "agri24_1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "agri24_2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "agri24_2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "agri24_2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "agri24_2018"], 2)) # from winter 18/19 to winter 19/20
roadLengthS <- cbind(rep.col(covScale[, "roadLength2012"], 22), # from winter 93/94 to winter 14/15
                       rep.col(covScale[, "roadLength2015"], 3), # from winter 15/16 to winter 17/18
                       rep.col(covScale[, "roadLength2018"], 2)) # from winter 18/19 to winter 19/20
preyS <- cbind(rep.col(covScale[, "prey1993"], 5), # from winter 93/94 to winter 97/98
                 rep.col(covScale[, "prey1998"], 4), # from winter 98/99 to winter 01/02
                 rep.col(covScale[, "prey2002"], 5), # from winter 02/03 to winter 06/07
                 rep.col(covScale[, "prey2007"], 5), # from winter 07/08 to winter 11/12
                 rep.col(covScale[, "prey2012"], 5), # from winter 12/13 to winter 16/17
                 rep.col(covScale[, "prey2017"], 2), # from winter 17/18 to winter 18/19
                 rep.col(covScale[, "prey2019"], 1)) #  winter 19/20
hfiS <- cbind(rep.col(covScale[, "hfi_1993"], 16), # from winter 93/94 to winter 08/09
                rep.col(covScale[, "hfi_2009"], 11)) # from winter 09/10 to winter 19/20

###################
# Occupancy model #
###################
# List of data to run the model
dataFull <- list(nSites = nSites, 
              nVisits = nVisits, 
              nYears = nYears, 
              forest = forestS, 
              connectForest = connectForestS,
              shrub = shrubS,
              openLand = openLandS,
              agri21 = agri21S,
              agri22 = agri22S,
              agri23 = agri23S,
              agri24 = agri24S,
              distHighway = covScale[, "distHgwCov[, -1]"],
              roadLength = roadLengthS,
              hDens = covScale[, "hDensCov[, -1]"],
              tri = covScale[, "triCov[, -1]"],
              prey = preyS,
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
  for(i in 1:13) {
    alpha[i] ~ dnorm(0, 1)
  }
  for(i in 1:nSites) {
    for(j in 1:nYears) {
      logit(phi[i, j]) <- alpha[1] * forest[i, j] + alpha[2] * connectForest[i, j] + alpha[3] * shrub[i, j] + 
        alpha[4] * openLand[i, j] + alpha[5] * agri21[i, j] + alpha[6] * agri22[i, j] + alpha[7] * agri23[i, j] + 
        alpha[8] * agri24[i, j] + alpha[9] * distHighway[i] + alpha[10] * hDens[i] + alpha[11] * tri[i] + 
        alpha[12] * prey[i, j] + alpha[13] * hfi[i, j] + b[j] + u[i]
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

# Initial value
zst <- dataFull$y
zst[zst > 0] <- 1
init1 <- list(alpha = runif(13, -2, 2),
              beta = runif(2, -2, 2),
              z = zst,
              b = rep(0, nYears),
              mu.b = 0,
              u = rep(0, nSites),
              v = matrix(ncol = nYears, nrow = nSites, data = rep(0, nSites * nYears)),
              sd.b = 1,
              sd.u = 1,
              sd.v = 1)
init2 <- list(alpha = runif(13, -2, 2),
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
                "beta",
                "b",
                "mu.b",
                "v",
                "sd.b",
                "sd.u",
                "sd.v",
                "z")

# Run the model
# lynxSim <- jags(data = dataFull,
#                 inits = inits,
#                 parameters = parameters,
#                 n.iter = 20000,
#                 model.file = "modelOcc.txt",
#                 n.chains = 2,
#                 n.burnin = 10000)
# save(lynxSim, file="outputs/lynxZ.RData")
load("outputs/lynxZ.RData")

# Summary results
round(lynxSim$summary, 2)

# Slope of the effects of the covariates on the occupancy probability (alpha) to see if the convergence is satisfied 
MCMCtrace(lynxSim, params = 'alpha', ISB = TRUE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, pdf = FALSE)
# Time effects on the occupancy probability (b) to see if the convergence is satisfied 
MCMCtrace(lynxSim, params = 'b', ISB = TRUE, ind = TRUE, Rhat = TRUE, n.eff = TRUE, pdf = FALSE)

# Effect size:
MCMCplot(object = lynxSim, 
         params = 'alpha', 
         rank = TRUE,
         labels = c("forest", "connectForest", "shrub", "openLand", "agri21", "agri22", "agri23", "agri24",
         "distHighway", "hDens", "tri", "prey", "hfi"))


#####################
# Model predictions #
#####################

# Maps of occupancy probability based on the habitat only

# Function to recreate occupancy probabilities
occEstimate <- function(matrixAlpa, vectorBeta, lengthCov, 
                        covFor, covConnF, covShrub, covOpenL, covA21, covA22, covA23, 
                        covA24, covDistH, covHDens, covTri, covPrey, covHfi){
  
  mapsValues <- matrix(nrow = lengthCov, ncol = nrow(matrixAlpa))
  
  for(i in 1:nrow(matrixAlpa)){
    
    logitVal <- matrixAlpa[i, 1] * covFor + matrixAlpa[i, 2] * covConnF + matrixAlpa[i, 3] * covShrub + 
      matrixAlpa[i, 4] * covOpenL + matrixAlpa[i, 5] * covA21 + matrixAlpa[i, 6] * covA22 + 
      matrixAlpa[i, 7] * covA23 + matrixAlpa[i, 8] * covA24 + matrixAlpa[i, 9] * covDistH + 
      matrixAlpa[i, 10] * covHDens + matrixAlpa[i, 11] * covTri + matrixAlpa[i, 12] * covPrey + 
      matrixAlpa[i, 13] * covHfi + vectorBeta[i]
    
    mapsValues[, i] <- exp(logitVal) / (1 + exp(logitVal))
    #print(i)
  }
  
  return(mapsValues)
}

# lynxSim$sims.list$alpha where each line is an iteration with the 13 alpha values estimated
# lynxSim$sims.list$b for the estimated b
alphaEstim <- lynxSim$sims.list$alpha
bEstim <- lynxSim$sims.list$b[, 27]

# Let's use the more recent values for the covariates
# We need to do the same scaling as for the model
# So we need the mean and sd on covariates where the cells never sampled were removed
# But we compute occupancy probabilities over the whole landscape, even the cells that were never sampled
allFor <- forestCov[, ncol(forestCov)]
allForS <- as.numeric((allFor - mean(covSampl[, "forest2018"])) / sd(covSampl[, "forest2018"]))
allConnF <- connectForCov[, ncol(connectForCov)]
allConnFS <- as.numeric((allConnF - mean(covSampl[, "connectFor18"])) / sd(covSampl[, "connectFor18"]))
allShrub <- shrubCov[, ncol(shrubCov)]
allShrubS <- as.numeric((allShrub - mean(covSampl[, "shrub2018"])) / sd(covSampl[, "shrub2018"]))
allOpenL <- openLandCov[, ncol(openLandCov)]
allOpenLS <- as.numeric((allOpenL - mean(covSampl[, "openLand2018"])) / sd(covSampl[, "openLand2018"]))
allA21 <- agri21Cov[, ncol(agri21Cov)]
allA21S <- as.numeric((allA21 - mean(covSampl[, "agri21_2018"])) / sd(covSampl[, "agri21_2018"]))
allA22 <- agri22Cov[, ncol(agri22Cov)]
allA22S <- as.numeric((allA22 - mean(covSampl[, "agri22_2018"])) / sd(covSampl[, "agri22_2018"]))
allA23 <- agri23Cov[, ncol(agri23Cov)]
allA23S <- as.numeric((allA23 - mean(covSampl[, "agri23_2018"])) / sd(covSampl[, "agri23_2018"]))
allA24 <- agri24Cov[, ncol(agri24Cov)]
allA24S <- as.numeric((allA24 - mean(covSampl[, "agri24_2018"])) / sd(covSampl[, "agri24_2018"]))
allDistH <- distHgwCov[, ncol(distHgwCov)]
allDistHS <- as.numeric((allDistH - mean(covSampl[, "distHgwCov[, -1]"])) / sd(covSampl[, "distHgwCov[, -1]"]))
allDensH <- hDensCov[, ncol(hDensCov)]
allDensHS <- as.numeric((allDensH - mean(covSampl[, "hDensCov[, -1]"])) / sd(covSampl[, "hDensCov[, -1]"]))
allTri <- triCov[, ncol(triCov)]
allTriS <- as.numeric((allTri - mean(covSampl[, "triCov[, -1]"])) / sd(covSampl[, "triCov[, -1]"]))
allPrey <- preyCov[, "prey2020"]
allPreyS <- as.numeric((allPrey - mean(covSampl[, "prey2020"])) / sd(covSampl[, "prey2020"]))
allHfi <- hfiCov[, ncol(hfiCov)]
allHfiS <- as.numeric((allHfi - mean(covSampl[, "hfi_2009"])) / sd(covSampl[, "hfi_2009"]))

gridFrCompleteCovAll <- gridFrCompleteCov
gridFrCompleteCovAll@data <- data.frame(allForS = allForS, allConnFS = allConnFS, allShrubS = allShrubS,
                                        allOpenLS = allOpenLS, allA21S = allA21S, allA22S = allA22S,
                                        allA23S = allA23S, allA24S = allA24S, allDistHS = allDistHS,
                                        allDensHS = allDensHS, allTriS = allTriS, allPreyS = allPreyS,
                                        allHfiS = allHfiS)
# Plot each covariate with and without the point presence
dataP <- st_as_sf(lynxData_spdfTr)
covSF <- st_as_sf(gridFrCompleteCovAll)
# Forest
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allForS)) +
  scale_fill_viridis_c() +
  ggtitle("Forest")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allForS)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Forest")
print(plotWithPoints)
# Forest connectivity
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allConnFS)) +
  scale_fill_viridis_c() +
  ggtitle("Forest connectivity")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allConnFS)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Forest connectivity")
print(plotWithPoints)
# Shrub
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allShrubS)) +
  scale_fill_viridis_c() +
  ggtitle("Shrub")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allShrubS)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Shrub")
print(plotWithPoints)
# Open land
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allOpenLS)) +
  scale_fill_viridis_c() +
  ggtitle("Open land")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allOpenLS)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Open land")
print(plotWithPoints)
# Agri 21
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allA21S)) +
  scale_fill_viridis_c() +
  ggtitle("Agr21")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allA21S)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Agri 21")
print(plotWithPoints)
# Agri 22
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allA22S)) +
  scale_fill_viridis_c() +
  ggtitle("Agri 22")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allA22S)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Agri 22")
print(plotWithPoints)
# Agri 23
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allA23S)) +
  scale_fill_viridis_c() +
  ggtitle("Agri 23")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allA23S)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Agri 23")
print(plotWithPoints)
# Agri 24
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allA24S)) +
  scale_fill_viridis_c() +
  ggtitle("Agri 24")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allA24S)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Agri 24")
print(plotWithPoints)
# Distance to highways
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allDistHS)) +
  scale_fill_viridis_c() +
  ggtitle("Distance to highways")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allDistHS)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Distance to highways")
print(plotWithPoints)
# Human density
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allDensHS)) +
  scale_fill_viridis_c() +
  ggtitle("Human density")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allDensHS)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Human density")
print(plotWithPoints)
# TRI
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allTriS)) +
  scale_fill_viridis_c() +
  ggtitle("Terrain ruggeness")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allTriS)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Terrain ruggeness")
print(plotWithPoints)
# Prey
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allPreyS)) +
  scale_fill_viridis_c() +
  ggtitle("Preys")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allPreyS)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Preys")
print(plotWithPoints)
# HFI
plotWithoutPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allHfiS)) +
  scale_fill_viridis_c() +
  ggtitle("Human footprint")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = covSF, aes(fill = allHfiS)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Human footprint")
print(plotWithPoints)

allMapsValues <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = length(allForS), 
                             covFor = allForS, covConnF = allConnFS, covShrub = allShrubS, covOpenL = allOpenLS, 
                             covA21 = allA21S, covA22 = allA22S, covA23 = allA23S, covA24 = allA24S, covDistH = allDistHS, 
                             covHDens = allDensHS, covTri = allTriS, covPrey = allPreyS, covHfi = allHfiS)

# Compute the mean and quantiles for each cell over all iterations
meanOcc <- rowMeans(allMapsValues, na.rm = TRUE)
quantileOcc <- apply(allMapsValues, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))

mapPotentialOcc <- gridFrCompleteCov
mapPotentialOcc@data <- cbind.data.frame(meanOcc = meanOcc, lowCIocc = quantileOcc[1, ], upCIocc = quantileOcc[2, ])
mapPotentialOcc <- st_as_sf(mapPotentialOcc)
plot(mapPotentialOcc)

# Plot predicted occupancy probability with data points
plotWithoutPoints <- ggplot() + 
  geom_sf(data = mapPotentialOcc, aes(fill = meanOcc)) +
  scale_fill_viridis_c() +
  ggtitle("Mean occupancy predicted")
print(plotWithoutPoints)
plotWithPoints <- ggplot() + 
  geom_sf(data = mapPotentialOcc, aes(fill = meanOcc)) + 
  scale_fill_viridis_c() +
  geom_sf(data = dataP, shape = 1) +
  ggtitle("Mean occupancy predicted")
print(plotWithPoints)

# Create figures of mean occupancy as a function of each covariate while the other covariates are at their mean value

## Tri
# Define the range of value for the covariate
covVal <- as.numeric(allTri)
covValForScaling <- as.numeric(covSampl[, "triCov[, -1]"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      covFor = mean(allForS), covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 = mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covHDens = mean(allDensHS), covTri = covValS, 
                      covPrey = mean(allPreyS), covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(min(quantileOcc), max(quantileOcc)), lwd = 2,
     xlab = "Terrain Ruggeness Index", ylab = "Occupancy probability")
lines(x = covValTest, y = quantileOcc[1,], lty = 2)
lines(x = covValTest, y = quantileOcc[2,], lty = 2)

## Forest
# Define the range of value for the covariate
covVal <- as.numeric(allFor)
covValForScaling <- as.numeric(covSampl[, "forest2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      covFor = covValS, covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 = mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covPrey = mean(allPreyS), covHfi = mean(allHfiS))
# Compute the mean and CI occupancy probabilities
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(min(quantileOcc), max(quantileOcc)), lwd = 2,
     xlab = "Forest cover", ylab = "Occupancy probability")
lines(x = covValTest, y = quantileOcc[1,], lty = 2)
lines(x = covValTest, y = quantileOcc[2,], lty = 2)

## Prey
# Define the range of value for the covariate
covVal <- as.numeric(allPrey)
covValForScaling <- as.numeric(covSampl[, "prey2020"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      covFor = mean(allForS), covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 = mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covPrey = covValS, covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(min(quantileOcc), max(quantileOcc)), lwd = 2,
     xlab = "Preys", ylab = "Occupancy probability")
lines(x = covValTest, y = quantileOcc[1,], lty = 2)
lines(x = covValTest, y = quantileOcc[2,], lty = 2)

## Shrub
# Define the range of value for the covariate
covVal <- as.numeric(allShrub)
covValForScaling <- as.numeric(covSampl[, "shrub2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      covFor = mean(allForS), covConnF = mean(allConnFS), covShrub = covValS, covOpenL = mean(allOpenLS), 
                      covA21 = mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covPrey = mean(allPreyS), covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(min(quantileOcc), max(quantileOcc)), lwd = 2,
     xlab = "Shrub cover", ylab = "Occupancy probability")
lines(x = covValTest, y = quantileOcc[1,], lty = 2)
lines(x = covValTest, y = quantileOcc[2,], lty = 2)

## Agri21
covVal <- as.numeric(allA21)
covValForScaling <- as.numeric(covSampl[, "agri21_2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      covFor = mean(allForS), covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 = covValS, covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covPrey = mean(allPreyS), covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(min(quantileOcc), max(quantileOcc)), lwd = 2,
     xlab = "Agri21 cover", ylab = "Occupancy probability")
lines(x = covValTest, y = quantileOcc[1,], lty = 2)
lines(x = covValTest, y = quantileOcc[2,], lty = 2)

## Agri22
covVal <- as.numeric(allA22)
covValForScaling <- as.numeric(covSampl[, "agri22_2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      covFor = mean(allForS), covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = covValS, covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covPrey = mean(allPreyS), covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(min(quantileOcc), max(quantileOcc)), lwd = 2,
     xlab = "Agri22 cover", ylab = "Occupancy probability")
lines(x = covValTest, y = quantileOcc[1,], lty = 2)
lines(x = covValTest, y = quantileOcc[2,], lty = 2)

## Distance to highways
covVal <- as.numeric(allDistH)
covValForScaling <- as.numeric(covSampl[, "distHgwCov[, -1]"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      covFor = mean(allForS), covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = covValS, covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covPrey = mean(allPreyS), covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(min(quantileOcc), max(quantileOcc)), lwd = 2,
     xlab = "Distance to highways", ylab = "Occupancy probability")
lines(x = covValTest, y = quantileOcc[1,], lty = 2)
lines(x = covValTest, y = quantileOcc[2,], lty = 2)

## Forest connectivity
covVal <- as.numeric(allConnF)
covValForScaling <- as.numeric(covSampl[, "connectFor18"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      covFor = mean(allForS), covConnF = covValS, covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covPrey = mean(allPreyS), covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(min(quantileOcc), max(quantileOcc)), lwd = 2,
     xlab = "Forest connectivity", ylab = "Occupancy probability")
lines(x = covValTest, y = quantileOcc[1,], lty = 2)
lines(x = covValTest, y = quantileOcc[2,], lty = 2)
