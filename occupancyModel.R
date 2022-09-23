setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_lynxOccupancy/")
# Appendix Bauduin, ... and Gimenez

library(raster)
library(lubridate)
library(sf)
library(jagsUI)

# Occupancy model

# Load covariates
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
load("data/gridFrComplete.RData")
# Crop the grid to the study area extent (east of France)
gridFrComplete <- crop(gridFrComplete, extent(c(720000, 1090000, 6260000, 6920000)))

# # Look at the covariates
# gridFrCompleteCov <- gridFrComplete
# gridFrCompleteCov@data <- cbind.data.frame(gridFrCompleteCov@data, cov, effort)
# gridFrCompleteCovSf <- st_as_sf(gridFrCompleteCov)
# for(i in 1:ncol(gridFrCompleteCovSf)){
#   plot(gridFrCompleteCovSf[,i])
# }

# Lynx presence data
lynxData1 <- read.csv("data/lynx/lynx_INDallR.csv", header = TRUE, sep = ";")
lynxData2 <- read.csv("data/lynx/Export_Requete_pour_export.csv", header = TRUE, sep = ";")
lynxData2 <- lynxData2[lynxData2$Fiabilité == "Retenu", ] # confirmed
lynxData3 <- read.csv("data/lynx/Lynx_DOMLNE_1989-2018.csv", header = TRUE, sep = ";")
lynxData4 <- read.csv("data/lynx/Lynx_DOMLNE_2019-2020.csv", header = TRUE, sep = ";")

allLynxData <- cbind.data.frame(date = c(as.Date(lynxData1$Date.observation), 
                                         format(as.Date(lynxData2$Date, format = "%d/%m/%y"), "%Y-%m-%d"), 
                                         as.Date(lynxData3$DateAttaque), 
                                         as.Date(lynxData4$Date.attaque)),
                                X = c(lynxData1$X, lynxData2$x_l93, lynxData3$X, lynxData4$X),
                                Y = c(lynxData1$Y, lynxData2$y_l93, lynxData3$Y, lynxData4$Y))

allLynxData <- allLynxData[!is.na(allLynxData$date),]
allLynxData <- allLynxData[!is.na(allLynxData$X),]
allLynxData <- allLynxData[!is.na(allLynxData$Y),]

# Add the coordinates system
lynxData_spdf <- SpatialPointsDataFrame(coords = allLynxData[,c("X", "Y")], data = allLynxData, 
                                         proj4string = CRS(as.character("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")),
                                         match.ID = TRUE, bbox = NULL)
lynxData_spdf$date <- as.Date(allLynxData$date, format = "%Y-%m-%d")
lynxData_spdfTr <- spTransform(lynxData_spdf, gridFrComplete@proj4string)

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
  
  plot(gridFrComplete, main = paste0("Winter ", year, "-", year + 1))
  plot(novYear, add = TRUE, col = "red", pch = 16)
  plot(decYear, add = TRUE, col = "blue", pch = 16)
  plot(janYearPlus1, add = TRUE, col = "green", pch = 16)
  plot(febYearPlus1, add = TRUE, col = "yellow", pch = 16)
  plot(marYearPlus1, add = TRUE, col = "orange", pch = 16)
  plot(aprYearPlus1, add = TRUE, col = "purple", pch = 16)
  
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
dim(obsLynx) # 2131 (n cells in the grid) x 6 months x 1993 to 2019
effort <- effort[,(1993 - 1993 + 1):(1993 - 1993 + 1 + 2019 - 1993)]

# Identify cells which were never sampled
# Transform the effort in a binary variable
samplEffort <- ifelse(effort > 0, 1, 0)
# Place NA in effort and in obsLynx for cells not sampled
samplEffort[samplEffort == 0] <- NA
for(i in 1:nrow(samplEffort)){
  for(j in 1:ncol(samplEffort)){
    if(is.na(samplEffort[i, j])){
      obsLynx[i, 1:6, j] <- NA 
    } 
  }
}

# Which cell were never sampled across all years
noSampl <- apply(obsLynx, 1, function(x) all(is.na(x)))
sum(noSampl)
# Remove cells never samples
covSampl <- cov[!noSampl,]
samplEffortSampl <- samplEffort[!noSampl,]
obsLynxSampl <- obsLynx[!noSampl,,]

# Dataset characteristics
nSites <- dim(obsLynxSampl)[1]
nVisits <- dim(obsLynxSampl)[2]
nYears <- dim(obsLynxSampl)[3]

# Scaling the covariates
covScale <- scale(covSampl) # scaling by column
binaryEffort <- ifelse(samplEffortSampl > 0, 1, 0)
effortScale <- scale(samplEffortSampl)
# Duplicate values in-between years for the missing years
rep.col <- function(x, n){
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}
forestCov <- cbind(rep.col(covScale[, "forest1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "forest2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "forest2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "forest2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "forest2018"], 2)) # from winter 18/19 to winter 19/20
connectForestCov <- cbind(rep.col(covScale[, "connectFor90"], 7), # from winter 93/94 to winter 99/00
                          rep.col(covScale[, "connectFor00"], 6), # from winter 00/01 to winter 05/06
                          rep.col(covScale[, "connectFor06"], 6), # from winter 06/07 to winter 11/12
                          rep.col(covScale[, "connectFor12"], 6), # from winter 12/13 to winter 17/18
                          rep.col(covScale[, "connectFor18"], 2)) # from winter 18/19 to winter 19/20
shrubCov <- cbind(rep.col(covScale[, "shrub1990"], 7), # from winter 93/94 to winter 99/00
                  rep.col(covScale[, "shrub2000"], 6), # from winter 00/01 to winter 05/06
                  rep.col(covScale[, "shrub2006"], 6), # from winter 06/07 to winter 11/12
                  rep.col(covScale[, "shrub2012"], 6), # from winter 12/13 to winter 17/18
                  rep.col(covScale[, "shrub2018"], 2)) # from winter 18/19 to winter 19/20
openLandCov <- cbind(rep.col(covScale[, "openLand1990"], 7), # from winter 93/94 to winter 99/00
                     rep.col(covScale[, "openLand2000"], 6), # from winter 00/01 to winter 05/06
                     rep.col(covScale[, "openLand2006"], 6), # from winter 06/07 to winter 11/12
                     rep.col(covScale[, "openLand2012"], 6), # from winter 12/13 to winter 17/18
                     rep.col(covScale[, "openLand2018"], 2)) # from winter 18/19 to winter 19/20
agri21Cov <- cbind(rep.col(covScale[, "agri21_1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "agri21_2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "agri21_2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "agri21_2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "agri21_2018"], 2)) # from winter 18/19 to winter 19/20
agri22Cov <- cbind(rep.col(covScale[, "agri22_1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "agri22_2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "agri22_2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "agri22_2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "agri22_2018"], 2)) # from winter 18/19 to winter 19/20
agri23Cov <- cbind(rep.col(covScale[, "agri23_1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "agri23_2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "agri23_2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "agri23_2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "agri23_2018"], 2)) # from winter 18/19 to winter 19/20
agri24Cov <- cbind(rep.col(covScale[, "agri24_1990"], 7), # from winter 93/94 to winter 99/00
                   rep.col(covScale[, "agri24_2000"], 6), # from winter 00/01 to winter 05/06
                   rep.col(covScale[, "agri24_2006"], 6), # from winter 06/07 to winter 11/12
                   rep.col(covScale[, "agri24_2012"], 6), # from winter 12/13 to winter 17/18
                   rep.col(covScale[, "agri24_2018"], 2)) # from winter 18/19 to winter 19/20
distHighwayCov <- cbind(rep.col(covScale[, "distHgw2012"], 22), # from winter 93/94 to winter 14/15
                        rep.col(covScale[, "distHgw2015"], 3), # from winter 15/16 to winter 17/18
                        rep.col(covScale[, "distHgw2018"], 2)) # from winter 18/19 to winter 19/20
roadLengthCov <- cbind(rep.col(covScale[, "roadLength2012"], 22), # from winter 93/94 to winter 14/15
                       rep.col(covScale[, "roadLength2015"], 3), # from winter 15/16 to winter 17/18
                       rep.col(covScale[, "roadLength2018"], 2)) # from winter 18/19 to winter 19/20
preyCov <- cbind(rep.col(covScale[, "prey1993"], 5), # from winter 93/94 to winter 97/98
                 rep.col(covScale[, "prey1998"], 4), # from winter 98/99 to winter 01/02
                 rep.col(covScale[, "prey2002"], 5), # from winter 02/03 to winter 06/07
                 rep.col(covScale[, "prey2007"], 5), # from winter 07/08 to winter 11/12
                 rep.col(covScale[, "prey2012"], 5), # from winter 12/13 to winter 16/17
                 rep.col(covScale[, "prey2017"], 2), # from winter 17/18 to winter 18/19
                 rep.col(covScale[, "prey2019"], 1)) #  winter 19/20
hfiCov <- cbind(rep.col(covScale[, "hfi_1993"], 16), # from winter 93/94 to winter 08/09
                rep.col(covScale[, "hfi_2009"], 11)) # from winter 09/10 to winter 19/20

# List of data to run the model
dataFull <- list(nSites = nSites, 
              nVisits = nVisits, 
              nYears = nYears, 
              effort = effortScale, 
              forest = forestCov, 
              connectForest = connectForestCov,
              shrub = shrubCov,
              openLand = openLandCov,
              agri21 = agri21Cov,
              agri22 = agri22Cov,
              agri23 = agri23Cov,
              agri24 = agri24Cov,
              distHighway = distHighwayCov,
              roadLengthCov = roadLengthCov,
              hDens = covScale[, "hDensCov[, -1]"],
              tri = covScale[, "triCov[, -1]"],
              prey = preyCov,
              hfi = hfiCov,
              binaryEffort = binaryEffort,
              nSitesYear = apply(binaryEffort, 2, sum),
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
  
  for (i in 1:nSites){
    u[i] ~ dnorm(0, tau.u) # Random site effect      
  } 
  tau.u <- 1 / (sd.u * sd.u)
  sd.u ~ dunif(0, 5) # Half-uniform hyperpriors
  
  # Priors for alpha and beta coeffs
  for(i in 1:2) {
    beta[i]  ~ dnorm(0, 1)
  }
  for(i in 1:14) {
    alpha[i] ~ dnorm(0, 1)
  }
  for(i in 1:nSites) {
    for(j in 1:nYears) {
      logit(phi[i, j]) <- alpha[1] * forest[i, j] + alpha[2] * connectForest[i, j] + alpha[3] * shrub[i, j] + 
        alpha[4] * openLand[i, j] + alpha[5] * agri21[i, j] + alpha[6] * agri22[i, j] + alpha[7] * agri23[i, j] + 
        alpha[8] * agri24[i, j] + alpha[9] * distHighway[i, j] + alpha[10] * roadLengthCov[i, j] + 
        alpha[11] * hDens[i] + alpha[12] * tri[i] + alpha[13] * prey[i, j] + alpha[14] * hfi[i, j] +
        b[j] + u[i]
      z[i, j] ~ dbern(phi[i, j]) # True presence/absences states
      lp[i, j] <- beta[1] + beta[2] * effort[i, j]
      p[i, j] <- (1 / (1 + exp(- lp[i, j]))) * (1 - step(- binaryEffort[i, j]))
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
init1 <- list(alpha = runif(10, -2, 2),
              beta = runif(2, -2, 2),
              z = zst,
              b = rep(0, nYears),
              mu.b = 0,
              u = rep(0, nSites),
              sd.b = 1,
              sd.u = 1)
init2 <- list(alpha = runif(10, -2, 2),
              beta = runif(2, -2, 2),
              z = zst,
              b = rep(0, nYears),
              mu.b = 0,
              u = rep(0, nSites),
              sd.b = 1,
              sd.u = 1)
inits <- list(init1, init2)

# Parameters to be monitored
parameters <- c("psi.fs",
                "alpha",
                "beta",
                "b",
                "mu.b",
                "sd.u",
                "sd.b",
                "z")

# Run the model
lynxSim <- jags(data = dataFull, 
                inits = inits, 
                parameters = parameters,
                n.iter = 10000,
                model.file = "ModelOcc.txt",
                n.chains = 2,
                n.burnin = 2500)
save(lynxSim, file="outputs/lynxZ.RData")

# Summary results
round(lynxSim$summary, 2)



