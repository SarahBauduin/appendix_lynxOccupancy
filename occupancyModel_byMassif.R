# Adding an effect of the massif in the model
library(rgdal)

# Load the limits of the different massifs
alps <- readOGR("C:/Users/sarah.bauduin/Documents/SaveOFB/Projets/13.OccupancyLynx/limites massifs/proposition_limites_ecologiques_lynx/alpes_zonage_presence reguliere.shp")
alps@data$Zone <- "Alps"
bfc <- readOGR("C:/Users/sarah.bauduin/Documents/SaveOFB/Projets/13.OccupancyLynx/limites massifs/proposition_limites_ecologiques_lynx/BFC_zonage_presence reguliere.shp")
bfc@data$Zone <- "BFC"
jura <- readOGR("C:/Users/sarah.bauduin/Documents/SaveOFB/Projets/13.OccupancyLynx/limites massifs/proposition_limites_ecologiques_lynx/jura_zonage_presence reguliere.shp")
jura@data$Zone <- "Jura"
massif <- readOGR("C:/Users/sarah.bauduin/Documents/SaveOFB/Projets/13.OccupancyLynx/limites massifs/proposition_limites_ecologiques_lynx/massif central_zonage_presence reguliere.shp")
massif@data$Zone <- "MassifC"
vosges <- readOGR("C:/Users/sarah.bauduin/Documents/SaveOFB/Projets/13.OccupancyLynx/limites massifs/proposition_limites_ecologiques_lynx/vosges_zonage_presence reguliere.shp")
vosges@data$Zone <- "Vosges"

# Assign 1 or 0 in the grid cell if part of the massif or not
# Combine all massifs together
limits <- rbind(alps, bfc, jura, massif, vosges)
library(rgeos)
ptsGrids <- gCentroid(gridFrCompleteCov, byid = TRUE)
ptsGrids$ID <- 1:length(ptsGrids)
zonesGrid <- intersect(limits, ptsGrids)
missingPts <- setdiff(1:length(ptsGrids), zonesGrid$ID)
zonesData <- rbind(zonesGrid@data[,c("Zone", "ID")], cbind.data.frame(Zone = NA, ID = missingPts))
zonesDataOrder <- zonesData[order(zonesData$ID),]

gridFrCompleteCov$ZoneAlps <- ifelse(zonesDataOrder$Zone == "Alps", 1, 0)
gridFullAlps <- gridFrCompleteCov$ZoneAlps
gridFrCompleteCov$ZoneJura <- ifelse(zonesDataOrder$Zone == "Jura", 1, 0)
gridFullJura <- gridFrCompleteCov$ZoneJura
gridFrCompleteCov$ZoneVosges <- ifelse(zonesDataOrder$Zone == "Vosges", 1, 0)
gridFullVosges <- gridFrCompleteCov$ZoneVosges
gridFrCompleteCov$ZoneBFC <- ifelse(zonesDataOrder$Zone == "BFC", 1, 0)
gridFullBFC <- gridFrCompleteCov$ZoneBFC
gridFrCompleteCovSf <- st_as_sf(gridFrCompleteCov)

## MCP
lynxMCP <- mcp(lynxData_spdfTr, percent = 100)
#lynxMCP <- raster::buffer(SpatialPoints(lynxData_spdfTr@coords,proj4string = lynxData_spdfTr@proj4string), width = 50000)
#lynxMCP <- LoCoH.k(lynxData_spdfTr)
cellSelected <- st_intersection(gridFrCompleteCovSf, st_as_sf(lynxMCP))
gridFrCompleteCovSf <- gridFrCompleteCovSf %>% filter(ID %in% cellSelected$ID)
gridFrCompleteCov <- gridFrCompleteCov[gridFrCompleteCov$ID %in% cellSelected$ID,]
effort <- st_drop_geometry(gridFrCompleteCovSf[,58:85])
cov <- st_drop_geometry(gridFrCompleteCovSf[,c(2,4:89)])
gridFrComplete <- gridFrComplete[gridFrComplete$ID %in% cellSelected$ID,]
forestCov <- st_drop_geometry(gridFrCompleteCovSf[,4:8])
connectForCov <- st_drop_geometry(gridFrCompleteCovSf[,9:13])
shrubCov <- st_drop_geometry(gridFrCompleteCovSf[,14:18])
openLandCov <- st_drop_geometry(gridFrCompleteCovSf[,19:23])
agri21Cov <- st_drop_geometry(gridFrCompleteCovSf[,24:28])
agri22Cov  <- st_drop_geometry(gridFrCompleteCovSf[,29:33])
agri23Cov <- st_drop_geometry(gridFrCompleteCovSf[,34:38])
agri24Cov <- st_drop_geometry(gridFrCompleteCovSf[,39:43])
distHgwCov <- st_drop_geometry(gridFrCompleteCovSf[,44])
distHgwCov2 <- st_drop_geometry(gridFrCompleteCovSf[,45])
roadLengthCov <- st_drop_geometry(gridFrCompleteCovSf[,46:48])
pathLengthCov <- st_drop_geometry(gridFrCompleteCovSf[,49:51])
hDensCov <- st_drop_geometry(gridFrCompleteCovSf[,52])
elevCov <- st_drop_geometry(gridFrCompleteCovSf[,53:54])
triCov <- st_drop_geometry(gridFrCompleteCovSf[,55])
hfiCov <- st_drop_geometry(gridFrCompleteCovSf[,56:57])

####

# List of data to run the model
dataFull <- list(nSites = nSites, 
                 nVisits = nVisits, 
                 nYears = nYears, 
                 ZoneAlps = gridFrCompleteCovSf$ZoneAlps,
                 ZoneJura = gridFrCompleteCovSf$ZoneJura,
                 ZoneVosges = gridFrCompleteCovSf$ZoneVosges,
                 ZoneBFC = gridFrCompleteCovSf$ZoneBFC,
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
  for(i in 1:74) {
    alpha[i] ~ dnorm(0, 1)
  }
  for(i in 1:nSites) {
    for(j in 1:nYears) {
      logit(phi[i, j]) <- 
        
        alpha[1] * ZoneAlps[i] + alpha[2] * ZoneJura[i] + alpha[3] * ZoneVosges[i] + alpha[4] * ZoneBFC[i] +
        
        alpha[5] * connectForest[i, j] + alpha[6] * connectForest[i, j] * ZoneAlps[i] +
                                         alpha[7] * connectForest[i, j] * ZoneJura[i] +
                                         alpha[8] * connectForest[i, j] * ZoneVosges[i] +
                                         alpha[9] * connectForest[i, j] * ZoneBFC[i] +

        alpha[10] * shrub[i, j] + alpha[11] * shrub[i, j] * ZoneAlps[i] +
                                  alpha[12] * shrub[i, j] * ZoneJura[i] +
                                  alpha[13] * shrub[i, j] * ZoneVosges[i] +
                                  alpha[14] * shrub[i, j] * ZoneBFC[i] +

        alpha[15] * openLand[i, j] + alpha[16] * openLand[i, j] * ZoneAlps[i] +
                                     alpha[17] * openLand[i, j] * ZoneJura[i] +
                                     alpha[18] * openLand[i, j] * ZoneVosges[i] +
                                     alpha[19] * openLand[i, j] * ZoneBFC[i] +

        alpha[20] * agri21[i, j] + alpha[21] * agri21[i, j] * ZoneAlps[i] +
                                   alpha[22] * agri21[i, j] * ZoneJura[i] +
                                   alpha[23] * agri21[i, j] * ZoneVosges[i] +
                                   alpha[24] * agri21[i, j] * ZoneBFC[i] +

        alpha[25] * agri22[i, j] + alpha[26] * agri22[i, j] * ZoneAlps[i] +
                                   alpha[27] * agri22[i, j] * ZoneJura[i] +
                                   alpha[28] * agri22[i, j] * ZoneVosges[i] +
                                   alpha[29] * agri22[i, j] * ZoneBFC[i] +

        alpha[30] * agri23[i, j] + alpha[31] * agri23[i, j] * ZoneAlps[i] +
                                   alpha[32] * agri23[i, j] * ZoneJura[i] +
                                   alpha[33] * agri23[i, j] * ZoneVosges[i] +
                                   alpha[34] * agri23[i, j] * ZoneBFC[i] +

        alpha[35] * agri24[i, j] + alpha[36] * agri24[i, j] * ZoneAlps[i] +
                                   alpha[37] * agri24[i, j] * ZoneJura[i] +
                                   alpha[38] * agri24[i, j] * ZoneVosges[i] +
                                   alpha[39] * agri24[i, j] * ZoneBFC[i] +

        alpha[40] * distHighway[i] + alpha[41] * distHighway[i] * ZoneAlps[i] +
                                     alpha[42] * distHighway[i] * ZoneJura[i] +
                                     alpha[43] * distHighway[i] * ZoneVosges[i] +
                                     alpha[44] * distHighway[i] * ZoneBFC[i] +

        alpha[45] * distHighway2[i] + alpha[46] * distHighway2[i] * ZoneAlps[i] +
                                      alpha[47] * distHighway2[i] * ZoneJura[i] +
                                      alpha[48] * distHighway2[i] * ZoneVosges[i] +
                                      alpha[49] * distHighway2[i] * ZoneBFC[i] +

        alpha[50] * roadLength[i, j] + alpha[51] * roadLength[i, j] * ZoneAlps[i] +
                                       alpha[52] * roadLength[i, j] * ZoneJura[i] +
                                       alpha[53] * roadLength[i, j] * ZoneVosges[i] +
                                       alpha[54] * roadLength[i, j] * ZoneBFC[i] +

        alpha[55] * pathLength[i, j] + alpha[56] * pathLength[i, j] * ZoneAlps[i] +
                                       alpha[57] * pathLength[i, j] * ZoneJura[i] +
                                       alpha[58] * pathLength[i, j] * ZoneVosges[i] +
                                       alpha[59] * pathLength[i, j] * ZoneBFC[i] +

        alpha[60] * hDens[i] + alpha[61] * hDens[i] * ZoneAlps[i] +
                               alpha[62] * hDens[i] * ZoneJura[i] +
                               alpha[63] * hDens[i] * ZoneVosges[i] +
                               alpha[64] * hDens[i] * ZoneBFC[i] +

        alpha[65] * tri[i] + alpha[66] * tri[i] * ZoneAlps[i] +
                             alpha[67] * tri[i] * ZoneJura[i] +
                             alpha[68] * tri[i] * ZoneVosges[i] +
                             alpha[69] * tri[i] * ZoneBFC[i] +

        alpha[70] * hfi[i, j] + alpha[71] * hfi[i, j] * ZoneAlps[i] +
                             alpha[72] * hfi[i, j] * ZoneJura[i] +
                             alpha[73] * hfi[i, j] * ZoneVosges[i] +
                             alpha[74] * hfi[i, j] * ZoneBFC[i] +
        b[j] + u[i]
  
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
init1 <- list(alpha = runif(74, -1, 1),
              beta = runif(2, -2, 2),
              z = zst,
              b = rep(0, nYears),
              mu.b = 0,
              u = rep(0, nSites),
              v = matrix(ncol = nYears, nrow = nSites, data = rep(0, nSites * nYears)),
              sd.b = 1,
              sd.u = 1,
              sd.v = 1)
init2 <- list(alpha = runif(74, -1, 1),
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
lynxSim <- jags(data = dataFull,
                inits = inits,
                parameters = parameters,
                n.iter = 20000,
                model.file = "modelOcc.txt",
                n.chains = 2,
                n.burnin = 10000)
save(lynxSim, file="outputs/lynxMCP100_byMassif.RData")
load("outputs/lynxMCP100_byMassif.RData")

#####################
# Model predictions #
#####################

# Maps of occupancy probability based on the habitat only

# Function to recreate occupancy probabilities
occEstimate <- function(matrixAlpa, vectorBeta, lengthCov, 
                        coveZoneA, covZoneJ, covZoneV, covZoneB,
                        covConnF, covShrub, covOpenL, covA21, covA22, covA23, covA24, 
                        covDistH, covDistH2, covRoadLength, covPathLength, covHDens, 
                        covTri, covHfi){
  
  mapsValues <- matrix(nrow = lengthCov, ncol = nrow(matrixAlpa))
  
  for(i in 1:nrow(matrixAlpa)){
    
    logitVal <- 
      
      matrixAlpa[i, 1] * coveZoneA + matrixAlpa[i, 2] * covZoneJ + 
      matrixAlpa[i, 3] * covZoneV + matrixAlpa[i, 4] * covZoneB +
      
      matrixAlpa[i, 5] * covConnF + 
      matrixAlpa[i, 6] * covConnF* coveZoneA +
      matrixAlpa[i, 7] * covConnF * covZoneJ+
      matrixAlpa[i, 8] * covConnF * covZoneV +
      matrixAlpa[i, 9] * covConnF * covZoneB +
      
      matrixAlpa[i, 10] * covShrub + 
      matrixAlpa[i, 11] * covShrub * coveZoneA +
      matrixAlpa[i, 12] * covShrub * covZoneJ +
      matrixAlpa[i, 13] * covShrub * covZoneV +
      matrixAlpa[i, 14] * covShrub * covZoneB +
      
      matrixAlpa[i, 15] * covOpenL + 
      matrixAlpa[i, 16] * covOpenL * coveZoneA +
      matrixAlpa[i, 17] * covOpenL * covZoneJ +
      matrixAlpa[i, 18] * covOpenL * covZoneV +
      matrixAlpa[i, 19] * covOpenL * covZoneB +
      
      matrixAlpa[i, 20] * covA21 + 
      matrixAlpa[i, 21] * covA21 * coveZoneA +
      matrixAlpa[i, 22] * covA21 * covZoneJ +
      matrixAlpa[i, 23] * covA21 * covZoneV +
      matrixAlpa[i, 24] * covA21 * covZoneB +
      
      matrixAlpa[i, 25] * covA22 + 
      matrixAlpa[i, 26] * covA22 * coveZoneA +
      matrixAlpa[i, 27] * covA22 * covZoneJ +
      matrixAlpa[i, 28] * covA22 * covZoneV +
      matrixAlpa[i, 29] * covA22 * covZoneB +
      
      matrixAlpa[i, 30] * covA23 +
      matrixAlpa[i, 31] * covA23 * coveZoneA +
      matrixAlpa[i, 32] * covA23 * covZoneJ +
      matrixAlpa[i, 33] * covA23 * covZoneV +
      matrixAlpa[i, 34] * covA23 * covZoneB +
      
      matrixAlpa[i, 35] * covA24 + 
      matrixAlpa[i, 36] * covA24 * coveZoneA +
      matrixAlpa[i, 37] * covA24 * covZoneJ +
      matrixAlpa[i, 38] * covA24 * covZoneV +
      matrixAlpa[i, 39] * covA24 * covZoneB +
      
      matrixAlpa[i, 40] * covDistH + 
      matrixAlpa[i, 41] * covDistH * coveZoneA +
      matrixAlpa[i, 42] * covDistH * covZoneJ +
      matrixAlpa[i, 43] * covDistH * covZoneV +
      matrixAlpa[i, 44] * covDistH * covZoneB +
      
      matrixAlpa[i, 45] * covDistH2 +
      matrixAlpa[i, 46] * covDistH2 * coveZoneA +
      matrixAlpa[i, 47] * covDistH2 * covZoneJ +
      matrixAlpa[i, 48] * covDistH2 * covZoneV +
      matrixAlpa[i, 49] * covDistH2 * covZoneB +
      
      matrixAlpa[i, 50] * covRoadLength + 
      matrixAlpa[i, 51] * covRoadLength * coveZoneA +
      matrixAlpa[i, 52] * covRoadLength * covZoneJ +
      matrixAlpa[i, 53] * covRoadLength * covZoneV +
      matrixAlpa[i, 54] * covRoadLength * covZoneB +
      
      matrixAlpa[i, 55] * covPathLength +
      matrixAlpa[i, 56] * covPathLength * coveZoneA +
      matrixAlpa[i, 57] * covPathLength * covZoneJ +
      matrixAlpa[i, 58] * covPathLength * covZoneV +
      matrixAlpa[i, 59] * covPathLength * covZoneB +
      
      matrixAlpa[i, 60] * covHDens +
      matrixAlpa[i, 61] * covHDens * coveZoneA +
      matrixAlpa[i, 62] * covHDens * covZoneJ +
      matrixAlpa[i, 63] * covHDens * covZoneV +
      matrixAlpa[i, 64] * covHDens * covZoneB +
      
      matrixAlpa[i, 65] * covTri + 
      matrixAlpa[i, 66] * covTri * coveZoneA +
      matrixAlpa[i, 67] * covTri * covZoneJ +
      matrixAlpa[i, 68] * covTri * covZoneV +
      matrixAlpa[i, 69] * covTri * covZoneB +
      
      matrixAlpa[i, 70] * covHfi + 
      matrixAlpa[i, 71] * covHfi * coveZoneA +
      matrixAlpa[i, 72] * covHfi * covZoneJ +
      matrixAlpa[i, 73] * covHfi * covZoneV +
      matrixAlpa[i, 74] * covHfi * covZoneB +
      
      vectorBeta[i]
    
    mapsValues[, i] <- exp(logitVal) / (1 + exp(logitVal))
    #print(i)
  }
  
  return(mapsValues)
}

alphaEstim <- lynxSim$sims.list$alpha
bEstim <- lynxSim$sims.list$b[, 27]

# Predictions of the effect of each covariate as a function of the area

## Distance to highways
covVal <- as.numeric(allDistH)
covValForScaling <- as.numeric(covSampl[, "distHgwCov[, -1]"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
covValForScaling2 <- as.numeric(covSampl[, "distHgwCov2[, -1]"])
covValS2 <- ((covValTest * covValTest) - mean(covValForScaling2)) / sd(covValForScaling2)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                      covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = covValS, covDistH2 = covValS2, covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                      covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 1), lwd = 2,
     xlab = "Distance to highways", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                         coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                         covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                         covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                         covDistH = covValS, covDistH2 = covValS2, covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                         covHDens = mean(allDensHS), covTri = mean(allTriS), 
                         covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Tri
# Define the range of value for the covariate
covVal <- as.numeric(allTri)
covValForScaling <- as.numeric(covSampl[, "triCov[, -1]"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                       coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                       covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                       covA21 = mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                       covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), 
                       covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                       covHDens = mean(allDensHS), covTri = covValS, 
                       covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 1), lwd = 2,
     xlab = "Terrain Ruggeness Index", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                         coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                         covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                         covA21 = mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                         covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), 
                         covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                         covHDens = mean(allDensHS), covTri = covValS, 
                         covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Shrub
# Define the range of value for the covariate
covVal <- as.numeric(allShrub)
covValForScaling <- as.numeric(covSampl[, "shrub2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                       coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                       covConnF = mean(allConnFS), covShrub = covValS, covOpenL = mean(allOpenLS), 
                       covA21 = mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                       covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                       covHDens = mean(allDensHS), covTri = mean(allTriS), 
                       covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 1), lwd = 2,
     xlab = "Shrub cover", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                         coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                         covConnF = mean(allConnFS), covShrub = covValS, covOpenL = mean(allOpenLS), 
                         covA21 = mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                         covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                         covHDens = mean(allDensHS), covTri = mean(allTriS), 
                         covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Agri21
covVal <- as.numeric(allA21)
covValForScaling <- as.numeric(covSampl[, "agri21_2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                       coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                       covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                       covA21 = covValS, covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                       covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                       covHDens = mean(allDensHS), covTri = mean(allTriS), 
                       covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 0.8), lwd = 2,
     xlab = "Agri21 cover", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                         coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                         covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                         covA21 = covValS, covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                         covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                         covHDens = mean(allDensHS), covTri = mean(allTriS), 
                         covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Agri22
covVal <- as.numeric(allA22)
covValForScaling <- as.numeric(covSampl[, "agri22_2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                       coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                       covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                       covA21 =  mean(allA21S), covA22 = covValS, covA23 = mean(allA23S), covA24 = mean(allA24S), 
                       covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                       covHDens = mean(allDensHS), covTri = mean(allTriS), 
                       covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 1), lwd = 2,
     xlab = "Agri22 cover", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                         coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                         covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                         covA21 =  mean(allA21S), covA22 = covValS, covA23 = mean(allA23S), covA24 = mean(allA24S), 
                         covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                         covHDens = mean(allDensHS), covTri = mean(allTriS), 
                         covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Agri23
covVal <- as.numeric(allA23)
covValForScaling <- as.numeric(covSampl[, "agri23_2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                      covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                       covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = covValS, covA24 = mean(allA24S), 
                       covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                       covHDens = mean(allDensHS), covTri = mean(allTriS), 
                       covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 0.7), lwd = 2,
     xlab = "Agri23 cover", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                         coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                         covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                         covA21 =  mean(allA21S), covA22 = mean(allA23S), covA23 = covValS, covA24 = mean(allA24S), 
                         covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                         covHDens = mean(allDensHS), covTri = mean(allTriS), 
                         covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Agri24
covVal <- as.numeric(allA24)
covValForScaling <- as.numeric(covSampl[, "agri24_2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                      covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = covValS, 
                      covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                      covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 0.8), lwd = 2,
     xlab = "Agri24 cover", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                        coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                        covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                        covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = covValS, 
                        covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                        covHDens = mean(allDensHS), covTri = mean(allTriS), 
                        covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Open land
covVal <- as.numeric(allOpenL)
covValForScaling <- as.numeric(covSampl[, "agri24_2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                      covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = covValS, 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                      covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 1), lwd = 2,
     xlab = "Open land", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                        coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                        covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = covValS, 
                        covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                        covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                        covHDens = mean(allDensHS), covTri = mean(allTriS), 
                        covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Human density
covVal <- as.numeric(allDensH)
covValForScaling <- as.numeric(covSampl[, "hDensCov[, -1]"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                      covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                      covHDens = covValS, covTri = mean(allTriS), 
                      covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 1), lwd = 2,
     xlab = "Human density", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                        coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                        covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                        covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                        covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                        covHDens = covValS, covTri = mean(allTriS), 
                        covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Human footprint
covVal <- as.numeric(allHfi)
covValForScaling <- as.numeric(covSampl[, "hDensCov[, -1]"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                      covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                      covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covHfi = covValS)
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 0.7), lwd = 2,
     xlab = "Human footprint", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                        coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                        covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                        covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                        covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                        covHDens = mean(allDensHS), covTri = mean(allTriS), 
                        covHfi = covValS)
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Forest connectivity
covVal <- as.numeric(allConnF)
covValForScaling <- as.numeric(covSampl[, "connectFor18"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                      covConnF = covValS, covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                      covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 0.8), lwd = 2,
     xlab = "Forest connectivity", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                        coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                        covConnF = covValS, covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                        covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                        covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = mean(allPathLS),
                        covHDens = mean(allDensHS), covTri = mean(allTriS), 
                        covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Road length
covVal <- as.numeric(allRdL)
covValForScaling <- as.numeric(covSampl[, "roadLength2018"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                      covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = covValS, covPathLength = mean(allPathLS),
                      covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 1), lwd = 2,
     xlab = "Road length", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                        coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                        covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                        covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                        covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = covValS, covPathLength = mean(allPathLS),
                        covHDens = mean(allDensHS), covTri = mean(allTriS), 
                        covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}

## Path length
covVal <- as.numeric(allPathL)
covValForScaling <- as.numeric(covSampl[, "pathLength2023"])
rangeCov <- c(min(covVal, na.rm = TRUE), max(covVal, na.rm = TRUE))
covValTest <- seq(from = rangeCov[1], to = rangeCov[2], length.out = 100) # test 100 values for the covariates in the natural range in the landscape
covValS <- (covValTest - mean(covValForScaling)) / sd(covValForScaling)
# Estimate occupancy for this range of values for this covariate with all the others at their mean value
# ALPS
valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                      coveZoneA = 1, covZoneJ = 0, covZoneV = 0, covZoneB = 0,
                      covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                      covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                      covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = covValS,
                      covHDens = mean(allDensHS), covTri = mean(allTriS), 
                      covHfi = mean(allHfiS))
meanOcc <- rowMeans(valFor, na.rm = TRUE)
quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
# Compute the mean and CI occupancy probabilities
plot(x = covValTest, y = meanOcc, type = "l", ylim = c(0, 1), lwd = 2,
     xlab = "Path length", ylab = "Occupancy probability", col = "red")
lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = "red")
lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = "red")
legend("topleft",legend = c("Alps", "Jura", "Vosges", "BFC", "Massif C."), col = c("red", "blue", "green", "black", "orange"), lwd = c(2, 2, 2, 2, 2))

for(i in 1:4){
  zoneA <- 0
  zoneJ <- 0
  zoneV <- 0
  zoneB <- 0
  if(i == 1){
    zoneJ <- 1
    colZone <- "blue"
  } else if(i == 2){
    zoneV <- 1
    colZone <- "green"
  } else if(i == 3){
    zoneB <- 1
    colZone <- "black"
  } else if(i == 4){
    colZone <- "orange"
  }
  
  valFor <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = 100,
                        coveZoneA = zoneA, covZoneJ = zoneJ, covZoneV = zoneV, covZoneB = zoneB,
                        covConnF = mean(allConnFS), covShrub = mean(allShrubS), covOpenL = mean(allOpenLS), 
                        covA21 =  mean(allA21S), covA22 = mean(allA22S), covA23 = mean(allA23S), covA24 = mean(allA24S), 
                        covDistH = mean(allDistHS), covDistH2 = mean(allDistHS2), covRoadLength = mean(allRdLS), covPathLength = covValS,
                        covHDens = mean(allDensHS), covTri = mean(allTriS), 
                        covHfi = mean(allHfiS))
  meanOcc <- rowMeans(valFor, na.rm = TRUE)
  quantileOcc <- apply(valFor, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
  # Compute the mean and CI occupancy probabilities
  lines(x = covValTest, y = meanOcc, type = "l", lwd = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[1,], lty = 2, col = colZone)
  lines(x = covValTest, y = quantileOcc[2,], lty = 2, col = colZone)
}


###

# Plot predicted occupancy probability on the map
allMapsValues <- occEstimate(matrixAlpa = alphaEstim, vectorBeta = bEstim, lengthCov = length(allForS), 
                             coveZoneA = gridFullAlps, covZoneJ = gridFullJura, covZoneV = gridFullVosges, 
                             covZoneB = gridFullBFC,
                             covConnF = allConnFS, covShrub = allShrubS, covOpenL = allOpenLS, 
                             covA21 = allA21S, covA22 = allA22S, covA23 = allA23S, covA24 = allA24S, covDistH = allDistHS, 
                             covDistH2 = allDistHS2, covRoadLength = allRdLS, covPathLength = allPathLS,
                             covHDens = allDensHS, covTri = allTriS, covHfi = allHfiS)

# Compute the mean and quantiles for each cell over all iterations
meanOcc <- rowMeans(allMapsValues, na.rm = TRUE)
quantileOcc <- apply(allMapsValues, 1, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))

mapPotentialOcc <- covSF %>% 
  mutate(meanOcc = meanOcc, lowCIocc = quantileOcc[1, ], upCIocc = quantileOcc[2, ])

# Plot predicted occupancy probability with data points
plotWithoutPoints <- ggplot() + 
  geom_sf(data = mapPotentialOcc, aes(fill = meanOcc), color = NA) +
  scale_fill_viridis_c(name = "Occupancy prob.") +
  ggtitle("Mean occupancy predicted")
print(plotWithoutPoints)

plotWithPoints <- ggplot() + 
  geom_sf(data = mapPotentialOcc, aes(fill = meanOcc), color = NA) + 
  scale_fill_viridis_c(name = "Occupancy prob.") +
  geom_sf(data = dataP, shape = 4, color = alpha("black", 0.2)) +
  ggtitle("Mean occupancy predicted")
print(plotWithPoints)
