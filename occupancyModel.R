setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_lynxOccupancy/")
# Appendix Bauduin, ... and Gimenez

library(raster)
library(lubridate)

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

covScale <- scale(cov) # scaling by column

# Grid cell on which lynx presence data are identified
load("data/gridFrComplete.RData")
# Crop the grid to the study area extent (east of France)
gridFrComplete <- crop(gridFrComplete, extent(c(720000, 1090000, 6260000, 6920000)))

# Lynx presence data
dataSCALP <- read.csv("C:/Users/sarah.bauduin/Documents/SaveOFB/Projets/8.SCALP/a classifier/2019-2020/codesSCALP.csv",
                      header = TRUE, sep = ";")
dataSCALP <- dataSCALP[!is.na(dataSCALP$X),]
dataSCALP <- dataSCALP[!is.na(dataSCALP$Y),]
# Add the coordinates system
lynxData_spdf <- SpatialPointsDataFrame(coords = dataSCALP[,c("X", "Y")], data = dataSCALP, 
                                         proj4string = CRS(as.character("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")),
                                         match.ID = TRUE, bbox = NULL)
lynxData_spdf$DateObservation <- as.Date(lynxData_spdf$DateObservation)
lynxData_spdfTr <- spTransform(lynxData_spdf, gridFrComplete@proj4string)

# Which cells are occupied during which year
yearStart <- min(unique(lynxData_spdfTr$Annee))
yearEnd <- max(unique(lynxData_spdfTr$Annee)) - 1
obsLynx <- array(data = 0, dim = c(length(gridFrComplete$ID), 6, yearEnd - yearStart + 1))

for(year in yearStart:yearEnd){
  
  # Extract the lynx observations for the different time period
  novYear <- lynxData_spdfTr[lynxData_spdfTr$DateObservation >= as.Date(paste(year, "11", "01", sep = "-")) &
                               lynxData_spdfTr$DateObservation < (as.Date(paste(year, "11", "01", sep = "-")))%m+% months(1),]
  decYear <- lynxData_spdfTr[lynxData_spdfTr$DateObservation >= as.Date(paste(year, "12", "01", sep = "-")) &
                               lynxData_spdfTr$DateObservation < (as.Date(paste(year, "12", "01", sep = "-")))%m+% months(1),]
  janYearPlus1 <- lynxData_spdfTr[lynxData_spdfTr$DateObservation >= as.Date(paste(year + 1, "01", "01", sep = "-")) &
                                    lynxData_spdfTr$DateObservation < (as.Date(paste(year + 1, "01", "01", sep = "-")))%m+% months(1),]
  febYearPlus1 <- lynxData_spdfTr[lynxData_spdfTr$DateObservation >= as.Date(paste(year + 1, "02", "01", sep = "-")) &
                                    lynxData_spdfTr$DateObservation < (as.Date(paste(year + 1, "02", "01", sep = "-")))%m+% months(1),]
  marYearPlus1 <- lynxData_spdfTr[lynxData_spdfTr$DateObservation >= as.Date(paste(year + 1, "03", "01", sep = "-")) &
                                    lynxData_spdfTr$DateObservation < (as.Date(paste(year + 1, "03", "01", sep = "-")))%m+% months(1),]
  aprYearPlus1 <- lynxData_spdfTr[lynxData_spdfTr$DateObservation >= as.Date(paste(year + 1, "04", "01", sep = "-")) &
                                    lynxData_spdfTr$DateObservation < (as.Date(paste(year + 1, "04", "01", sep = "-")))%m+% months(1),]
  
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


# Change format of obsLynx
obsLynx2 <- obsLynx[, , 1]
for (i in 2:dim(obsLynx)[3]){
  obsLynx2 <- cbind(obsLynx2, obsLynx[, , i])
}
dim(obsLynx2) # 2131 (n cells in the grid) x (6 months x yearStart to yearEnd)

# Which cell were never sampled across all years
noSampl <- apply(obsLynx2, 1, function(x) all(is.na(x)))
sum(noSampl)
# Remove cells never samples
covScaleSampl <- covScale[!noSampl,]
samplEffortSampl <- samplEffort[!noSampl,]
obsLynx2Sampl <- obsLynx2[!noSampl,]



