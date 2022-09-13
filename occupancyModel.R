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



