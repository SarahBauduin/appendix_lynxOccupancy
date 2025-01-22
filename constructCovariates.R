setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_lynxOccupancy/")
# Appendix Bauduin, ... and Gimenez

library(raster)
library(sf)
library(dplyr)
library(terra)
library(spatialEco)
library(osmextract)

# # Covariates to run the lynx occupancy model
# # Grid cell on which covariates must be computed
# load("data/gridFrComplete.RData")
# # Crop the grid to the study area extent (east of France)
# gridFrComplete <- crop(gridFrComplete, extent(c(720000, 1090000, 6260000, 6920000)))
# # Remove cells that are not completely in France
# load("data/franceShape.RData")
# gridOnlyFr <- raster::intersect(gridFrComplete, franceShape)
# gridOnlyFr$area <- gArea(gridOnlyFr, byid = TRUE)
# gridFrComplete <- gridOnlyFr[gridOnlyFr$area >= 1e+08, ]
# #save(gridFrComplete, file = "data/gridFrCompleteCrop.RData")
load("data/gridFrCompleteCrop.RData")
gridFrComplete <- st_as_sf(gridFrComplete)


####################
####################
## LANDSCAPE DATA ##
####################
####################


#####################
## LANDSCAPE COVER ##
#####################

# CLC data from 1990, 2000, 2006, 2012 and 2018
clc1990 <- st_read("data/CLC/CLC1990/CLC90_FR_RGF.shp")
clc2000 <- st_read("data/CLC/CLC2000/CLC00_FR_RGF.shp")
clc2006 <- st_read("data/CLC/CLC2006/CLC06_FR_RGF.shp")
clc2012 <- st_read("data/CLC/CLC2012/CLC12_FR_RGF.shp")
clc2018 <- st_read("data/CLC/CLC2018/clc2018.shp")
clc2018 <- clc2018[,c(1,3)] # keep the CLC code as the 2nd column
listCLC <- list()
listCLC[[1]] <- clc1990
listCLC[[2]] <- clc2000
listCLC[[3]] <- clc2006
listCLC[[4]] <- clc2012
listCLC[[5]] <- clc2018

# Tranform the gridded study area in a sf format 
# with the same coordinates as the CLC
gridFr_sf <- gridFrComplete %>%
  st_transform(crs = st_crs(clc1990))

# Function to extract the proportion of certain CLC covers (given by the "codesCLC")
# for the 5 CLC files (different year)
extractPropCLC <- function(codesCLC){ # codesCLC as a vector of characters
  
  coverCov <- data.frame(ID = gridFrComplete$ID)

  for(i in 1:length(listCLC)){
    
    # Extract only the desired CLC covers
    coverYear <- listCLC[[i]] %>%
      filter(.[[2]] %in% codesCLC)
    # Compute the proportion of these covers
    cellCoverYear <- st_intersection(gridFr_sf , coverYear) %>%
      mutate(areaInter = st_area(.)) %>%
      group_by(ID) %>%
      summarise(areaCoverCell = sum(areaInter)) %>%
      mutate(propCoverYear = areaCoverCell / 1e+08) %>%
      as_tibble() %>%
      select(ID, propCoverYear)
    
    # Merge the results to the main df with the cell ID
    # Rename the column according to the CLC year
    coverCov <- merge(coverCov, as.data.frame(cellCoverYear), all = TRUE)
    if(i == 1){
      colnames(coverCov)[i+1] <- "propCover1990"
    } else if(i == 2){
      colnames(coverCov)[i+1] <- "propCover2000"
    } else if(i == 3){
      colnames(coverCov)[i+1] <- "propCover2006"
    } else if(i == 4){
      colnames(coverCov)[i+1] <- "propCover2012"
    } else if(i == 5){
      colnames(coverCov)[i+1] <- "propCover2018"
    }
  }
  coverCov[is.na(coverCov)] <- 0
  return(coverCov)
}


#############################
## Proportion of forest cover
forestCov <- extractPropCLC(codesCLC = c("311", "312", "313"))
colnames(forestCov) <- c("ID", "forest1990", "forest2000", "forest2006", "forest2012", "forest2018")


######################
## Forest connectivity 
# The number of "forested cells" (i.e., cells with more than 50% of forest cover),
# in a 10 km buffer
forest <- list()
forest[[1]] <- clc1990 %>%
  filter(CODE_90 %in% c("311", "312", "313"))
forest[[2]] <- clc2000 %>%
  filter(CODE_00 %in% c("311", "312", "313"))
forest[[3]] <- clc2006 %>%
  filter(CODE_06 %in% c("311", "312", "313"))
forest[[4]] <- clc2012%>%
  filter(CODE_12 %in% c("311", "312", "313"))
forest[[5]] <- clc2018%>%
  filter(code_18 %in% c("311", "312", "313"))

gridFr_connectFor <- gridFr_sf

for(i in 1:length(forest)){
  
  # Which cells have 50% of more forest cover
  gridForest <- st_intersection(gridFr_sf, forest[[i]])
  gridForest$area <- st_area(gridForest)
  # Forested cells
  gridForest2 <- gridForest %>% 
    group_by(ID) %>%
    summarise(area = sum(area)) %>% 
    mutate(forest = ifelse(as.numeric(area) >= 50000000, 1, 0))
  
  # 9 km buffer around all cells
  gridForest3 <- gridForest2 %>% 
    st_buffer(dist = 9000) %>% 
    rename(IDBuffer = ID)
  
  # How many in buffers
  cellsInBuffer <- gridForest3 %>% 
    st_intersection(gridFr_sf) %>% 
    group_by(IDBuffer) %>%
    summarise(nCells = n())

  # Forested cells in the buffers
  gridForest4 <- gridForest2 %>% 
    st_intersection(gridForest3) %>% 
    group_by(IDBuffer) %>% 
    summarise(sumForest1 = sum(forest))
  gridForest4 <- gridForest4 %>% 
    mutate(sumForest = sumForest1 / cellsInBuffer$nCells)
  
  gridForest5 <- gridFr_sf %>% 
    left_join(as.data.frame(gridForest4)[,c("IDBuffer", "sumForest")], by = c("ID" = "IDBuffer"))
  gridForest5[is.na(gridForest5$sumForest), "sumForest"] <- 0
  
  # Put back the data in the grid
  if(i == 1){
    gridFr_connectFor <- gridFr_connectFor %>% 
      mutate(connectFor90 = gridForest5$sumForest)
  }
  if(i == 2){
    gridFr_connectFor <- gridFr_connectFor %>% 
      mutate(connectFor00 = gridForest5$sumForest)
  }
  if(i == 3){
    gridFr_connectFor <- gridFr_connectFor %>% 
      mutate(connectFor06 = gridForest5$sumForest)
  }
  if(i == 4){
    gridFr_connectFor <- gridFr_connectFor %>% 
      mutate(connectFor12 = gridForest5$sumForest)
  }
  if(i == 5){
    gridFr_connectFor <- gridFr_connectFor %>% 
      mutate(connectFor18 = gridForest5$sumForest)
  }
  print(i)
  
}

connectForCov <- gridFr_connectFor %>% 
  st_drop_geometry() %>% 
  as.data.frame() %>% 
  select(-1)


############################
## Proportion of shrub cover 
shrubCov <- extractPropCLC(codesCLC = c("321", "322", "323", "324"))
colnames(shrubCov) <- c("ID", "shrub1990", "shrub2000", "shrub2006", "shrub2012", "shrub2018")


################################
## Proportion of open land cover 
openLandCov <- extractPropCLC(codesCLC = c("331", "332", "333", "334", "335"))
colnames(openLandCov) <- c("ID", "openLand1990", "openLand2000", "openLand2006", "openLand2012", "openLand2018")


##############################
## Proportion of agri 21 cover
agri21Cov <- extractPropCLC(codesCLC = c("211", "213"))
colnames(agri21Cov) <- c("ID", "agri21_1990", "agri21_2000", "agri21_2006", "agri21_2012", "agri21_2018")


##############################
## Proportion of agri 22 cover
agri22Cov <- extractPropCLC(codesCLC = c("221", "222", "223"))
colnames(agri22Cov) <- c("ID", "agri22_1990", "agri22_2000", "agri22_2006", "agri22_2012", "agri22_2018")


############
## Proportion of agri 23 cover
agri23Cov <- extractPropCLC(codesCLC = "231")
colnames(agri23Cov) <- c("ID", "agri23_1990", "agri23_2000", "agri23_2006", "agri23_2012", "agri23_2018")


############
## Proportion of agri 24 cover
agri24Cov <- extractPropCLC(codesCLC = c("241", "242", "243", "244"))
colnames(agri24Cov) <- c("ID", "agri24_1990", "agri24_2000", "agri24_2006", "agri24_2012", "agri24_2018")



###############
## ELEVATION ##
###############

# Elevation data on 25 m x 25 m resolution
elev1 <- rast("data/elevation/032ab314564b9cb72c98fbeb093aeaf69720fbfd/eu_dem_v11_E30N20.TIF") 
gridFrCompleteTr <- gridFrComplete %>%
  st_transform(crs = st_crs(elev1))
elev1 <- terra::crop(elev1, st_sf(st_union(gridFrCompleteTr)))
elev2 <- rast("data/elevation/97824c12f357f50638d665b5a58707cd82857d57/eu_dem_v11_E40N20.TIF")
elev2 <- terra::crop(elev2, st_sf(st_union(gridFrCompleteTr)))
elevation <- terra::merge(elev1, elev2)


############
## Elevation

# Compute mean and sd elevation on 100 km2 cells
elevMean <- elevation %>% 
  terra::extract(gridFrCompleteTr) %>% 
  group_by(ID) %>% 
  summarise(meanElev = mean(Band_1, na.rm = TRUE), sdElev = sd(Band_1, na.rm = TRUE))
elevCov <- as.data.frame(elevMean)


#############
## Ruggedness

# Compute TRI (terrain ruggedness index) on 1 km2 cells
elevAggr <- terra::aggregate(elevation, fact = 40, fun = "mean", na.rm = TRUE)
elevTRI <- tri(elevAggr)

# Compute the number of 1 km2 "rugged cells" into 100 km2 cells
elevTRI[elevTRI < 162] <- 0 # 162 = "intermediately rugged surface"
elevTRI[elevTRI >= 162] <- 0.01
elevTRI_100km2 <- aggregate(elevTRI, fact = 10, fun = sum, na.rm = TRUE)
triCells <- elevTRI_100km2 %>% 
  terra::extract(., st_centroid(gridFrCompleteTr))
triCov <- cbind.data.frame(ID = gridFrCompleteTr$ID, tri = triCells[ ,2])

# Other index of ruggedness: VRM (Vector Ruggedness Measure; Sappington 2007)
# Compute VRM on 1 km2 cells (same resolution as TRI)
elevVRM <- vrm(elevAggr, s = c(3, 3)) # use of a 3x3 window to be similar to TRI

# Compute the number of 1 km2 "rugged cells" into 100 km2 cells
elevVRM[elevVRM < 0.005] <- 0 # 0.005 = medium rugged surface as mean values of mountains = 0.01 (Sappington)
elevVRM[elevVRM >= 0.005] <- 0.01
elevVRM_100km2 <- aggregate(elevVRM, fact = 10, fun = sum, na.rm = TRUE)
vrmCells <- elevVRM_100km2 %>% 
  terra::extract(., st_centroid(gridFrCompleteTr))
vrmCov <- cbind.data.frame(ID = gridFrCompleteTr$ID, vrm005 = vrmCells[ ,2])

# Trying with another (lower) threshold to obtain a more heterogeous surface
elevVRM[elevVRM < 0.0001] <- 0 
elevVRM[elevVRM >= 0.0001] <- 0.01
elevVRM_100km2 <- aggregate(elevVRM, fact = 10, fun = sum, na.rm = TRUE)
vrmCells <- elevVRM_100km2 %>% 
  terra::extract(., st_centroid(gridFrCompleteTr))
vrmCov <- cbind.data.frame(vrmCov, vrm0001 = vrmCells[ ,2])


###############
## PREY DATA ##
###############

#########
# Chamois
chamData <- st_read("data/ungulates/chamois/chamois_effectif.shp")
gridFrCompleteTr <- st_transform(gridFrComplete, st_crs(chamData))

# Estimate chamois density per m2
chamData <- chamData %>% 
  mutate(area = st_area(.)) %>% 
  mutate(densCham = effectif_m / area) %>% 
  st_intersection(gridFrCompleteTr)

# Sum the density inside each grid cell
cham1993 <- chamData %>% 
  filter(annee_effe == 1993) %>% 
  group_by(ID) %>% 
  summarise(dens1993 = sum(densCham, na.rm = TRUE))

cham2005 <- chamData %>% 
  filter(annee_effe == 2005) %>% 
  group_by(ID) %>% 
  summarise(dens2005 = sum(densCham, na.rm = TRUE))

cham2010 <- chamData %>% 
  filter(annee_effe == 2010) %>% 
  group_by(ID) %>% 
  summarise(dens2010 = sum(densCham, na.rm = TRUE))

# Replace NA by 0 (no chamois there)
gridFrCompleteTr <- merge(gridFrCompleteTr, st_drop_geometry(cham1993[,c("ID", "dens1993")]), all = TRUE) %>% 
  mutate(dens1993 = ifelse(is.na(dens1993), 0, dens1993))
gridFrCompleteTr <- merge(gridFrCompleteTr, st_drop_geometry(cham2005[,c("ID", "dens2005")]), all = TRUE) %>% 
  mutate(dens2005 = ifelse(is.na(dens2005), 0, dens2005))
gridFrCompleteTr <- merge(gridFrCompleteTr, st_drop_geometry(cham2010[,c("ID", "dens2010")]), all = TRUE) %>% 
  mutate(dens2010 = ifelse(is.na(dens2010), 0, dens2010))

# Add spatial autocorrelation by averaging the 9 cells in each 3x3 cells square
gridFrCompleteRef <- gridFrCompleteTr
gridFrCompleteBuff <- gridFrCompleteTr[,"ID"] %>% 
  st_buffer(dist = 10000) %>% 
  st_intersection(gridFrCompleteRef) %>% 
  group_by(ID) %>% 
  summarise(dens1993AS = mean(dens1993), dens2005AS = mean(dens2005), dens2010AS = mean(dens2010))

chamCov <- cbind.data.frame(ID = gridFrCompleteBuff$ID, 
                            cham1993 = gridFrCompleteBuff$dens1993AS,
                            cham2005 = gridFrCompleteBuff$dens2005AS,
                            cham2010 = gridFrCompleteBuff$dens2010AS)


######
# Deer

# Deer data (hunting bags) per french districts
deerDistrictsData <- read.csv("data/ungulates/deer/enq_comm_chevreuil_chamois_oct_2024.csv", header = TRUE, sep = ",")
deerDistrictsData <- deerDistrictsData[deerDistrictsData$espece == "Chevreuil                ",]
# French districts
frDistricts <- st_read("data/ungulates/deer/communes-20220101-shp/communes-20220101.shp")
# Merge the deer data with the district shapefile based on the insee number for each year
for(year in unique(deerDistrictsData$annee)){
  deerYear <- deerDistrictsData[deerDistrictsData$annee == year, c("insee", "realisations")]
  colnames(deerYear) <- c("insee", paste0("realisations_", year))
  frDistricts <- merge(frDistricts, deerYear[,1:2], all = TRUE)
}
gridFrCompleteTr <- st_transform(gridFrComplete, st_crs(frDistricts))

# Divide the deer data into two time periods
deerData <- frDistricts %>% 
  mutate(meanAllYears = rowMeans(st_drop_geometry(frDistricts[,5:27]), na.rm = TRUE)) %>% 
  mutate(meanFirstPeriod = rowMeans(st_drop_geometry(frDistricts[,6:11]), na.rm = TRUE)) %>%  # 1993-2007
  mutate(meanScndPeriod = rowMeans(st_drop_geometry(frDistricts[,12:27]), na.rm = TRUE)) %>%  # 2012-2023
  st_crop(., gridFrCompleteTr) %>% 
  mutate(area = st_area(.)) %>% 
  mutate(densAllYears = meanAllYears / area) %>% 
  mutate(densFP = meanFirstPeriod / area) %>% 
  mutate(densSP = meanScndPeriod / area) 

# Grid intersected with the french districts
disGrid <- st_read("data/ungulates/deer/grid_comm_intersect/intersection_commShp2_gridFrComplete2.shp")
disGridFull <- disGrid %>% 
  merge(st_drop_geometry(deerData[,c(1, 33:35)]), all = TRUE) %>% 
  mutate(newArea = st_area(.)) %>% 
  mutate(partDensAllY = densAllYears * newArea) %>% 
  mutate(partDensFP = densFP * newArea) %>% 
  mutate(partDensSP = densSP * newArea) %>% 
  group_by(ID) %>% 
  summarise(meanDensAllYPerCell = sum(partDensAllY), meanDensFPPerCell = sum(partDensFP), meanDensSPPerCell = sum(partDensSP), sumArea = sum(newArea))

# For incomplete cell (i.e., less than 90000000 for the area - missing 10%) or those with NA, we take the mean
ID10percent <- st_drop_geometry(disGridFull[as.numeric(disGridFull$sumArea) < 90000000, "ID"])
IDNA_allY <- unique(c(ID10percent$ID, st_drop_geometry(disGridFull[is.na(disGridFull$meanDensAllYPerCell), "ID"]$ID)))
IDNA_FP <-  unique(c(ID10percent$ID, st_drop_geometry(disGridFull[is.na(disGridFull$meanDensFPPerCell), "ID"]$ID)))
IDNA_SP <-  unique(c(ID10percent$ID, st_drop_geometry(disGridFull[is.na(disGridFull$meanDensSPPerCell), "ID"]$ID)))
missingAllY <- disGrid %>% 
  merge(st_drop_geometry(deerData[,c(1, 33)]), all = TRUE) %>% 
  filter(ID %in% IDNA_allY[!is.na(IDNA_allY)]) %>% 
  group_by(ID) %>% 
  summarise(meanDensAllYPerCell = (mean(densAllYears, na.rm = TRUE)) * 100000000)
missingFP <- disGrid %>% 
  merge(st_drop_geometry(deerData[,c(1, 34)]), all = TRUE) %>% 
  filter(ID %in% IDNA_FP[!is.na(IDNA_FP)]) %>% 
  group_by(ID) %>% 
  summarise(meanDensFPPerCell = (mean(densFP, na.rm = TRUE)) * 100000000)
missingSP <- disGrid %>% 
  merge(st_drop_geometry(deerData[,c(1, 35)]), all = TRUE) %>% 
  filter(ID %in% IDNA_SP[!is.na(IDNA_SP)]) %>% 
  group_by(ID) %>% 
  summarise(meanDensSPPerCell = (mean(densSP, na.rm = TRUE)) * 100000000)

# Put these new estimates
disGridFull <- disGridFull %>% 
  mutate(meanDensAllYPerCell = ifelse(ID %in% IDNA_allY[!is.na(IDNA_allY)], missingAllY$meanDensAllYPerCell, meanDensAllYPerCell)) %>% 
  mutate(meanDensFPPerCell = ifelse(ID %in% IDNA_FP[!is.na(IDNA_FP)], missingFP$meanDensFPPerCell, meanDensFPPerCell)) %>% 
  mutate(meanDensSPPerCell = ifelse(ID %in% IDNA_SP[!is.na(IDNA_SP)], missingSP$meanDensSPPerCell, meanDensSPPerCell))

# Remove the few outliers
tail(disGridFull$meanDensAllYPerCell[order(disGridFull$meanDensAllYPerCell)])
disGridFull[!is.na(disGridFull$meanDensAllYPerCell) & disGridFull$meanDensAllYPerCell > 407, "meanDensAllYPerCell"] <- 407
tail(disGridFull$meanDensFPPerCell[order(disGridFull$meanDensFPPerCell)][!is.na(disGridFull$meanDensFPPerCell[order(disGridFull$meanDensFPPerCell)])])
disGridFull[!is.na(disGridFull$meanDensFPPerCell) & disGridFull$meanDensFPPerCell > 409, "meanDensFPPerCell"] <- 409
tail(disGridFull$meanDensSPPerCell[order(disGridFull$meanDensSPPerCell)][!is.na(disGridFull$meanDensSPPerCell[order(disGridFull$meanDensSPPerCell)])])
disGridFull[!is.na(disGridFull$meanDensSPPerCell) & disGridFull$meanDensSPPerCell > 426, "meanDensSPPerCell"] <- 426

# Reassign to the original grid
gridDeer <- merge(gridFrComplete, st_drop_geometry(disGridFull), all = TRUE)
gridDeer <- gridDeer %>% 
  filter(!is.na(ID))

# Add spatial autocorrelation by averaging the 9 cells in each 3x3 cells square
gridFrRef <- gridFrComplete
gridFrCompleteBuff <- gridFrComplete[,"ID"] %>% 
  st_buffer(dist = 10000) %>% 
  st_intersection(gridDeer[,c("meanDensAllYPerCell", "meanDensFPPerCell", "meanDensSPPerCell")]) %>% 
  group_by(ID) %>% 
  summarise(deerAllYears = mean(meanDensAllYPerCell, na.rm = TRUE), deerFirstPeriod = mean(meanDensFPPerCell, na.rm = TRUE), deerSecondPeriod = mean(meanDensSPPerCell, na.rm = TRUE))

deerCov <- cbind.data.frame(ID = gridFrCompleteBuff$ID, 
                            deerAllY = gridFrCompleteBuff$deerAllYears,
                            deerFP = gridFrCompleteBuff$deerFirstPeriod,
                            deerSP = gridFrCompleteBuff$deerSecondPeriod)



##################
## Save covariates
save(forestCov, connectForCov, shrubCov, openLandCov, agri21Cov, agri22Cov, agri23Cov, agri24Cov, 
     elevCov, triCov, vrmCov,
     chamCov, deerCov,
     file = "outputs/covLandscape.RData")

##########################################################



#################
#################
## LINEAR DATA ##
#################
#################

# Route500 data from 2012, 2015, and 2018 for the main roads
roads2012 <- st_read("data/route500/2012/TRONCON_ROUTE.shp")
roads2015 <- st_read("data/route500/2015/TRONCON_ROUTE.shp")
roads2018 <- st_read("data/route500/2018/TRONCON_ROUTE.shp")
# OpenStreetMap from 2014, 2019, 2023 for the paths
france2014 <- oe_read("data/openStreetMap/france-140101.osm.pbf")
france2019 <- oe_read("data/openStreetMap/france-190101.osm.pbf")
france2023 <- oe_read("data/openStreetMap/france-230101.osm.pbf")

# Highways from 2014
highways <- st_read("data/highways/cropHighways.shp") 
# not to the exact extent of the grid to calculate closest highway which may be outside of the grid extent

# Rivers
rivers <- st_read("data/river/TronconHydrograElt_FXX.shp/TronconHydrograElt_FXX.shp")

#######################
## Distance to highways
highwaysTr <- highways %>%
  st_transform(crs = st_crs(gridFrComplete)) %>% 
  st_union(.)

# Cell centroids
gridCentroid <- st_centroid(gridFrComplete)
# Distance for each centroid
distHgwsCells <- st_distance(gridCentroid, highwaysTr)
distHgwCov <- cbind.data.frame(ID = gridFrComplete$ID, distHgws = as.numeric(distHgwsCells))


#####################
## Distance to rivers
riversTr <- rivers %>%
  st_transform(crs = st_crs(gridFrComplete)) %>% 
  # Removing dried and underground rivers
  filter(Etat %in% c("En attente de mise à jour", "Inconnu", "Intermittent", "Permanent")) %>% 
  st_union(.)
# Distance for each centroid
distRiverCells <- st_distance(gridCentroid, riversTr)
distRiverCov <- cbind.data.frame(ID = gridFrComplete$ID, distRiver = as.numeric(distRiverCells))


##############
## Road length 
roads2012Tr <- roads2012 %>%
  st_transform(crs = st_crs(gridFrComplete)) %>% 
  st_crop(., gridFrComplete) %>% 
  st_union(.)
roads2015Tr <- roads2015 %>%
  st_transform(crs = st_crs(gridFrComplete)) %>% 
  st_crop(., gridFrComplete) %>% 
  st_union(.)
roads2018Tr <- roads2018 %>%
  st_transform(crs = st_crs(gridFrComplete)) %>% 
  st_crop(., gridFrComplete) %>% 
  st_union(.)
listRoads <- list()
listRoads[[1]] <- roads2012Tr
listRoads[[2]] <- roads2015Tr
listRoads[[3]] <- roads2018Tr

# Total road length in each cell
roadLengthCov <- data.frame(ID = gridFrComplete$ID)

for(i in 1:length(listRoads)){
  
  roadsIn <- st_intersection(gridFrComplete, listRoads[[i]]) %>% # extract all the roads inside each cell
    mutate(x = st_length(.)) # calculate the length of each type of road per cell ID
  lengthRoadsCells <- st_drop_geometry(roadsIn)[,c("ID","x")]
  lengthRoadsCells$x <- as.numeric(lengthRoadsCells$x)
  
  # Not all cells had roads in it
  withoutRoads <- setdiff(gridFrComplete$ID, lengthRoadsCells$ID) # cells without roads in it
  lengthRoadsID <- c(lengthRoadsCells$ID, withoutRoads)
  lengthRoads <- c(lengthRoadsCells$x, rep(0, length(withoutRoads))) # put 0 as the length of roads in the empty cells
  lengthRoadsOrder <- cbind(lengthRoadsID, lengthRoads)
  lengthRoadsOrderID <- lengthRoadsOrder[order(match(lengthRoadsOrder[,"lengthRoadsID"] ,gridFrComplete$ID)),]
  
  roadLengthCov <- cbind.data.frame(roadLengthCov, as.data.frame(lengthRoadsOrderID)[,"lengthRoads"])
  if(i == 1){
    colnames(roadLengthCov)[i+1] <- "roadLength2012"
  } else if(i == 2){
    colnames(roadLengthCov)[i+1] <- "roadLength2015"
  } else if(i == 3){
    colnames(roadLengthCov)[i+1] <- "roadLength2018"
  }
  print(i)
}


##############
## Path length 
path2014Tr <- france2014 %>%
  st_transform(crs = st_crs(gridFrComplete)) %>% 
  filter(highway %in% c("track", "footway", "bridleway", "path", "cylceway")) %>% 
  st_crop(., gridFrComplete) %>% 
  st_union(.)
path2019Tr <- france2019 %>%
  st_transform(crs = st_crs(gridFrComplete)) %>% 
  filter(highway %in% c("track", "footway", "bridleway", "path", "cylceway")) %>% 
  st_crop(., gridFrComplete) %>% 
  st_union(.)
path2023Tr <- france2023 %>%
  st_transform(crs = st_crs(gridFrComplete)) %>% 
  filter(highway %in% c("track", "footway", "bridleway", "path", "cylceway")) %>% 
  st_crop(., gridFrComplete) %>% 
  st_union(.)
listPaths <- list()
listPaths[[1]] <- path2014Tr
listPaths[[2]] <- path2019Tr
listPaths[[3]] <- path2023Tr

# Total path length in each cell
pathLengthCov <- data.frame(ID = gridFrComplete$ID)

for(i in 1:length(listPaths)){
  
  pathIn <- st_intersection(gridFrComplete, listPaths[[i]]) %>% # extract all the path inside each cell
    mutate(x = st_length(.)) # calculate the length of the paths per cell ID
  lengthPathCells <- st_drop_geometry(pathIn)[,c("ID","x")]
  lengthPathCells$x <- as.numeric(lengthPathCells$x)
  
  # Not all cells had path in it
  withoutPath <- setdiff(gridFrComplete$ID, lengthPathCells$ID) # cells without path in it
  lengthPathID <- c(lengthPathCells$ID, withoutPath)
  lengthPaths <- c(lengthPathCells$x, rep(0, length(withoutPath))) # put 0 as the length of path in the empty cells
  lengthPathsOrder <- cbind(lengthPathID, lengthPaths)
  lengthPathsOrderID <- lengthPathsOrder[order(match(lengthPathsOrder[,"lengthPathID"] ,gridFrComplete$ID)),]
  
  pathLengthCov <- cbind.data.frame(pathLengthCov, as.data.frame(lengthPathsOrderID)[,"lengthPaths"])
  if(i == 1){
    colnames(pathLengthCov)[i+1] <- "pathLength2014"
  } else if(i == 2){
    colnames(pathLengthCov)[i+1] <- "pathLength2019"
  } else if(i == 3){
    colnames(pathLengthCov)[i+1] <- "pathLength2023"
  }
  print(i)
}


##################
## Save covariates
save(distHgwCov, distRiverCov, roadLengthCov, pathLengthCov, 
     file = "outputs/covLinear.RData")

###############################################################



########################
########################
## HUMAN DISTURBANCES ##
########################
########################


################
## Human density

# Human density for cells of 1km2 in 2017
human_sf <- st_read("data/humanPop/2017/Filosofi2017_carreaux_1km_met.shp")

gridFr_sf <- gridFrComplete %>%
  st_transform(crs = st_crs(human_sf))

# Compute the mean human density (1km2) into the 100km2 cells
humanDens <- st_intersection(gridFr_sf , human_sf) %>%
  group_by(ID) %>%
  summarise(humanDens1km2 = mean(Ind, na.rm = TRUE)) %>%
  as_tibble() %>%
  dplyr::select(ID, humanDens1km2)
# Fill with 0 the missing data
hDensCov <- data.frame(ID = gridFrComplete$ID)
hDensCov <- merge(hDensCov, as.data.frame(humanDens), all = TRUE)
hDensCov[is.na(hDensCov)] <- 0


#######################
# Human Footprint Index

hfi_1993 <- rast("data/HFI/1993/wildareas-v3-1993-human-footprint.tif")
hfi_2009 <- rast("data/HFI/2009/wildareas-v3-2009-human-footprint.tif")
# 50 is the maximum value on land, then 128 is given for water areas
hfi_1993[hfi_1993 == 128] <- 25 # gives a medium value
hfi_2009[hfi_2009 == 128] <- 25

gridFr_sf <- gridFrComplete %>%
  st_transform(crs = crs(hfi_1993))
# Mean HFI per grid cells
hfiCells_1993 <- terra::extract(hfi_1993, gridFr_sf, fun = mean)
hfiCells_2009 <- terra::extract(hfi_2009, gridFr_sf, fun = mean)

hfiCov <- cbind.data.frame(ID = gridFrComplete$ID, 
                           hfi_1993 = hfiCells_1993[, 2],
                           hfi_2009 = hfiCells_2009[, 2])


############
# Nightlight
light <- rast("data/nightlight/VNL_npp_2023_global_vcmslcfg_v2_c202402081600.average.dat.tif/VNL_npp_2023_global_vcmslcfg_v2_c202402081600.average.dat.tif")
gridFr_sf <- gridFrComplete %>%
  st_transform(crs = crs(light))

# Crop the nightlight data to the grid extent
lightFr <- crop(light, extent(gridFr_sf))
# Mean light per grid cells
lightCells <- terra::extract(lightFr, gridFr_sf, fun = mean)

lightCov <- cbind.data.frame(ID = gridFrComplete$ID, 
                             light = lightCells[, 2])


########
# Strava
# Strava records outdoor activities
strava <- rast("data/strava/heatmap_strava_20220414_100_all.tiff")
# 4 layers
# layer 1 = all
# layer 2 = run
# layer 3 = ride
# layer 4 = winter

gridFrCompleteTr <- gridFrComplete %>%
  st_transform(crs = st_crs(strava))

stravaMean <- strava %>% 
  terra::extract(gridFrCompleteTr) %>% 
  mutate(strava3cat = heatmap_strava_20220414_100_all_2 + heatmap_strava_20220414_100_all_3 + 
           heatmap_strava_20220414_100_all_4) %>% 
  group_by(ID) %>% 
  summarise(meanStravaAll = mean(heatmap_strava_20220414_100_all_1), meanStrava3cat = mean(strava3cat),
            meanStravaRun = mean(heatmap_strava_20220414_100_all_2), meanStravaRide = mean(heatmap_strava_20220414_100_all_3),
            meanStravaWinter = mean(heatmap_strava_20220414_100_all_4))

stravaCov <- as.data.frame(stravaMean)


################
# Outdoor vision
# Another type of data recording outdoor activities

outVisCov


##################
## Save covariates
save(hDensCov, hfiCov, lightCov, stravaCov, outVisCov,
     file = "outputs/covHuman.RData")

########################################################################


