setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_lynxOccupancy/")
# Appendix Bauduin, ... and Gimenez

library(raster)
library(sf)
library(dplyr)
library(terra)
library(spatialEco)
library(osmextract)

# Covariates to run the lynx occupancy model

# # Grid cell on which covariates must be computed
# load("data/gridFrComplete.RData")
# # Crop the grid to the study area extent (east of France)
# gridFrComplete <- crop(gridFrComplete, extent(c(720000, 1090000, 6260000, 6920000)))
# 
# # Create a grid of 1km2
# grid1 <- st_make_grid(st_as_sf(gridFrComplete), cellsize = 1000, square = TRUE)
# # Remove cells that are not completely in France
# load("data/franceShape.RData")
# grid1 <- grid1 %>%
#   st_sf(.) %>%
#   mutate(ID = 1:nrow(.)) %>%
#   st_intersection(., st_as_sf(franceShape)) %>%
#   mutate(area = st_area(.)) %>%
#   subset(as.numeric(area) == 1000000)
# save(grid1, file = "data/grid1.RData")
load("data/grid1.RData")


##########################
## LANDSCAPE COVER DATA ##
##########################
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
gridFr_sf <- grid1 %>%
  st_transform(crs = st_crs(clc1990))

# Function to extract the proportion of certain CLC covers (given by the "codesCLC")
# for the 5 CLC files (different year)
extractPropCLC <- function(codesCLC){ # codesCLC as a vector of characters
  
  coverCov <- data.frame(ID = grid1$ID)

  for(i in 1:length(listCLC)){
    
    # Extract only the desired CLC covers
    coverYear <- listCLC[[i]] %>%
      filter(.[[2]] %in% codesCLC)
    # Compute the proportion of these covers
    cellCoverYear <- st_intersection(gridFr_sf , coverYear) %>%
      mutate(areaInter = st_area(.)) %>%
      group_by(ID) %>%
      summarise(areaCoverCell = sum(areaInter)) %>%
      mutate(propCoverYear = areaCoverCell / 1e+06) %>%
      as_tibble() %>%
      dplyr::select(ID, propCoverYear)
    
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
forest[[4]] <- clc2012 %>%
  filter(CODE_12 %in% c("311", "312", "313"))
forest[[5]] <- clc2018 %>%
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
    mutate(forest = ifelse(as.numeric(area) >= 500000, 1, 0))
  
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
  dplyr::select(-1)

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

##################
## Save covariates
save(forestCov, connectForCov, shrubCov, openLandCov, agri21Cov, agri22Cov, 
     agri23Cov, agri24Cov, file = "outputs/covariatesCover1km.RData")

##########################################################


###############
## ROAD DATA ##
###############
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
  st_transform(crs = st_crs(grid1)) %>% 
  st_union(.)

# Cell centroids
gridCentroid <- st_centroid(grid1)
# Distance for each centroid
distHgwsCells <- st_distance(gridCentroid, highwaysTr)
distHgwCov <- cbind.data.frame(ID = grid1$ID, distHgws = as.numeric(distHgwsCells))

#####################
## Distance to rivers
riversTr <- rivers %>%
  st_transform(crs = st_crs(grid1)) %>% 
  # Removing dried and underground rivers
  filter(Etat %in% c("En attente de mise à jour", "Inconnu", "Intermittent", "Permanent")) %>% 
  st_union(.)
# Distance for each centroid
distRiverCells <- st_distance(gridCentroid, riversTr)
distRiverCov <- cbind.data.frame(ID = grid1$ID, distRiver = as.numeric(distRiverCells))

##############
## Road length 
roads2012Tr <- roads2012 %>%
  st_transform(crs = st_crs(grid1)) %>% 
  st_crop(., grid1) %>% 
  st_union(.)
roads2015Tr <- roads2015 %>%
  st_transform(crs = st_crs(grid1)) %>% 
  st_crop(., grid1) %>% 
  st_union(.)
roads2018Tr <- roads2018 %>%
  st_transform(crs = st_crs(grid1)) %>% 
  st_crop(., grid1) %>% 
  st_union(.)
listRoads <- list()
listRoads[[1]] <- roads2012Tr
listRoads[[2]] <- roads2015Tr
listRoads[[3]] <- roads2018Tr

# Total road length in each cell
roadLengthCov <- data.frame(ID = grid1$ID)

for(i in 1:length(listRoads)){

  roadsIn <- st_intersection(grid1, listRoads[[i]]) %>% # extract all the roads inside each cell
    mutate(x = st_length(.)) # calculate the length of each type of road per cell ID
  lengthRoadsCells <- st_drop_geometry(roadsIn)[,c("ID","x")]
  lengthRoadsCells$x <- as.numeric(lengthRoadsCells$x)
  
  # Not all cells had roads in it
  withoutRoads <- setdiff(grid1$ID, lengthRoadsCells$ID) # cells without roads in it
  lengthRoadsID <- c(lengthRoadsCells$ID, withoutRoads)
  lengthRoads <- c(lengthRoadsCells$x, rep(0, length(withoutRoads))) # put 0 as the length of roads in the empty cells
  lengthRoadsOrder <- cbind(lengthRoadsID, lengthRoads)
  lengthRoadsOrderID <- lengthRoadsOrder[order(match(lengthRoadsOrder[,"lengthRoadsID"] ,grid1$ID)),]

  roadLengthCov <- cbind.data.frame(roadLengthCov, as.data.frame(lengthRoadsOrderID)[,"lengthRoads"])
  if(i == 1){
    colnames(roadLengthCov)[i+1] <- "roadLength2012"
  } else if(i == 2){
    colnames(roadLengthCov)[i+1] <- "roadLength2015"
  } else if(i == 3){
    colnames(roadLengthCov)[i+1] <- "roadLength2018"
  }
}

##############
## Path length 
path2014Tr <- france2014 %>%
  st_transform(crs = st_crs(grid1)) %>% 
  filter(highway %in% c("track", "footway", "bridleway", "path", "cylceway")) %>% 
  st_crop(., grid1) %>% 
  st_union(.)
path2019Tr <- france2019 %>%
  st_transform(crs = st_crs(grid1)) %>% 
  filter(highway %in% c("track", "footway", "bridleway", "path", "cylceway")) %>% 
  st_crop(., grid1) %>% 
  st_union(.)
path2023Tr <- france2023 %>%
  st_transform(crs = st_crs(grid1)) %>% 
  filter(highway %in% c("track", "footway", "bridleway", "path", "cylceway")) %>% 
  st_crop(., grid1) %>% 
  st_union(.)
listPaths <- list()
listPaths[[1]] <- path2014Tr
listPaths[[2]] <- path2019Tr
listPaths[[3]] <- path2023Tr

# Total path length in each cell
pathLengthCov <- data.frame(ID = grid1$ID)

for(i in 1:length(listPaths)){
  
  pathIn <- st_intersection(grid1, listPaths[[i]]) %>% # extract all the path inside each cell
    mutate(x = st_length(.)) # calculate the length of the paths per cell ID
  lengthPathCells <- st_drop_geometry(pathIn)[,c("ID","x")]
  lengthPathCells$x <- as.numeric(lengthPathCells$x)
  
  # Not all cells had path in it
  withoutPath <- setdiff(grid1$ID, lengthPathCells$ID) # cells without path in it
  lengthPathID <- c(lengthPathCells$ID, withoutPath)
  lengthPaths <- c(lengthPathCells$x, rep(0, length(withoutPath))) # put 0 as the length of path in the empty cells
  lengthPathsOrder <- cbind(lengthPathID, lengthPaths)
  lengthPathsOrderID <- lengthPathsOrder[order(match(lengthPathsOrder[,"lengthPathID"] ,grid1$ID)),]
  
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
     file = "outputs/covariatesRoads.RData")

###############################################################


###########################
## HUMAN POPULATION DATA ##
###########################
# Human density for cells of 1km2 in 2017
human_sf <- st_read("data/humanPop/2017/Filosofi2017_carreaux_1km_met.shp")

gridFr_sf <- grid1 %>%
  st_transform(crs = st_crs(human_sf))

################
## Human density
# Compute the mean human density in 1km2
humanDens <- gridFr_sf %>% 
  st_join(human_sf) %>% 
  group_by(ID) %>%
  summarise(humanDens1km2 = mean(Ind, na.rm = TRUE))
humanDens[is.nan(humanDens$humanDens1km2),"humanDens1km2"] <- NA

humanDensTiblble <- humanDens %>%
  as_tibble() %>%
  dplyr::select(ID, humanDens1km2)

# There is no data on every cells
# Compute the mean from the surrounding cells
# Identify the 8 neighboring cells
buffer <- gridFr_sf %>% 
  st_buffer(dist = 900) %>% # one cell buffer (a bit less)
  rename(IDbuffer = ID) %>% 
  dplyr::select(IDbuffer)

bufferCells <- st_intersection(buffer, humanDens)
humanDensCov <- merge(as.data.frame(humanDensTiblble), data.frame(ID = grid1$ID), all = TRUE)

# As long as there are NA, fill the cells
while(sum(is.na(humanDensCov$humanDens1km2)) != 0){
  # Cells which are NA
  cellYearNAID <- humanDensCov %>% 
    filter(is.na(.[, "humanDens1km2"])) %>% 
    pull(ID)
  
  # Buffer of cells which are NA
  bufferYearNA <- bufferCells %>% 
    dplyr::select(IDbuffer, humanDens1km2) %>%
    filter(IDbuffer %in% cellYearNAID) %>% 
    st_drop_geometry() %>% 
    group_by(IDbuffer) %>% 
    summarise(humanDens1km2 := mean(humanDens1km2, na.rm = TRUE))
  
  # Replace the mean value in the NA cells
  humanDensCov[humanDensCov$ID %in% bufferYearNA$IDbuffer, "humanDens1km2"] <- pull(bufferYearNA, 2)
  mergeBuff <- merge(humanDensCov[, c("ID", "humanDens1km2")], data.frame(ID = pull(bufferCells, "ID")), all.y = TRUE)
  newValue <- mergeBuff[match(pull(bufferCells, "ID"), mergeBuff$ID), "humanDens1km2"]
  bufferCells[, "humanDens1km2"] <- newValue
  
  print(sum(is.na(humanDensCov$humanDens1km2))) 
}

hDensCov <- humanDensCov


####################
## ELEVATION DATA ##
####################
# Elevation data on 25 m x 25 m resolution
elev1 <- rast("data/elevation/032ab314564b9cb72c98fbeb093aeaf69720fbfd/eu_dem_v11_E30N20.TIF")# 
grid1Tr <- grid1 %>%
  st_transform(crs = st_crs(elev1))
elev1 <- terra::crop(elev1, st_sf(st_union(grid1Tr)))
elev2 <- rast("data/elevation/97824c12f357f50638d665b5a58707cd82857d57/eu_dem_v11_E40N20.TIF")
elev2 <- terra::crop(elev2, st_sf(st_union(grid1Tr)))
elevation <- terra::merge(elev1, elev2)

# Compute mean and sd elevation on 100 km2 cells
elevMean <- elevation %>% 
  terra::extract(grid1Tr) %>% 
  group_by(ID) %>% 
  summarise(meanElev = mean(Band_1), sdElev = sd(Band_1))
elevCov <- as.data.frame(elevMean)


#############
## Ruggedness
# Compute TRI (terrain ruggedness index) on 1km2 cells
elevAggr <- terra::aggregate(elevation, fact = 40, fun = "mean", na.rm = TRUE)
elevTRI <- tri(elevAggr)

# Define "rugged cells"
elevTRI[elevTRI < 162] <- 0 # 162 = "intermediately rugged surface"
elevTRI[elevTRI >= 162] <- 1
triCells <- elevTRI %>% 
  terra::extract(., st_centroid(grid1Tr))
triCov <- cbind.data.frame(ID = grid1Tr$ID, tri = triCells[, 2])


################
## TALLY DATA ##
################
# Chamois data
cha1985 <- st_read("data/ungulates/CHA_1985_departement_L93_1698223650_4441/CHA_1985_departement_L93_1698223650_4441/CHA_1985_departement_L93.shp")
cha1990 <- st_read("data/ungulates/CHA_1991_departement_L93_1698224637_3806/CHA_1991_departement_L93_1698224637_3806/CHA_1991_departement_L93.shp")
cha1995 <- st_read("data/ungulates/CHA_1995_departement_L93_1698223650_3646/CHA_1995_departement_L93_1698223650_3646/CHA_1995_departement_L93.shp")
cha2000 <- st_read("data/ungulates/CHA_2000_departement_L93_1698223650_1913/CHA_2000_departement_L93_1698223650_1913/CHA_2000_departement_L93.shp")
cha2005 <- st_read("data/ungulates/CHA_2005_departement_L93_1698223650_3713/CHA_2005_departement_L93_1698223650_3713/CHA_2005_departement_L93.shp")
cha2010 <- st_read("data/ungulates/CHA_2010_departement_L93_1698223650_2960/CHA_2010_departement_L93_1698223650_2960/CHA_2010_departement_L93.shp")
cha2015 <- st_read("data/ungulates/CHA_2015_departement_L93_1698223650_3241/CHA_2015_departement_L93_1698223650_3241/CHA_2015_departement_L93.shp")
cha2018 <- st_read("data/ungulates/CHA_2018_departement_L93_1698235389_8043/CHA_2018_departement_L93_1698235389_8043/CHA_2018_departement_L93.shp")
# Deer data
deer1985 <- st_read("data/ungulates/CHE_1985_departement_L93_1698223584_8629/CHE_1985_departement_L93_1698223584_8629/CHE_1985_departement_L93.shp")
deer1990 <- st_read("data/ungulates/CHE_1990_departement_L93_1698223583_9203/CHE_1990_departement_L93_1698223583_9203/CHE_1990_departement_L93.shp")
deer1995 <- st_read("data/ungulates/CHE_1995_departement_L93_1698223583_1133/CHE_1995_departement_L93_1698223583_1133/CHE_1995_departement_L93.shp")
deer2000 <- st_read("data/ungulates/CHE_2000_departement_L93_1698223583_4856/CHE_2000_departement_L93_1698223583_4856/CHE_2000_departement_L93.shp")
deer2005 <- st_read("data/ungulates/CHE_2005_departement_L93_1698223583_3019/CHE_2005_departement_L93_1698223583_3019/CHE_2005_departement_L93.shp")
deer2010 <- st_read("data/ungulates/CHE_2010_departement_L93_1698223583_2115/CHE_2010_departement_L93_1698223583_2115/CHE_2010_departement_L93.shp")
deer2015 <- st_read("data/ungulates/CHE_2015_departement_L93_1698223583_4561/CHE_2015_departement_L93_1698223583_4561/CHE_2015_departement_L93.shp")
deer2018 <- st_read("data/ungulates/CHE_2018_departement_L93_1698235371_2084/CHE_2018_departement_L93_1698235371_2084/CHE_2018_departement_L93.shp")

####################
## Prey availability
grid1Tr <- grid1 %>%
  st_transform(crs = st_crs(cha1985))

# Compute preys density per km2
preyDensFun <- function(dataUngulates, nameColYear){
  dataUngulates <- dataUngulates %>% 
    mutate(area = st_area(.)) %>% 
    mutate(preyDens = Realstn / area) %>% 
    st_intersection(grid1Tr) %>% 
    mutate(areaInter = st_area(.)) %>% 
    mutate(preyDensInter = preyDens * (areaInter / area.1)) %>% # ratio of the cell
    group_by(ID) %>% 
    summarise(preyDensMean = sum(preyDensInter)) %>% 
    select(ID, preyDensMean) %>% 
    st_drop_geometry() %>% 
    as.data.frame()
  # Fill with 0 the missing data
  ungCov <- data.frame(ID = grid1$ID)
  ungCov <- merge(ungCov, dataUngulates, all = TRUE)
  ungCov$preyDensMean <- as.numeric(ungCov$preyDensMean)
  ungCov[is.na(ungCov)] <- 0
  colnames(ungCov)[2] <- nameColYear
  return(ungCov)
}

# Chamois
cha1985Cov <- preyDensFun(dataUngulates = cha1985, nameColYear = "cha1985")
cha1990Cov <- preyDensFun(dataUngulates = cha1990, nameColYear = "cha1990")
cha1995Cov <- preyDensFun(dataUngulates = cha1995, nameColYear = "cha1995")
cha2000Cov <- preyDensFun(dataUngulates = cha2000, nameColYear = "cha2000")
cha2005Cov <- preyDensFun(dataUngulates = cha2005, nameColYear = "cha2005")
cha2010Cov <- preyDensFun(dataUngulates = cha2010, nameColYear = "cha2010")
cha2015Cov <- preyDensFun(dataUngulates = cha2015, nameColYear = "cha2015")
cha2018Cov <- preyDensFun(dataUngulates = cha2018, nameColYear = "cha2018")
chaCov <- cbind.data.frame(cha1985Cov, cha1990Cov = cha1990Cov[, 2], cha1995Cov = cha1995Cov[, 2], 
                           cha2000Cov = cha2000Cov[, 2], cha2005Cov = cha2005Cov[, 2], cha2010Cov = cha2010Cov[, 2], 
                           cha2015Cov = cha2015Cov[, 2], cha2018Cov = cha2018Cov[, 2])
# Deer
deer1985Cov <- preyDensFun(dataUngulates = deer1985, nameColYear = "deer1985")
deer1990Cov <- preyDensFun(dataUngulates = deer1990, nameColYear = "deer1990")
deer1995Cov <- preyDensFun(dataUngulates = deer1995, nameColYear = "deer1995")
deer2000Cov <- preyDensFun(dataUngulates = deer2000, nameColYear = "deer2000")
deer2005Cov <- preyDensFun(dataUngulates = deer2005, nameColYear = "deer2005")
deer2010Cov <- preyDensFun(dataUngulates = deer2010, nameColYear = "deer2010")
deer2015Cov <- preyDensFun(dataUngulates = deer2015, nameColYear = "deer2015")
deer2018Cov <- preyDensFun(dataUngulates = deer2018, nameColYear = "deer2018")
deerCov <- cbind.data.frame(deer1985Cov, deer1990Cov = deer1990Cov[, 2], deer1995Cov = deer1995Cov[, 2], 
                            deer2000Cov = deer2000Cov[, 2], deer2005Cov = deer2005Cov[, 2], deer2010Cov = deer2010Cov[, 2], 
                            deer2015Cov = deer2015Cov[, 2], deer2018Cov = deer2018Cov[, 2])


###########################
## HUMAN FOOTPRINT INDEX ##
###########################
# HFI data
hfi_1993 <- rast("data/HFI/1993/wildareas-v3-1993-human-footprint.tif")
hfi_2009 <- rast("data/HFI/2009/wildareas-v3-2009-human-footprint.tif")
# 50 is the maximum value on land, then 128 is given for water areas
hfi_1993[hfi_1993 == 128] <- 25 # gives a medium value
hfi_2009[hfi_2009 == 128] <- 25

####################
## Human disturbance 
gridFr_sf <- grid1 %>%
  st_transform(crs = crs(hfi_1993))
# HFI per grid cells
hfiCells_1993 <- hfi_1993 %>% 
  terra::extract(., st_centroid(gridFr_sf))
hfiCells_2009 <- hfi_2009 %>% 
  terra::extract(., st_centroid(gridFr_sf))

hfiCov <- cbind.data.frame(ID = grid1$ID, 
                           hfi_1993 = hfiCells_1993[, 2],
                           hfi_2009 = hfiCells_2009[, 2])


############
## STRAVA ##
############
# Strava
strava <- rast("data/strava/heatmap_strava_20220414_100_all.tiff")
# 4 layers
# layer 1 = all
# layer 2 = run
# layer 3 = ride
# layer 4 = winter
grid1Tr <- grid1 %>%
  st_transform(crs = st_crs(strava))
stravaMean <- strava %>% 
  terra::extract(grid1Tr) %>% 
  mutate(strava3cat = heatmap_strava_20220414_100_all_2 + heatmap_strava_20220414_100_all_3 + 
           heatmap_strava_20220414_100_all_4) %>% 
  group_by(ID) %>% 
  summarise(meanStravaAll = mean(heatmap_strava_20220414_100_all_1), meanStrava3cat = mean(strava3cat),
            meanStravaRun = mean(heatmap_strava_20220414_100_all_2), meanStravaRide = mean(heatmap_strava_20220414_100_all_3),
            meanStravaWinter = mean(heatmap_strava_20220414_100_all_4))
stravaCov <- as.data.frame(stravaMean)

##################
## Save covariates
save(hDensCov, elevCov, triCov, chaCov, deerCov, hfiCov, stravaCov,
     file = "outputs/covariatesOthers1km.RData")

########################################################################


