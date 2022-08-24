setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_lynxOccupancy/")
# Appendix Bauduin, ... and Gimenez

library(raster)
library(sf)
library(dplyr)
library(rgeos)

# Covariates to run the lynx occupancy model
# Grid cell on which covariates must be computed
load("data/gridFrComplete.RData")
# Crop the grid to the study area extent (east of France)
gridFrComplete <- crop(gridFrComplete, extent(c(720000, 1090000, 6260000, 6920000)))

##########################
## Landscape cover data ##
##########################
# CLC data from 1990, 2000, 2006, 2012 and 2018
clc1990 <- shapefile("data/CLC/CLC1990/CLC90_FR_RGF.shp")
clc2000 <- shapefile("data/CLC/CLC2000/CLC00_FR_RGF.shp")
clc2006 <- shapefile("data/CLC/CLC2006/CLC06_FR_RGF.shp")
clc2012 <- shapefile("data/CLC/CLC2012/CLC12_FR_RGF.shp")
clc2018 <- shapefile("data/CLC/CLC2018/clc2018.shp")
clc2018 <- clc2018[,c(1,3)] # keep the CLC code as the 2nd column
listCLC <- list()
listCLC[[1]] <- clc1990
listCLC[[2]] <- clc2000
listCLC[[3]] <- clc2006
listCLC[[4]] <- clc2012
listCLC[[5]] <- clc2018

# Tranform the gridded study area in a sf format 
# with the same coordinates as the CLC
gridFr_sf <- st_as_sf(gridFrComplete) %>%
  st_transform(crs = st_crs(st_as_sf(clc1990)))

# Function to extract the proportion of certain CLC covers (given by the "codesCLC")
# for the 5 CLC files (different year)
extractPropCLC <- function(codesCLC){ # codesCLC as a vector of characters
  
  coverCov <- data.frame(ID = gridFrComplete$ID)

  for(i in 1:length(listCLC)){
    
    # Extract only the desired CLC covers
    coverYear <- st_as_sf(listCLC[[i]]) %>%
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


######################
## Forest connectivity 
# From "forested cells" (i.e., cells with more than 50% of forest cover),
# the number of "forested cells" in a 10 km buffer
forest <- list()
forest[[1]] <- st_as_sf(clc1990) %>%
  filter(CODE_90 %in% c("311", "312", "313"))
forest[[2]] <- st_as_sf(clc2000) %>%
  filter(CODE_00 %in% c("311", "312", "313"))
forest[[3]] <- st_as_sf(clc2006) %>%
  filter(CODE_06 %in% c("311", "312", "313"))
forest[[4]] <- st_as_sf(clc2012)%>%
  filter(CODE_12 %in% c("311", "312", "313"))
forest[[5]] <- st_as_sf(clc2018)%>%
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
  
  # 10 km buffer around these cells
  gridForest3 <- gridForest2 %>% 
    filter(forest == 1) %>% 
    st_buffer(dist = 10000) %>% 
    rename(IDBuffer = ID)
  
  # Forested cells in the buffers
  gridForest4 <- gridForest2 %>% 
    filter(forest == 1) %>%
    st_intersection(gridForest3) %>% 
    group_by(IDBuffer) %>% 
    summarise(sumForest = sum(forest))
  
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


################################
## Proportion of open land cover 
openLandCov <- extractPropCLC(codesCLC = c("331", "332", "333", "334", "335"))


##############################
## Proportion of agri 21 cover
agri21Cov <- extractPropCLC(codesCLC = c("211", "213"))


##############################
## Proportion of agri 22 cover
agri22Cov <- extractPropCLC(codesCLC = c("221", "222", "223"))


############
## Proportion of agri 23 cover
agri23Cov <- extractPropCLC(codesCLC = "231")


############
## Proportion of agri 24 cover
agri24Cov <- extractPropCLC(codesCLC = c("241", "242", "243", "244"))


###############
## ROAD DATA ##
###############
# Route500 data from 2012, 2016, and 2020
roads2012 <- shapefile("data/route500/2012/TRONCON_ROUTE.shp")
highways2012 <- roads2012[roads2012$VOCATION == "Type autoroutier", ]
roads2016 <- shapefile("data/route500/2016/TRONCON_ROUTE.shp")
highways2016 <- roads2016[roads2016$VOCATION == "Type autoroutier", ]
roads2020 <- shapefile("data/route500/2020/TRONCON_ROUTE.shp")
highways2020 <- roads2020[roads2020$VOCATION == "Type autoroutier", ]


#######################
## Distance to highways
highways2012Tr <- spTransform(highways2012, CRS(proj4string(gridFrComplete)))
highways2016Tr <- spTransform(highways2016, CRS(proj4string(gridFrComplete)))
highways2020Tr <- spTransform(highways2020, CRS(proj4string(gridFrComplete)))
listHighways <- list()
listHighways[[1]] <- highways2012Tr
listHighways[[2]] <- highways2016Tr
listHighways[[3]] <- highways2020Tr

# Cell centroids
gridCentroid <- gCentroid(gridFrComplete, byid = TRUE)
# Distance for each centroid
distHgwCov <- data.frame(ID = gridFrComplete$ID)

for(i in 1:length(listHighways)){
  
  distHgws <- rep(NA, length(gridFrComplete))
  
  for(j in 1:length(gridCentroid)){
    distHgws[j] <- gDistance(gridCentroid[j], listHighways[[i]])
  }
  
  distHgwCov <- cbind.data.frame(distHgwCov, as.data.frame(distHgws))
  if(i == 1){
    colnames(distHgwCov)[i+1] <- "distHgw2012"
  } else if(i == 2){
    colnames(distHgwCov)[i+1] <- "distHgw2016"
  } else if(i == 3){
    colnames(distHgwCov)[i+1] <- "distHgw2020"
  }
}


##############
## Road length 
roads2012Tr <- crop(spTransform(roads2012, CRS(proj4string(gridFrComplete))), extent(gridFrComplete))
roads2016Tr <- crop(spTransform(roads2016, CRS(proj4string(gridFrComplete))), extent(gridFrComplete))
roads2020Tr <- crop(spTransform(roads2020, CRS(proj4string(gridFrComplete))), extent(gridFrComplete))
listRoads <- list()
listRoads[[1]] <- roads2012Tr
listRoads[[2]] <- roads2016Tr
listRoads[[3]] <- roads2020Tr

# Total road length in each cell
roadLengthCov <- data.frame(ID = gridFrComplete$ID)

for(i in 1:length(listRoads)){
  
  roadsIn <- raster::intersect(listRoads[[i]], gridFrComplete) # extract all the roads inside each cell
  lengthAllRoads <- gLength(roadsIn, byid = TRUE) # calculate the length of each type of road per cell ID
  lengthRoadsCells <- aggregate(lengthAllRoads, by = list(ID = roadsIn$FID), FUN = function(x)sum(x, na.rm = TRUE)) # and sum the length per cell ID
  withoutRoads <- setdiff(out$FID, lengthRoadsCells$ID) # cells without roads in it
  lengthRoadsID <- c(as.numeric(as.character(lengthRoadsCells$ID)), as.numeric(as.character(withoutRoads)))
  lengthRoads <- c(lengthRoadsCells$x, rep(0, length(withoutRoads))) # put 0 as the length of roads in the empty cells
  lengthRoadsOrder <- cbind(lengthRoadsID, lengthRoads)
  lengthRoadsOrderID <- lengthRoadsOrder[order(match(lengthRoadsOrder[,"lengthRoadsID"] ,out@data$FID)),]

  roadLengthCov <- cbind.data.frame(roadLengthCov, as.data.frame(lengthRoadsOrderID))
  if(i == 1){
    colnames(roadLengthCov)[i+1] <- "roadLength2012"
  } else if(i == 2){
    colnames(roadLengthCov)[i+1] <- "roadLength2016"
  } else if(i == 3){
    colnames(roadLengthCov)[i+1] <- "roadLength2020"
  }
}


###################
## Human density ##
###################



##############
## Rugosity ##
##############

#######################
## Prey availability ##
#######################

#######################
## Human disturbance ##
#######################




save(forestCov, connectForCov, shrubCov, openLandCov, agri21Cov, agri22Cov, 
     agri23Cov, distHgwCov, roadLengthCov, file = "covariates.RData")