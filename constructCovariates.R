setwd("C:/Users/sarah.bauduin/Documents/GitHub/appendix_lynxOccupancy/")
# Appendix Bauduin, ... and Gimenez

library(raster)
library(sf)
library(dplyr)
library(rgeos)
library(spatialEco)
library(geojsonio)
library(stringr)
library(tidyr)
library(janitor)

# Covariates to run the lynx occupancy model
# Grid cell on which covariates must be computed
load("data/gridFrComplete.RData")
# Crop the grid to the study area extent (east of France)
gridFrComplete <- crop(gridFrComplete, extent(c(720000, 1090000, 6260000, 6920000)))


##########################
## LANDSCAPE COVER DATA ##
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
colnames(forestCov) <- c("ID", "forest1990", "forest2000", "forest2006", "forest2012", "forest2018")


######################
## Forest connectivity 
# The number of "forested cells" (i.e., cells with more than 50% of forest cover),
# in a 10 km buffer
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
  
  # 10 km buffer around all cells
  gridForest3 <- gridForest2 %>% 
    st_buffer(dist = 9000) %>% 
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
     agri23Cov, agri24Cov, file = "outputs/covariatesCover.RData")

##########################################################


###############
## ROAD DATA ##
###############
# Route500 data from 2012, 2015, and 2018
roads2012 <- shapefile("data/route500/2012/TRONCON_ROUTE.shp")
highways2012 <- roads2012[roads2012$VOCATION == "Type autoroutier", ]
roads2015 <- shapefile("data/route500/2015/TRONCON_ROUTE.shp")
highways2015 <- roads2015[roads2015$VOCATION == "Type autoroutier", ]
roads2018 <- shapefile("data/route500/2018/TRONCON_ROUTE.shp")
highways2018 <- roads2018[roads2018$VOCATION == "Type autoroutier", ]


#######################
## Distance to highways
highways2012Tr <- spTransform(highways2012, CRS(proj4string(gridFrComplete)))
highways2015Tr <- spTransform(highways2015, CRS(proj4string(gridFrComplete)))
highways2018Tr <- spTransform(highways2018, CRS(proj4string(gridFrComplete)))
listHighways <- list()
listHighways[[1]] <- highways2012Tr
listHighways[[2]] <- highways2015Tr
listHighways[[3]] <- highways2018Tr

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
    colnames(distHgwCov)[i+1] <- "distHgw2015"
  } else if(i == 3){
    colnames(distHgwCov)[i+1] <- "distHgw2018"
  }
}


##############
## Road length 
roads2012Tr <- crop(spTransform(roads2012, CRS(proj4string(gridFrComplete))), extent(gridFrComplete))
# Doesn't work with sp. Need to use sf
#roads2015Tr <- crop(spTransform(roads2015, CRS(proj4string(gridFrComplete))), extent(gridFrComplete))
roads2015TrSF <- st_crop(st_as_sf(spTransform(roads2015, CRS(proj4string(gridFrComplete)))), st_as_sf(gridFrComplete))
roads2015Tr <- sf::as_Spatial(roads2015TrSF)
roads2018Tr <- crop(spTransform(roads2018, CRS(proj4string(gridFrComplete))), extent(gridFrComplete))
listRoads <- list()
listRoads[[1]] <- roads2012Tr
listRoads[[2]] <- roads2015Tr
listRoads[[3]] <- roads2018Tr

# Total road length in each cell
roadLengthCov <- data.frame(ID = gridFrComplete$ID)

# Need to cut the grid into cell sequences otherwise the following functions fail with too many cells
cutSeq <- seq(1, 2131, 100) 
cutCells <- list()
for(i in 1:(length(cutSeq) - 1)){
  cutCells[[i]] <- seq(cutSeq[i], cutSeq[i + 1] - 1, by = 1)
}
cutCells[[22]] <- 2101:2131

for(i in 1:length(listRoads)){
  
  lengthRoadsCells <- data.frame(ID = numeric(), x = numeric())
  for(j in 1:length(cutCells)){
    
    roadsIn <- raster::intersect(listRoads[[i]], gridFrComplete[cutCells[[j]],]) # extract all the roads inside each cell
    lengthAllRoads <- gLength(roadsIn, byid = TRUE) # calculate the length of each type of road per cell ID
    lengthRoadsCellCut <- aggregate(lengthAllRoads, by = list(ID = roadsIn$ID), FUN = function(x)sum(x, na.rm = TRUE)) # and sum the length per cell ID
    lengthRoadsCells <- rbind.data.frame(lengthRoadsCells, lengthRoadsCellCut)
    save(lengthRoadsCells, file = "outputs/lengthRoadsCellsTemp.RData")
    print(j)
  }
  
  withoutRoads <- setdiff(gridFrComplete$ID, lengthRoadsCells$ID) # cells without roads in it
  lengthRoadsID <- c(as.numeric(as.character(lengthRoadsCells$ID)), as.numeric(as.character(withoutRoads)))
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
}


##################
## Save covariates
save(distHgwCov, roadLengthCov, file = "outputs/covariatesRoads.RData")

###############################################################


###########################
## HUMAN POPULATION DATA ##
###########################
# Human density for cells of 1km2 in 2017
human <- shapefile("data/humanPop/2017/Filosofi2017_carreaux_1km_met.shp")
human_sf <- st_as_sf(human)

gridFr_sf <- st_as_sf(gridFrComplete) %>%
  st_transform(crs = st_crs(human_sf))


################
## Human density
# Compute the mean human density (1km2) into the 10km2 cells
humanDens <- st_intersection(gridFr_sf , human_sf) %>%
  group_by(ID) %>%
  summarise(humanDens1km2 = mean(Ind, na.rm = TRUE)) %>%
  as_tibble() %>%
  select(ID, humanDens1km2)
# Fill with 0 the missing data
hDensCov <- data.frame(ID = gridFrComplete$ID)
hDensCov <- merge(hDensCov, as.data.frame(humanDens), all = TRUE)
hDensCov[is.na(hDensCov)] <- 0


####################
## ELEVATION DATA ##
####################
# Elevation data on 25 m x 25 m resolution
elev1 <- raster("data/elevation/032ab314564b9cb72c98fbeb093aeaf69720fbfd/eu_dem_v11_E30N20.TIF")# 
gridFrCompleteTr <- spTransform(gridFrComplete, elev1@crs)
elev1 <- crop(elev1, extent(gridFrCompleteTr))
elev2 <- raster("data/elevation/97824c12f357f50638d665b5a58707cd82857d57/eu_dem_v11_E40N20.TIF")
elev2 <- crop(elev2, extent(gridFrCompleteTr))
elevation <- merge(elev1, elev2)


#############
## Ruggedness
# Compute TRI (terrain ruggedness index) on 1 km2 cells
elevAggr <- aggregate(elevation, fact = 40)
elevTRI <- tri(elevAggr)

# Compute the number of 1km2 "rugged cells" into 10 km2 cells
elevTRI[elevTRI < 162] <- 0 # 162 = "intermediately rugged surface"
elevTRI[elevTRI >= 162] <- 0.01
elevTRI_10km2 <- aggregate(elevTRI, fact = 10, fun = sum)
triCells <- raster::extract(elevTRI_10km2, gCentroid(gridFrCompleteTr, byid = TRUE))
triCov <- cbind.data.frame(ID = gridFrCompleteTr@data$ID, tri = triCells)


################
## TALLY DATA ##
################
# Deer number allowed for hunting
raw_attrib <- read.csv("data/tally/attributionsChevreuil_v2.txt", sep = ";")
# Real number of deer hunted
raw_realiz <- read.csv("data/tally/realisationsChevreuil_v2.txt", sep = ";")
# "Communes" polygons
communes <- geojson_read("data/tally/communes.geojson",  what = "sp")


####################
## Prey availability
# Merge deer number with "communes"
communes@data$nom <- str_to_lower(communes@data$nom) # to lowercase
communes@data$nom <- stringi::stri_trans_general(communes@data$nom, "Latin-ASCII") # remove accent

# Clean number of deer allowed and hunted
attrib <- raw_attrib %>%
  janitor::clean_names() %>%
  mutate(nom = str_to_lower(nom)) %>% # to lowercase
  mutate(nom = stringi::stri_trans_general(nom, "Latin-ASCII")) %>% # remove accent
  pivot_longer(cols = x1985:x2020, names_to = "year", values_to = "attrib") %>%
  mutate(year = str_remove(year, "x"),
         year = as.integer(year),
         code_insee = as.character(code_insee))

realiz <- raw_realiz %>%
  janitor::clean_names() %>%
  mutate(nom = str_to_lower(nom)) %>% # to lowercase
  mutate(nom = stringi::stri_trans_general(nom, "Latin-ASCII")) %>% # remove accent
  pivot_longer(cols = x1985:x2020, names_to = "year", values_to = "realiz") %>%
  mutate(year = str_remove(year, "x"),
         year = as.integer(year),
         code_insee = as.character(code_insee))

# Remove data for which there is no commune polygon associated
attrib <- attrib[attrib$nom %in% communes$nom,]
attrib <- attrib[!is.na(attrib$attrib),]
realiz <- realiz[realiz$nom %in% communes$nom,]
realiz <- realiz[!is.na(realiz$realiz),]
# Merge number allowed and hunted
attribReal <- attrib %>% 
  full_join(realiz)
# And with the communes
communes_all <- communes %>%
  st_as_sf() %>%
  rename("code_insee" = "code")  %>%
  inner_join(attribReal)

# Cut commune polygons with the grid cells
gridFr_sf <- st_as_sf(gridFrComplete) %>%
  st_transform(crs = st_crs(communes_all))

# Need to cut the grid into cell sequences otherwise the following function fails with too many cells
cutSeq <- seq(1, 2131, 50) 
cutCells <- list()
for(i in 1:(length(cutSeq) - 1)){
  cutCells[[i]] <- seq(cutSeq[i], cutSeq[i + 1] - 1, by = 1)
}
cutCells[[43]] <- 2101:2131

gridPreys <- st_intersection(gridFr_sf[cutCells[[1]],], communes_all)
for(j in 2:length(cutCells)){
  gridPreysCut <- st_intersection(gridFr_sf[cutCells[[j]],], communes_all)
  gridPreys <- rbind(gridPreys, gridPreysCut)
  save(gridPreys, file = "outputs/gridPreysTemp.RData")
  print(j)
} 
gridPreys$area <- st_area(gridPreys)

# Compute the portion of the number of deer allowed/hunted in each commune portion in the grid cells
gridPreys2 <- gridPreys %>%
  mutate(portAttrib = as.numeric((area / 1e+08) * attrib)) %>% 
  mutate(portReal = as.numeric((area / 1e+08) * realiz)) %>% 
  group_by(ID, year) %>%
  summarise(sumAttrib = sum(portAttrib, na.rm = TRUE), sumReal = sum(portReal, na.rm = TRUE))
# Put back NA to cells which were all NA
allNA <- gridPreys %>% 
  group_by(ID, year) %>%
  summarise(propNAattrib = sum(is.na(attrib))/length(attrib), propNAreal = sum(is.na(realiz))/length(realiz)) %>% 
  st_drop_geometry() %>%
  as_tibble() 
gridPreys2NA <- gridPreys2 %>% 
  left_join(allNA)
gridPreys2NA[gridPreys2NA$propNAattrib == 1, "sumAttrib"] <- NA
gridPreys2NA[gridPreys2NA$propNAreal == 1, "sumReal"] <- NA

gridPreys3 <- gridFr_sf %>%
  left_join(as.data.frame(gridPreys2NA)[,c("ID", "year", "sumAttrib", "sumReal")])
gridPreys3 <- gridPreys3[!is.na(gridPreys3$year),]

# Compute how much of the allowed number of deer are actually hunted
portHunted <- gridPreys3[(gridPreys3$sumAttrib != 0 & gridPreys3$sumReal != 0 & 
                      gridPreys3$sumReal < gridPreys3$sumAttrib), "sumReal"] / 
          gridPreys3[(gridPreys3$sumAttrib != 0 & gridPreys3$sumReal != 0 & 
                        gridPreys3$sumReal < gridPreys3$sumAttrib), "sumAttrib"]
# Take this portion of allowed number of deer to fill the NA in hunted numbers
gridPreys3Fill <- gridPreys3 %>% 
  mutate(sumReal = ifelse(is.na(sumReal), sumAttrib * mean(portHunted[,1], na.rm = TRUE), sumReal))

# Convert each year data as a column and complete the cells to obtain the full grid
gridPivot <- gridPreys3Fill %>% 
  as_tibble() %>% 
  select(ID, year, sumReal) %>% 
  pivot_wider(names_from = year, values_from = sumReal)
# Remove the year 2015 et 2016, not enough data for the year
gridPivot <- gridPivot %>% 
  select(-"2015", -"2016")

# Merge to the grid cells
gridPreys4 <- gridFr_sf %>% 
  full_join(gridPivot)

# Fill the NA cells with the mean of all 8 surrounding cells
yearList <- c(1985, 1993, 1998, 2002, 2007, 2012, 2017, 2018, 2019, 2020)

# Identify the 8 neighboring cells
buffer <- gridPreys4 %>% 
  st_buffer(dist = 9000) %>% # one cell (a bit less)
  rename(IDbuffer = ID) %>% 
  select(IDbuffer)

bufferCells <- st_intersection(buffer[cutCells[[1]],], gridPreys4)
for(j in 2:length(cutCells)){
  bufferCellsCut <- st_intersection(buffer[cutCells[[j]],], gridPreys4)
  bufferCells <- rbind(bufferCells, bufferCellsCut)
  save(bufferCells, file = "outputs/bufferCellsTemp.RData")
  print(j)
} 

preyCov <- gridPreys4 %>% 
  st_drop_geometry() %>% 
  as.data.frame() %>% 
  select(-area)

# As long as there are NA, fill the cells
for(y in yearList){
  
  while(sum(is.na(preyCov[, as.character(y)])) != 0){
      
    # Cells which are NA
    cellYearNAID <- preyCov %>% 
      select(ID, as.character(y)) %>% 
      filter(is.na(.[, 2])) %>% 
      pull(ID)
    
    # Buffer of cells which are NA
    bufferYearNA <- bufferCells %>% 
      select(IDbuffer, paste0("X", as.character(y))) %>%
      filter(IDbuffer %in% cellYearNAID) %>% 
      st_drop_geometry() %>% 
      group_by(IDbuffer) %>% 
      summarise(!!paste0("mean", quo_name(as.character(y))) := 
               mean(!!as.name(quo_name(paste0("X", as.character(y)))), na.rm = TRUE))

    # Replace the mean value in the NA cells
    preyCov[preyCov$ID %in% bufferYearNA$IDbuffer, as.character(y)] <- pull(bufferYearNA, 2)
    mergeBuff <- merge(preyCov[, c("ID", as.character(y))], data.frame(ID = pull(bufferCells, "ID")), all.y = TRUE)
    newValue <- mergeBuff[match(pull(bufferCells, "ID"), mergeBuff$ID), as.character(y)]
    bufferCells[, paste0("X", as.character(y))] <- newValue
  }
  print(y)
}

colnames(preyCov)[-1] <- paste0("prey", colnames(preyCov)[-1])


###########################
## HUMAN FOOTPRINT INDEX ##
###########################
# HFI data
hfi_1993 <- raster("data/HFI/1993/wildareas-v3-1993-human-footprint.tif")
hfi_2009 <- raster("data/HFI/2009/wildareas-v3-2009-human-footprint.tif")
# 50 is the maximum value on land, then 128 is given for water areas
hfi_1993[hfi_1993 == 128] <- 25 # gives a medium value
hfi_2009[hfi_2009 == 128] <- 25


####################
## Human disturbance 
gridFr_sf <- st_as_sf(gridFrComplete) %>%
  st_transform(crs = st_crs(hfi_1993@crs))
# Mean HFI per grid cells
hfiCells_1993 <- raster::extract(hfi_1993, gridFr_sf)
hfiCells_2009 <- raster::extract(hfi_2009, gridFr_sf)

hfiCov <- cbind.data.frame(ID = gridFrComplete@data$ID, 
                           hfi_1993 = unlist(lapply(hfiCells_1993, function(x) mean(x, na.rm = TRUE))),
                           hfi_2009 = unlist(lapply(hfiCells_2009, function(x) mean(x, na.rm = TRUE))))


##################
## Save covariates
save(hDensCov, triCov, preyCov, hfiCov, file = "outputs/covariatesOthers.RData")

########################################################################


