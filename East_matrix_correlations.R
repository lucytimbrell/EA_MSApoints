# Project: Eastern African Middle Stone Age point variability
# Process: Multiple matrix regressions
# Author: L. J. Timbrell
# Date: 27/03/2023

library("raster")
library("rgdal")
library("factoextra")
library("NbClust")
library("dendextend")
library("ggplot2")
library("plyr")
library("gdistance")

setwd("/YOUR/WORKING/DIRECTORY")

rm(list = ls())

#### Set extent #### 
points <- c(30,55,-9,20) # eastern Africa

#### Import and process climate data ####
bio12_rasterbrick <- raster::brick("bio12_800ka.nc",varname="bio12") # download from https://www.nature.com/articles/s41597-021-01009-3
af_precip <- crop(bio12_rasterbrick, points) #crop to eastern africa
af_precip <- stack(af_precip) #reorder
names(af_precip@layers) <- 799:0 #reorder
af_precip@layers <- af_precip@layers[order(as.numeric(names(af_precip@layers)))] #reorder from present to past

bio01_rasterbrick <- raster::brick("bio01_800ka.nc",varname="bio01")  # download from https://www.nature.com/articles/s41597-021-01009-3
af_temp <- crop(bio01_rasterbrick, points) #crop to eastern africa
af_temp <- stack(af_temp) #reorder
names(af_temp@layers) <- 799:0 #reorder
af_temp@layers <- af_temp@layers[order(as.numeric(names(af_temp@layers)))] #reorder from present to past


#### Import and process site data ####
sites <- read.csv("final_meanshape_pcs_Dec2023.csv") #load data
sites <- subset(sites, !is.na(N)) #remove any sites without lat lon 
sites_sp <- SpatialPointsDataFrame(sites[,c("E", "N")], sites, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") ) #turn to spatial object
sites_sp<- spTransform(sites_sp, CRSobj = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") #projection
mid_ages <- sites$Mid.age #extract mid ages of sites
mid_ages[is.na(mid_ages)] <- 1 # give undated artefacts a date of 1 to maintain order and then change back to NA 


#### Import and process bathymetry and sea level  data ####
bathymetry <- raster("bathymetry_model.tif") #load data
bathymetry_2 <- crop(bathymetry, points) #crop to eastern Africa
bathymetry_3 <- projectRaster(bathymetry_2, 
                              crs = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ") #projection

sea_level <- read.delim("global_sea_level_0_800ka_data.txt", 
                        stringsAsFactor = FALSE)  #load data
sea_level <- sea_level[-1, ] #remove first row (present sea-Level)

#### Process temperature mid-age (occupational phases) layers ####
a <- stack()
for(i in 1:length(mid_ages)){ 
  date <- as.numeric(mid_ages[i]) 
  layer  <- af_temp[[date+1]] #extract layer from stack (adding 1 as climate stack starts at 0 not 1)
  a <- stack(a, layer) 
}
names(a) <- paste0('Temperature_', sites$Edited.photo.name) #assign layer names based on what site they represent
crs(a) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
a = projectRaster(a, crs="+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

## Increase resolution and mask sea level 
a1 <- stack()
for(i in 1:nlayers(a)){
  resample_layer <- disaggregate(a[[i]], fact = 12, method = "bilinear") #disaggregate cells to increase resolution
  bathymetry_4 <- resample(bathymetry_3, resample_layer) #match resolution to climate model
  bathymetry_5 <- bathymetry_4; bathymetry_5[bathymetry_5 <= sea_level[mid_ages[i],2]] <- NA #create mask according to sea level
  resample_layer <- mask(resample_layer, bathymetry_5) #mask
  a1 <- stack(a1, resample_layer)
}

#### Process precipitation mid-age (occupational phases) layers####
b <- stack()
for (i in 1:length(mid_ages)){
  date <- as.numeric(mid_ages[i]) 
  layer <- af_precip[[date+1]] #extract layer from stack (adding 1 as climate stack starts at 0 not 1)
  b <- stack(b, layer) 
}
names(b) <- paste0('Precipitation_', sites$X) # assign layer names based on what site they represent
crs(b) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
b = projectRaster(b, crs="+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

## Increase resolution and mask sea level 
b1 <- stack()
for(i in 1:nlayers(b)){
  resample_layer <- disaggregate(b[[i]], fact = 12, method = "bilinear") #disaggregate cells to increase resolution
  bathymetry_4 <- resample(bathymetry_3, resample_layer) #match resolution to climate model
  bathymetry_5 <- bathymetry_4; bathymetry_5[bathymetry_5 <= sea_level[mid_ages[i],2]] <- NA #create mask according to sea level
  resample_layer <- mask(resample_layer, bathymetry_5) #mask
  b1 <- stack(b1, resample_layer)
}

#### Extract mean and standard deviation for each site ####
mean_p1 <- list()

for(i in 1:length(mid_ages)){
  mean_p1[i] <- raster::extract(b1[[i]], t(sites_sp@coords[i,]), buffer=50000)}
mean_p1 <- lapply(mean_p1, na.omit)
mean_p1_mean <- lapply(mean_p1, mean)
mean_p1_sd <- lapply(mean_p1, sd)
mean_p1_mean_50km <- dist(mean_p1_mean)
names(mean_p1_mean_50km) <- sites$AssID

mean_t1 <- list()

for(i in 1:length(mid_ages)){
  mean_t1[i] <- raster::extract(a1[[i]], t(sites_sp@coords[i,]), buffer=50000)}
mean_t1 <- lapply(mean_t1, na.omit)
mean_t1_mean <- lapply(mean_t1, mean)
mean_t1_sd <- lapply(mean_t1, sd)
mean_t1_mean_50km <- dist(mean_t1_mean)
names(mean_t1_mean_50km) <- sites$AssID

#### Create table of climate values ####
mean_t1_matrix <- as.matrix(unlist(mean_t1_mean))
mean_p1_matrix <- as.matrix(unlist(mean_p1_mean))
sd_t1_matrix <- as.matrix(unlist(mean_t1_sd))
sd_p1_matrix <- as.matrix(unlist(mean_p1_sd))
mean_matrix <- cbind(mean_t1_matrix, sd_t1_matrix, mean_p1_matrix, sd_p1_matrix)
rownames(mean_matrix) <- sites$X
colnames(mean_matrix) <- c("Mean_temperature", "SD_temperature", "Mean_precipiation", "SD_precipitation")


final_dataset <- cbind(sites, mean_matrix)
write.csv(final_dataset, "final_dataset_mmr.csv")

# Multiple Matrix Regressions - RUN FROM HERE FROM NOW ON 

final_dataset <- read.csv("final_dataset_mmr.csv")

#### Set extent #### 
points <- c(30,55,-9,20) # eastern Africa

#### Import and process bathymetry  data ####

#### Import and process site data ####
sites <- read.csv("final_meanshape_pcs.csv") #load data
sites <- subset(sites, !is.na(N)) #remove any sites without lat lon 
sites_sp <- SpatialPointsDataFrame(sites[,c("E", "N")], sites, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") ) #turn to spatial object
sites_sp<- spTransform(sites_sp, CRSobj = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") #projection

bathymetry <- raster("bathymetry_model.tif") #load data
bathymetry_2 <- crop(bathymetry, points) #crop to eastern Africa
bathymetry_3 <- projectRaster(bathymetry_2, 
                              crs = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ") #projection

# Costdistance
altDiff <- function(x){x[2] - x[1]} # Function to calculate altitude difference between cells
hd <- transition(bathymetry_3, altDiff, 8, symm = TRUE) #  Create transition layer
slope <- geoCorrection(hd, type = "c") # Geocorrection
adj <- adjacencyFromTransition(slope) # Identify cells that are adjacent
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05)) # Calculate speed of travelling between adjacent cells
conductance <- geoCorrection(speed, type = "c") # Conductance transition layer

cd <- costDistance(conductance, sites_sp) 
cd <- as.matrix(cd)
rownames(cd) <- sites$Assemblage_2
colnames(cd) <- sites$Assemblage_2
cd <- dist(cd)

#### Distance matrices ####

Length <- dist(final_dataset$Length)
Width <- dist(final_dataset$Width)
Thickness <- dist(final_dataset$Thickness)

PC1 <- dist(final_dataset$PC1)
PC2 <- dist(final_dataset$PC2)
PC3 <- dist(final_dataset$PC3)
PC4 <- dist(final_dataset$PC4)
PC5 <- dist(final_dataset$PC5)
PC6 <- dist(final_dataset$PC6)


Rawmat <- dist(final_dataset[,18:23], method = "binary")

Time <- dist(final_dataset$Mid.age)
Prec <- dist(final_dataset$Mean_precipiation)
Temp <- dist(final_dataset$Mean_temperature)

#### scale dists ####

distances <- list(length = Length, width = Width, thickness = Thickness, pc1 = PC1, pc2 = PC2, pc3 = PC3, pc4 = PC4,pc5 = PC5, pc6 = PC6, rawmat = Rawmat, time = Time, space = cd, prec = Prec, temp = Temp) 

scale_fun <- function(x){(x-min(x))/(max(x)-min(x))} # Function to scale

scaled_distances <- lapply(distances, scale_fun)
scaled_distances_matrix <- lapply(scaled_distances, as.matrix)

names <- as.character(sites$Assemblage_2)
distances2 <- lapply(scaled_distances_matrix, "rownames<-", names)
distances2 <- lapply(distances2, "colnames<-", names)
distances3 <- lapply(distances2, function(x) as.dist(x))

### Multiple Matrix Regression - PC1 #### 


  behav <- distances3[[4]]
  distances4 <- distances3
  distances4 <- distances4[-(4:9)]
  PC1_multi <- phytools::multi.mantel(behav, distances4, nperm=999)

  ### Multiple Matrix Regression - PC2 #### 
  behav <- distances3[[5]]
  distances4 <- distances3
  distances4 <- distances4[-(4:9)]
  PC2_multi <- phytools::multi.mantel(behav, distances4, nperm=999)

  ### Multiple Matrix Regression - PC3 #### 
  behav <- distances3[[6]]
  distances4 <- distances3
  distances4 <- distances4[-(4:9)]
  PC3_multi <- phytools::multi.mantel(behav, distances4, nperm=999)

  ### Multiple Matrix Regression - PC4 #### 
  behav <- distances3[[7]]
  distances4 <- distances3
  distances4 <- distances4[-(4:9)]
  PC4_multi <- phytools::multi.mantel(behav, distances4, nperm=999)


  ### Multiple Matrix Regression - PC5 #### 
  behav <- distances3[[8]]
  distances4 <- distances3
  distances4 <- distances4[-(4:9)]
  PC5_multi<- phytools::multi.mantel(behav, distances4, nperm=999)
  
  ### Multiple Matrix Regression - PC6 #### 
  behav <- distances3[[9]]
  distances4 <- distances3
  distances4 <- distances4[-(4:9)]
  PC6_multi <- phytools::multi.mantel(behav, distances4, nperm=999)

PC1_multi
PC2_multi
PC3_multi
PC4_multi
PC5_multi
PC6_multi




### Multiple Matrix Regression - Length #### 
behav <- distances3[[1]]
distances4 <- distances3
distances4 <- distances4[-(4:9)]
distances4 <- distances4[-(1)]
length_multi <- phytools::multi.mantel(behav, distances4, nperm=999)

### Multiple Matrix Regression - Width #### 
behav <- distances3[[2]]
distances4 <- distances3
distances4 <- distances4[-(4:9)]
distances4 <- distances4[-(2)]
width_multi <- phytools::multi.mantel(behav, distances4, nperm=999)

### Multiple Matrix Regression - Thickness #### 
behav <- distances3[[3]]
distances4 <- distances3
distances4 <- distances4[-(4:9)]
distances4 <- distances4[-(3)]
thickness_multi <- phytools::multi.mantel(behav, distances4, nperm=999)





library(ggpubr)
par(mfrow=c(2,3))
a <- ggscatter(final_dataset, x = "PC1", y = "N", 
          add = "reg.line", conf.int = TRUE,
          cor.method = "pearson",
          xlab = "PC1", ylab = "Latitude")

b <- ggscatter(final_dataset, x = "PC1", y = "Mid.age", 
          add = "reg.line", conf.int = TRUE,
          cor.method = "pearson",
          xlab = "PC1", ylab = "Chronology")

c <- ggscatter(final_dataset, x = "PC2", y = "Mid.age", 
          add = "reg.line", conf.int = TRUE,
          cor.method = "pearson",
          xlab = "PC2", ylab = "Chronology")

d <- ggscatter(final_dataset, x = "PC2", y = "Mean_precipiation", 
          add = "reg.line", conf.int = TRUE,
          cor.method = "pearson",
          xlab = "PC2", ylab = "Precipitation")

e <- ggscatter(final_dataset, x = "PC6", y = "Mean_temperature", 
          add = "reg.line", conf.int = TRUE,
          cor.method = "pearson",
          xlab = "PC6", ylab = "Temperature")

(a|b|c)/(d|e) +
  plot_annotation(tag_levels = "a",
                  theme = theme(plot.title = element_text(size = 30))) 


ggscatter(final_dataset, x = "PC2", y = "Mean_precipiation", 
                add = "reg.line", conf.int = TRUE,
                cor.method = "pearson",
                xlab = "PC1", ylab = "Precipitation")

