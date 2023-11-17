library(rgdal)
library(rgeos)
library(secr)

###############################
#### Make secr Habitat Mask ###
###############################

#Read in Capture History
Tiogacapt.comb<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/Tiogacapt.comb_300MSegments.rds")

#Read in Tioga boundary and transform to UTMs.
TiogaHabitat<- readOGR("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/TiogaHabitatBoundary.shp")
TiogaHabitat<-spTransform(TiogaHabitat, CRS("+proj=utm +zone=10 +datum=NAD83")) #transform to UTMs

#scale-down habitat boundary from meters to kilometers
extractCoords <- function(sp.df)
{
  results <- list()
  for(i in 1:length(sp.df@polygons[[1]]@Polygons))
  {
    results[[i]] <- sp.df@polygons[[1]]@Polygons[[i]]@coords
  }
  results
}

vertices<-extractCoords(TiogaHabitat)

meters.to.km<-function(meters){
  km<-meters/1000
}

scaled.vertices<-lapply(vertices, meters.to.km)

Polys<-list()
for(i in 1:length(scaled.vertices)){
  Polys[i]<-sp::Polygon(scaled.vertices[[i]])
}
Polys.plural<-sp::Polygons(Polys, ID = "0")
Polys.sp<-sp::SpatialPolygons(list(Polys.plural), proj4string = CRS("+proj=utm +zone=10 +datum=NAD83"))
Tioga.spdf<-sp::SpatialPolygonsDataFrame(Polys.sp, data = TiogaHabitat@data)

#Playing with suggested buffer size
suggest.buffer(Tiogacapt.comb.telem, detectfn = 'HZ', RBtarget = 0.001)

##Creating a blank mask with boundary of Tioga and buffer of 8,200 metres. 
#Spacing is relative to sigma - resolution should be less than half the expected estimate of sigma (Sutherland et al. 2019).
#(1.108307+1.0083449)/2=1.0958 1.0958 km/2 = 0.547 km. So we will use 0.5 km. 
#Then plotting traps (transects)
TiogaMask.comb.km<-make.mask(traps(Tiogacapt.comb), type = "trapbuffer", buffer=8, spacing = 0.5, poly=Tioga.spdf)
plot(TiogaMask.comb.km)
plot(traps(Tiogacapt.comb), detpar=list(pch=16, cex=0.8), add=TRUE)

##Extracting points from mask and convert them to meters (UTMs) to extract covariates from rasters
TiogaMaskcoords.2018<-as.data.frame(TiogaMask.comb.km$Tioga18_300Msegments)
TiogaMaskcoords.2019<-as.data.frame(TiogaMask.comb.km$Tioga19_300Msegments)
TiogaMaskcoords.2018<-TiogaMaskcoords.2018[,c(1,2)]*1000
TiogaMaskcoords.2019<-TiogaMaskcoords.2019[,c(1,2)]*1000
TiogaMaskcoords.2018<- SpatialPoints(coords=TiogaMaskcoords.2018, proj4string = CRS("+proj=utm +zone=10"))
TiogaMaskcoords.2019<- SpatialPoints(coords=TiogaMaskcoords.2019, proj4string = CRS("+proj=utm +zone=10"))

###Create data frame for extracted values
TiogaMaskextracts2018<-TiogaMaskcoords.2018
TiogaMaskextracts2019<-TiogaMaskcoords.2019

##Read in covariate layers and extract values to created data frame.
###DDE
TiogaDDE<-raster("TiogaDDE_Continuous.tif")
TiogaDDE<-projectRaster(TiogaDDE, crs = "+proj=utm +zone=10 +datum=NAD83") #project to UTMs
DDE.2018<-extract(TiogaDDE, TiogaMaskcoords.2018, buffer = 1000, fun=mean)
DDE.2019<-extract(TiogaDDE, TiogaMaskcoords.2019, buffer = 1000, fun=mean)
summary(DDE.2018)
summary(DDE.2019)

TiogaMaskextracts2018$DDE<-scale(DDE.2018)
TiogaMaskextracts2019$DDE<-scale(DDE.2019)

###Percent Slope
TiogaSlope<-raster("TIOGA_mean_slp1.tif")
Slope.2018<-extract(TiogaSlope, TiogaMaskcoords.2018, method = 'simple')
Slope.2019<-extract(TiogaSlope, TiogaMaskcoords.2019, method = 'simple')
summary(Slope.2018)
summary(Slope.2019)

TiogaMaskextracts2018$PercentSlope<-scale(Slope.2018)
TiogaMaskextracts2019$PercentSlope<-scale(Slope.2019)

###DFE
TiogaDFE<-raster("TIOGADFEint1.tif")
TiogaDFE<-projectRaster(TiogaDFE, crs = "+proj=utm +zone=10 +datum=NAD83") #project to UTMs
DFE.2018<-extract(TiogaDFE, TiogaMaskcoords.2018, buffer = 1000, fun = mean)
DFE.2019<-extract(TiogaDFE, TiogaMaskcoords.2019, buffer = 1000, fun = mean)
summary(DFE.2018)
summary(DFE.2019)

TiogaMaskextracts2018$DFE<-scale(DFE.2018)
TiogaMaskextracts2019$DFE<-scale(DFE.2019)

###Distance to Roads
#TiogaRoads<- readOGR("C:/Users/nelsonj7/Box/secr/TiogaRoads.shp")
#TiogaRoads<-spTransform(TiogaRoads, CRS("+proj=utm +zone=10 +datum=NAD83")) #transform to UTMs
#Template<-TiogaDFE
#Template[]<-NA
#TiogaRoadsRaster<-rasterize(TiogaRoads, Template, field=1)
#saveRDS(TiogaRoadsRaster, "TiogaRoadsRaster.rds")
#TiogaRoads.dist<-distance(TiogaRoadsRaster)
#saveRDS(TiogaRoads.dist, file="TiogaRoadsDist.rds")
TiogaRoads.dist<-readRDS("C:/Users/nelsonj7/Box/secr/TiogaRoadsDist.rds")
#writeRaster(TiogaRoads.dist, file = "C:/Users/nelsonj7/Box/secr/TiogaRoadsDist.tif")
TiogaRoads.dist.scale<-scale(TiogaRoads.dist)
TiogaMaskextracts2018$Roads<-extract(TiogaRoads.dist.scale, TiogaMaskcoords.2018, method = 'simple')
TiogaMaskextracts2019$Roads<-extract(TiogaRoads.dist.scale, TiogaMaskcoords.2019, method = 'simple')

###Precipitation
TiogaPrecip.2019<-raster("C:/Users/nelsonj7/Box/secr/PRISM/TiogaPrecip_2019.sampling_avg.tif")
TiogaPrecip.2018<-raster("C:/Users/nelsonj7/Box/secr/PRISM/TiogaPrecip_2018.sampling_avg.tif")

TiogaPrecip.2018.extract<-extract(TiogaPrecip.2018, TiogaMaskcoords.2018, method = 'simple')
TiogaPrecip.2019.extract<-extract(TiogaPrecip.2019, TiogaMaskcoords.2019, method = 'simple')
summary(TiogaPrecip.2018.extract)
summary(TiogaPrecip.2019.extract)

TiogaMaskextracts2018$Precip<-scale(TiogaPrecip.2018.extract)
TiogaMaskextracts2019$Precip<-scale(TiogaPrecip.2019.extract)

###All covariates added
head(TiogaMaskextracts2018)
head(TiogaMaskextracts2019)

saveRDS(TiogaMaskextracts2018, file = "TiogaMaskextracts.2018.comb_8000buffer_500spacing.rds")
saveRDS(TiogaMaskextracts2019, file = "TiogaMaskextracts.2019.comb_8000buffer_500spacing.rds")

TiogaMaskextracts2018<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/TiogaMaskextracts.2018.comb_8000buffer_500spacing.rds")
TiogaMaskextracts2019<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/TiogaMaskextracts.2019.comb_8000buffer_500spacing.rds")

###Add extracted values as covariates to mask
covariates(TiogaMask.comb.km$Tioga18_300Msegments)<-TiogaMaskextracts2018
covariates(TiogaMask.comb.km$Tioga19_300Msegments)<-TiogaMaskextracts2019
saveRDS(TiogaMask.comb.km, file = "TiogaMask.comb.km_8000buffer500res.rds")
#TiogaMask.comb.km<-readRDS("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/TiogaMask.comb.km_15000buffer.rds")
#polyarea(TiogaMask.comb.km)


###Spatial Location covariate is within the mask as x and y.