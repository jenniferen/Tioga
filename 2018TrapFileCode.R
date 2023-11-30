library(dplyr)
library(raster)
library(sp)
library(sf)

setwd("~/OSU/ElkSCR/")
data<-read.csv("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2018TiogaTransectTracks.csv")
head(data)

#Look at distances among all detections
#Read in coordinates and sex of identified elk
TiogaElk2018<-read.csv("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2018TiogaElkCoordinatesUTMSex.csv", header=TRUE)
head(TiogaElk2018)
colnames(TiogaElk2018)<-c("X", "TransectID", "LabID", "Elk.Assignment", "UTM_E", "UTM_N", "Lon", "Lat")

#Read in coordinates and sex of identified elk
TiogaElk2019<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2019TiogaElkCoordinatesUTMSex.txt", header=TRUE)
head(TiogaElk2019)

#Combine coordinates of Tioga2018 and 2019 into one dataframe and compute distance matrix among all detections between both years
coords.combined<-rbind(TiogaElk2018[,c(5,6)], TiogaElk2019[,c(1,2)])

dst<-dist(coords.combined)
dst<-data.matrix(dst)
dim<-ncol(dst)
dst<-dst[!dst==0]
hist(dst, breaks = 7637, xlim = c(0,1500), main = "Distances among all detected elk \n in 2018 and 2019 - Tioga",
     xlab = "Distance (meters)")
range(dst)

###What do I want to do? Group each transect, calculate distances between points,
##choose distance to separate by and extract those points. Choose segment distance 
#using histograms above. Going with 400 meters.
data.sp<-SpatialPointsDataFrame(coords = data[,8:9], 
                                   data = data,
                                   proj4string = CRS(as.character("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))

data.sp.utm<- spTransform(data.sp, CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"))

data.sp.utm$UTM<-coordinates(data.sp.utm)
data.sp.utm$X<-data.sp.utm$UTM[,1]
data.sp.utm$Y<-data.sp.utm$UTM[,2]
data.sp.utm$dist<-NA

#For each group of transects, Calculate distance between consecutive points, place 0 if first point.
for(i in 1:nrow(data.sp.utm)) {
  if(i > 1 && data.sp.utm$Transect[i] == data.sp.utm$Transect[i-1]) {
    data.sp.utm$dist[i] <- sqrt(((data.sp.utm$X[i] - data.sp.utm$X[i-1]) ^ 2) + ((data.sp.utm$Y[i] - data.sp.utm$Y[i-1]) ^ 2))
  } else {
    data.sp.utm$dist[i] <- NA
  }
}

data.sp.utm$dist.sum<-NA
for(i in 1:nrow(data.sp.utm)) {
  if(i > 1 && data.sp.utm$Transect[i] == data.sp.utm$Transect[i-1]) {
      data.sp.utm$dist.sum[i] <- data.sp.utm$dist.sum[i-1] + data.sp.utm$dist[i]
  } else {
    data.sp.utm$dist.sum[i] <- 0
  }
}

data.sp.utm$dist.sum
max(data.sp.utm$dist.sum)

#Now, I want to extract points that occur in increments of x-specified meters

seg.200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 200))), ]))
seg.300<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 300))), ]))
seg.400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 400))), ]))
seg.600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 600))), ]))
seg.800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 800))), ]))
seg.900<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 900))), ]))
seg.1000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1000))), ]))
seg.1200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1200))), ]))
seg.1400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1400))), ]))
seg.1500<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1500))), ]))
seg.1600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1600))), ]))
seg.1800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 1800))), ]))
seg.2000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2000))), ]))
seg.2100<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2100))), ]))
seg.2200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2200))), ]))
seg.2400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2400))), ]))
seg.2600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2600))), ]))
seg.2700<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2700))), ]))
seg.2800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 2800))), ]))
seg.3000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3000))), ]))
seg.3200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3200))), ]))
seg.3300<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3200))), ]))
seg.3400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3400))), ]))
seg.3600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3600))), ]))
seg.3800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3800))), ]))
seg.3900<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 3900))), ]))
seg.4000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4000))), ]))
seg.4200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4200))), ]))
seg.4400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4400))), ]))
seg.4500<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4500))), ]))
seg.4600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4600))), ]))
seg.4800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 4800))), ]))
seg.5000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5000))), ]))
seg.5100<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5100))), ]))
seg.5200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5200))), ]))
seg.5400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5400))), ]))
seg.5600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5600))), ]))
seg.5700<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5700))), ]))
seg.5800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 5800))), ]))
seg.6000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 6000))), ]))
seg.6200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 6200))), ]))
seg.6300<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 6300))), ]))
seg.6400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 6400))), ]))
seg.6600<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 6600))), ]))
seg.6800<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 6800))), ]))
seg.6900<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 6900))), ]))
seg.7000<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 7000))), ]))
seg.7200<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 7200))), ]))
seg.7400<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 7400))), ]))
seg.7500<-do.call('rbind', by(data.sp.utm, data.sp.utm$Transect, function(x) x[c(which.min(abs(x$dist.sum - 7500))), ]))
#Create 200 meter segments
#seg.comb<-rbind(seg.200, seg.400, seg.600, seg.800, seg.1000, seg.1200, seg.1400, seg.1600, seg.1800, seg.2000, seg.2200,
#            seg.2400, seg.2600, seg.2800, seg.3000, seg.3200, seg.3400, seg.3600, seg.3800, seg.4000, seg.4200, seg.4400,
#            seg.4600, seg.4800, seg.5000, seg.5200, seg.5400, seg.5600, seg.5800, seg.6000, seg.6200, seg.6400, seg.6600,
#            seg.6800, seg.7000, seg.7200, seg.7400)

#Create 300 meter segments
seg.comb<-rbind(seg.300, seg.600, seg.900, seg.1200, seg.1500, seg.1800, seg.2100, seg.2400,
                seg.2700, seg.3000, seg.3300, seg.3600, seg.3900, seg.4200, seg.4500, seg.4800,
                seg.5100, seg.5400, seg.5700, seg.6000, seg.6300, seg.6600, seg.6900, seg.7200,
                seg.7500)

#Create 400 meter segments
#seg.comb<-rbind(seg.400, seg.800, seg.1200, seg.1600, seg.2000,
#            seg.2400, seg.2800, seg.3200, seg.3600, seg.4000, seg.4400,
#            seg.4800, seg.5200, seg.5600, seg.6000, seg.6400,
#            seg.6800, seg.7200)

seg.comb.trim<-seg.comb[!duplicated(seg.comb$OBJECTID),]




#Create trap file
traps<-seg.comb.trim[,-c(2:12,14,17:22)]
traps<-traps[,c(1,3:4,2)]
head(traps)

#Need TrapID, X, Y, character, effort, Julian, precipitation, and # of observers.
TrapID<-seq(1:nrow(traps)) #Trap ID
traps$TrapID<-TrapID
traps$char<-rep("/", nrow(traps))

#effort covariate: want to take the max dist of each transect, scale and add to each point.
Effort<-NULL
Transects<-as.list(unique(seg.comb.trim$Transect))
for(i in Transects){
  sub<-subset(seg.comb.trim, seg.comb.trim$Transect == i)
  Effort$Transect[i]<- i
  Effort$dist[i]<-max(sub$dist.sum)
}

Effort$dist.scale<-scale(Effort$dist)

traps$effort <- Effort$dist.scale[match(traps$Transect, Effort$Transect)]

#Precipitation

Precip.julian.end<-traps$JulianDay+15
Precip.julian.start<-traps$JulianDay+1
Precip.df<-as.data.frame(cbind(traps$X, traps$Y, Precip.julian.end, Precip.julian.start))
#
coordinates(Precip.df)<-Precip.df[,c(1,2)]
precip.cov<-rep(NA, nrow(Precip.df))
test.precip<-NULL
x1_proj.2018<-raster("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/PRISM/TiogaPrecip_2018_daily.tif")

for(i in 1:nrow(Precip.df)){
  test.precip<-calc(x1_proj.2018[[Precip.df$Precip.julian.start[i]:Precip.df$Precip.julian.end[i]]], fun = mean, na.rm = T)
  precip.cov[i]<-extract(test.precip, Precip.df[i,], method = 'simple')
}

str(x1_proj.2018)
summary(precip.cov)
Precip.df
traps$precip<-scale(precip.cov)

#Organize and export traps
traps<-traps[,c(5,2,3,6,7,8,4)]
head(traps)

###Scale Julian Day
traps$JulianDay<-scale(traps$JulianDay)

write.table(traps@data, file = "2018TiogaTraps_nosegments.txt", quote = FALSE, row.names = FALSE)
#Make sure to add # in front of column names in .txt file
