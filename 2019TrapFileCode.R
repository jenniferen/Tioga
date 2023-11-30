library(dplyr)
library(raster)
library(sp)
library(sf)

#setwd("~/OSU/ElkSCR/")
data<-read.csv("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2019TiogaTracksCombined.csv")
head(data)

###choose distance to separate by and extract those points. Choose segment distance 
#using histograms above. Going with 300 meters.
data.sp<-SpatialPointsDataFrame(coords = data[,c(4:5)], 
                                data = data,
                                proj4string = CRS(as.character("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")))

data.sp$dist<-NA

#Calculate the distance between consecutive points by transect, place NA if first point of transect.
for(i in 1:nrow(data.sp)) {
  if(i > 1 && data.sp$Transect[i] == data.sp$Transect[i-1]) {
    data.sp$dist[i] <- sqrt(((data.sp$UTM_X[i] - data.sp$UTM_X[i-1]) ^ 2) + ((data.sp$UTM_Y[i] - data.sp$UTM_Y[i-1]) ^ 2))
  } else {
    data.sp$dist[i] <- NA
  }
}

#Calculate cumulative distance by transect.
data.sp$dist.sum<-NA
for(i in 1:nrow(data.sp)) {
  if(i > 1 && data.sp$Transect[i] == data.sp$Transect[i-1]) {
    data.sp$dist.sum[i] <- data.sp$dist.sum[i-1] + data.sp$dist[i]
  } else {
    data.sp$dist.sum[i] <- 0
  }
}

data.sp$dist.sum
max(data.sp$dist.sum)

#What do I want to do?
#I want to produce a data frame with two vectors: 1. segment size of transects. 2. ratio of recaptures that were spatial recaptures.
#What are the values of segment size? The max transect distance -> max(data.sp$dist.sum) <- is 6468 meters and minimum is 771 meters 
#-> max(Effort$dist) <- Every km. would be too big for the smaller ones. For a maximum let's use half of the minimum distance: 385.5,
#then decrease in increments of min(Effor$dist/)

seg.size<-seq(0,min(Effort$dist/2), 50)
seg.size<-c(seg.size, min(Effort$dist/2))
as.integer(seg.size)
#[1]   0  50 100 150 200 250 300 350 385

data.sp.utm<-data.sp
#Now, I want to extract points that occur in increments of x-specified meters. Using wich.min(abs(x$dist.sum)) may cause duplicates. 
#Will delete those which will result in less traps per transect when it occurs. Keep in mind when examining ratios.

#Want to make a data frame with four columns: 1. Transect Name. 2. segment size. 3/4. coordinates of point that occurs closest to the segment length 
#(50 meter segment size, trap point that occurs closest to 50, 100, 150, 200 meters walked...that occur in each increment.
Transect<-unique(data.sp$Transect)
#For each subset of transects, calculate max dist.sum, divide by seg.size, get integer value as number of times segment size can go in.
#extract values = seq(0, max dist.sum, seg.size)
#make a dataframe with columns: 1. Transect, 3. segment point, 2. extract points, 3. X, 4. Y

#seg.size=c(50,100,150,200,250,300,350,385)
seg.size = 300
df <- data.frame(matrix(ncol = 6, nrow = 0))
#seg.list<-list()

#for(a in 1:length(seg.size)){
#  df <- data.frame(matrix(ncol = 6, nrow = 0))
#######ADD seg.size[a]
for(i in 1:length(Transect)){
  Transect.name<-Transect[i]
  data<-subset(data.sp@data, data.sp@data$Transect == Transect.name)
  extract.points<-seq(0, max(data$dist.sum), seg.size) 
  for(j in 1:length(extract.points)){
    X<-data$UTM_X[which.min(abs(data$dist.sum-extract.points[j]))]
    Y<-data$UTM_Y[which.min(abs(data$dist.sum-extract.points[j]))]
    ObjectID<-data$OBJECTID[which.min(abs(data$dist.sum-extract.points[j]))]
    dist.sum<-data$dist.sum[which.min(abs(data$dist.sum-extract.points[j]))]
    df.j<-data.frame(Transect.name, extract.points[j], X, Y, ObjectID, dist.sum)
    colnames(df.j)<-c("Transect", "segment", "X", "Y", "ObjectID", "dist.sum")
    df<-rbind(df, df.j)
  }
}

df<-df[!duplicated(df$ObjectID),]
#  seg.list[[a]]<-df
#}
traps<-df[,c(1,3,4)]
head(traps)

#Need TrapID, X, Y, character, effort, Julian, precipitation, and # of observers.
TrapID<-seq(1:nrow(traps)) #Trap ID
traps$TrapID<-TrapID
traps$char<-rep("/", nrow(traps))

#effort covariate: want to take the max dist of each transect, scale and add to each point.
Effort<-data.frame(matrix(ncol = 2, nrow = length(Transects)))
colnames(Effort)<-c("Transect", "dist")
Transects<-unique(data.sp$Transect)
for(i in 1:length(Transects)){
  sub<-subset(data.sp, data.sp$Transect == Transects[i])
  Effort$Transect[i]<- as.character(Transects[i])
  Effort$dist[i]<-as.numeric(max(sub$dist.sum))
}

min(Effort$dist)
max(Effort$dist)
head(Effort)
Effort$dist.scale<-scale(Effort$dist)

traps$effort <- Effort$dist.scale[match(traps$Transect, Effort$Transect)]
head(traps)


#Precipitation
Precip.julian.end<-traps$JulianDay+15
Precip.julian.start<-traps$JulianDay+1
Precip.df<-as.data.frame(cbind(traps$X, traps$Y, Precip.julian.end, Precip.julian.start))
#
coordinates(Precip.df)<-Precip.df[,c(1,2)]
precip.cov<-rep(NA, nrow(Precip.df))
test.precip<-NULL
x1_proj.2019<-raster("C:/Users/jnelson/Documents/OSU/ElkSCR/PRISM/TiogaPrecip_2019_daily.tif")


for(i in 1:nrow(Precip.df)){
  test.precip<-calc(x1_proj.crop.2019[[Precip.df$Precip.julian.start[i]:Precip.df$Precip.julian.end[i]]], fun = mean, na.rm = T)
  precip.cov[i]<-extract(test.precip, Precip.df[i,], method = 'simple')
}

traps$precip<-scale(precip.cov)

#Number of observers
#Tioga2019obs<-c(1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,2,2,2,2,
#                2,2,2,2,2,2,2,2,2,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,
#                3,2,3,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
#                2,1,2)

#obs.df<-as.data.frame(cbind(unlist(Transects), Tioga2019obs))
#names(obs.df)<-c("Transect", "Tioga2019obs")
#traps$Observers<-obs.df$Tioga2019obs[match(traps$Transect, obs.df$Transect)]
#traps$Observers<-as.numeric(traps$Observers)
#traps$Observers<-scale(traps$Observers)
#head(traps)

#Organize and export traps
traps<-traps[,c(5,2,3,6,7,8,4)]
head(traps)

###Scale Julian Day
traps$JulianDay<-scale(traps$JulianDay)
head(traps)

write.table(traps, file = "2019TiogaTraps_nosegments.txt", quote = FALSE, row.names = FALSE)
#Make sure to add # in front of column names in .txt file