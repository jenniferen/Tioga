library(dplyr)
library(raster)
library(sp)
library(sf)

#setwd("~/OSU/ElkSCR/")
data<-read.csv("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2018TiogaTransectTracks.csv")
head(data)

#scale Julian Day column- do this later
#data$JulianDay<-scale(data$JulianDay)

###choose distance to separate by and extract those points. Choose segment distance 
#using histograms above. Going with 300 meters.
data.sp<-SpatialPointsDataFrame(coords = data[,8:9], 
                                data = data,
                                proj4string = CRS(as.character("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))

data.sp.utm<- spTransform(data.sp, CRS("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"))
data.sp<-data.sp.utm
data.sp$UTM_X<-coordinates(data.sp)[,1]
data.sp$UTM_Y<-coordinates(data.sp)[,2]

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
head(data.sp)

#Segment traps with a vector of segment lengths.
#seg.size<-seq(0,min(Effort$dist/2), 50)
##seg.size<-c(seg.size, min(Effort$dist/2))
#as.integer(seg.size)
#[1]   0  50 100 150 200 250 300 350 385

data.sp.utm<-data.sp

Transect<-unique(data.sp$Transect)

seg.size=c(50,100,150,200,250,300,350, 360)
seg.size.2018<-seg.size
df <- data.frame(matrix(ncol = 6, nrow = 0))
seg.list<-list()
traps.list<-list()

#Create list of transect segments for each segment length.
#Need to retain Julian day and effort covariates here.
for(a in seg.size){
  df <- data.frame(matrix(ncol = 6, nrow = 0))
  seg<-a
  for(i in 1:length(Transect)){
    Transect.name<-Transect[i]
    data<-subset(data.sp@data, data.sp@data$Transect == Transect.name)
    extract.points<-seq(seg, max(data$dist.sum), seg)-seg/2
    extract.points<-c(extract.points, max(data$dist.sum)-(max(data$dist.sum)-extract.points[length(extract.points)]+seg/2)/2)
    for(j in 1:length(extract.points)){
      X<-data$UTM_X[which.min(abs(data$dist.sum-extract.points[j]))]
      Y<-data$UTM_Y[which.min(abs(data$dist.sum-extract.points[j]))]
      ObjectID<-data$OBJECTID[which.min(abs(data$dist.sum-extract.points[j]))]
      dist.sum<-data$dist.sum[which.min(abs(data$dist.sum-extract.points[j]))]
      JulianDay<-min(data$JulianDay)
      Effort<-max(data$dist.sum)
      df.j<-data.frame(Transect.name, extract.points[j], X, Y, ObjectID, dist.sum, JulianDay, Effort)
      colnames(df.j)<-c("Transect", "segment", "X", "Y", "ObjectID", "dist.sum", "JulianDay", "Effort")
      df<-rbind(df, df.j)
    }
  }
  df$TrapID<-1:nrow(df)
  df$char<-rep("/", nrow(df))
  df<-df[!duplicated(df$ObjectID),]
  seg.list[[which(seg.size == a)]]<-df
  df$Effort<-scale(df$Effort)
  traps.list[[which(seg.size == a)]]<-df[,c(9,3,4,10,8,7)]
}

traps.list.2018<-traps.list
#Read in coordinates and sex of identified elk
TiogaElk2018<-read.csv("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2018TiogaElkCoordinatesUTMSex.csv", header=TRUE)
head(TiogaElk2018)
colnames(TiogaElk2018)[5:6]<-c("UTM_E", "UTM_N")
TiogaElk2018.sp<-SpatialPointsDataFrame(coords = TiogaElk2018[,5:6]/1000,
                                        data = TiogaElk2018,
                                        proj4string = CRS(as.character("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")))

#Calculate distances between individuals and nearest trap
library(sf)
seg.closest<-list()
distances<-list()
for(k in 1:length(traps.list)){
  a<-st_as_sf(TiogaElk2018, coords=c("UTM_E", "UTM_N"), crs="+proj=utm +zone=10 +datum=NAD83 
            +units=m +no_defs")
  b<-st_as_sf((traps.list[[k]]), coords=c("X", "Y"), crs="+proj=utm +zone=10 +datum=NAD83
            +units=m +no_defs")
  closest<-list() #creates a blank list to store data
  distance<-c()
  for(i in seq_len(nrow(a))){ #for each Elk in a
    closest[[i]]<-b[which.min( #write in each row of closest, b (Trap ID) which is the 
      st_distance(b, a[i,]) #minimum distance between all trap IDs and Elk i.
    ),]
    distance<-c(distance,which.min(st_distance(b, a[i,])))
  }
  seg.closest[[k]]<-closest
  distances[[k]]<-distance
}

seg.closest.2018<-seg.closest

#Extract vector of closest Trap IDs.
closestTrapID<-list()
for(i in 1:length(seg.closest.2018)){
  closestTrapID[[i]]<-c(sapply(seg.closest.2018[[i]], "[[", 1))
}

#Create capture history data frames and store them in list: TiogaElk2018.list
TiogaElk2018.list<-list()
for(i in 1:length(closestTrapID)){
  TiogaElk2018.list[[i]]<-TiogaElk2018
  TiogaElk2018.list[[i]]$TrapID<-closestTrapID[[i]]
  TiogaElk2018.list[[i]]$Session<-paste0("Tioga18", "_", seg.size[[i]], "Msegments")
  TiogaElk2018.list[[i]]$Occasion<-rep(1, nrow(TiogaElk2018.list[[i]]))
  TiogaElk2018.list[[i]]<-TiogaElk2018.list[[i]][,c(15,4,16,14,13)]
  colnames(TiogaElk2018.list[[i]])[2]<-"Individual"
}

#Calculate percentage of spatial recaptures in number of recaptures.
Ratio<-c()
for(i in 1:length(TiogaElk2018.list)){
  nonspatial.recaptures.sum<-0
  duplicate<-c()
  for(j in 1:length(unique(TiogaElk2018.list[[i]]$Individual))){
    nonspatial.recaptures<-duplicated(subset(TiogaElk2018.list[[i]], TiogaElk2018.list[[i]]$Individual == j)$TrapID)
    duplicate<-c(duplicate, nonspatial.recaptures)
    nonspatial.recaptures<-length(nonspatial.recaptures[nonspatial.recaptures==TRUE])
    nonspatial.recaptures.sum<-nonspatial.recaptures.sum+nonspatial.recaptures
  }
  recaptured.individuals<-duplicated(TiogaElk2018.list[[i]]$Individual) | duplicated(TiogaElk2018.list[[i]]$Individual, fromLast = TRUE)
  recaptures<-duplicated(TiogaElk2018.list[[i]]$Individual)
  recaptures.sum<-length(recaptures[recaptures==TRUE])
  spatial.recaptures<-recaptures.sum-nonspatial.recaptures.sum
  Ratio<-c(Ratio, spatial.recaptures/recaptures.sum)
  TiogaElk2018.list[[i]]$Recaptured.Individual<-recaptured.individuals
  TiogaElk2018.list[[i]]$Duplicate<-duplicate
}

#plot relationship between the percent of spatial recaptures and segment length.
ratios.2018<-data.frame(Seg.size = seg.size.2018, Ratio = Ratio)
plot(ratios.2018$Seg.size, ratios.2018$Ratio, main = "Proportion of Recaptures that were Spatial Recaptures \n vs. Transect Segment Size",
     xlab = "Transect Segment Size (meters)", ylab = "Spatial Recaptures/Total Recaptures",
     xaxt = "n", xlim = c(0, 400), ylim = c(0.5,1))
axis(1, at=seg.size.2018, labels=seg.size.2018)


dist.recaptures.list<-list()
#Calculate distance between spatial recaptures.
#Need vector of TRUE/FALSE, spatial recaptures.
for(i in 1:length(TiogaElk2018.list)){
  data<-subset(TiogaElk2018.list[[i]], TiogaElk2018.list[[i]]$Recaptured.Individual == "TRUE")
  dist.recaptures.vec<-c()
  for(j in unique(TiogaElk2018.list[[i]]$Individual)){
    ind.data<-subset(data, data$Individual==j)
    TrapIDs<-ind.data[!duplicated(ind.data$TrapID),]
    traps<-traps.list.2018[[i]][traps.list.2018[[i]]$TrapID %in% TrapIDs$TrapID,][,2:3]
    dist.recaptures<-c(dist(traps))
    dist.recaptures.vec<-c(dist.recaptures.vec,dist.recaptures)
  }
  dist.recaptures.list[[i]]<-dist.recaptures.vec[!dist.recaptures.vec==0]
}


hist(dist.recaptures.list[[1]], breaks = seq(0, 3000, 50), main = "Histogram of Distances between Spatial Recaptures", sub = "Transect Segment Length = 50 meters",
     xlab = "Distance (meters)", ylim = c(0,20))
hist(dist.recaptures.list[[2]], breaks = seq(0, 3000, 50), main = "Histogram of Distances between Spatial Recaptures", sub = "Transect Segment Length = 100 meters",
     xlab = "Distance (meters)", ylim = c(0,20))
hist(dist.recaptures.list[[3]], breaks = seq(0, 3000, 50), main = "Histogram of Distances between Spatial Recaptures", sub = "Transect Segment Length = 150 meters",
     xlab = "Distance (meters)", ylim = c(0,20))
hist(dist.recaptures.list[[4]], breaks = seq(0, 3000, 50), main = "Histogram of Distances between Spatial Recaptures", sub = "Transect Segment Length = 200 meters",
     xlab = "Distance (meters)", ylim = c(0,20))
hist(dist.recaptures.list[[5]], breaks = seq(0, 3000, 50), main = "Histogram of Distances between Spatial Recaptures", sub = "Transect Segment Length = 250 meters",
     xlab = "Distance (meters)", ylim = c(0,20))
hist(dist.recaptures.list[[6]], breaks = seq(0, 3000, 50), main = "Histogram of Distances between Spatial Recaptures", sub = "Transect Segment Length = 300 meters",
     xlab = "Distance (meters)", ylim = c(0,20))
hist(dist.recaptures.list[[7]], breaks = seq(0, 3000, 50), main = "Histogram of Distances between Spatial Recaptures", sub = "Transect Segment Length = 350 meters",
     xlab = "Distance (meters)", ylim = c(0,20))
hist(dist.recaptures.list[[8]], breaks = seq(0, 3000, 50), main = "Histogram of Distances between Spatial Recaptures", sub = "Transect Segment Length = 360 meters",
     xlab = "Distance (meters)", ylim = c(0,20))

#Delete non-spatial recaptures in TiogaElk2019.list
for(i in 1:length(TiogaElk2018.list)){
  TiogaElk2018.list[[i]]<-TiogaElk2018.list[[i]][TiogaElk2018.list[[i]]$Duplicate==FALSE,]
  TiogaElk2018.list[[i]]<-TiogaElk2018.list[[i]][,-ncol(TiogaElk2018.list[[i]])]
}

#Write out capture history files
for(i in 1:length(TiogaElk2018.list)){
  write.table(TiogaElk2018.list[[i]][,1:5], file = paste0("2018TiogaCH_",seg.size[i],"Msegments.txt"), row.names = FALSE, quote = FALSE)
}

#Add precipitation covariate to trap file
#Precipitation
#Use x1_proj.2018.crop from environment - created with PRISMScript

Precip.list<-list()
for(i in 1:length(traps.list.2018)){
  traps.list.2018[[i]]$Precip.julian.end<-traps.list.2018[[i]]$JulianDay+15
  traps.list.2018[[i]]$Precip.julian.start<-traps.list.2018[[i]]$JulianDay+1
  Precip.list[[i]]<-as.data.frame(cbind(traps.list.2018[[i]]$X, 
                                        traps.list.2018[[i]]$Y,
                                        traps.list.2018[[i]]$Precip.julian.end,
                                        traps.list.2018[[i]]$Precip.julian.start))
  names(Precip.list[[i]])<-c("X", "Y", "Precip.julian.end", "Precip.julian.start")
}

for(i in 1:length(Precip.list)){
  Precip.df<-Precip.list[[i]]
  coordinates(Precip.df)<-Precip.df[,c(1,2)]
  precip.cov<-rep(NA, nrow(Precip.df))
  test.precip<-NULL
  for(j in 1:nrow(Precip.df)){
    test.precip<-calc(x1_proj.2018.crop[[Precip.df$Precip.julian.start[j]:Precip.df$Precip.julian.end[j]]], fun = mean, na.rm = T)
    precip.cov[j]<-extract(test.precip, Precip.df[j,], method = 'simple')
  }
  traps.list.2018[[i]]$precip<-scale(precip.cov)
}


#Clean up trap file, Scale trap coordinates, scale Julian Day,
#and write out trap files.
for(i in 1:length(traps.list.2018)){
  traps.list.2018[[i]]<-traps.list.2018[[i]][,c(1:6,9)]
  traps.list.2018[[i]]$X<-traps.list.2018[[i]]$X/1000
  traps.list.2018[[i]]$Y<-traps.list.2018[[i]]$Y/1000
  traps.list.2018[[i]]$JulianDay<-scale(traps.list.2018[[i]]$JulianDay)
  write.table(traps.list.2018[[i]], file = paste0("2018TiogaTrapfile_", seg.size[i], "Msegments.txt"), row.names = FALSE, quote = TRUE)
}
