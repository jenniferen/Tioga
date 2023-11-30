##########################################
####Build Capture History File############
##########################################


setwd("~/OSU/ElkSCR/Tioga")

#Read in trap file
traps<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2018TiogaTraps_nosegments.txt")
head(traps)
#colnames(traps)<-c("TransectID", "X", "Y", "sep", "Effort", "Precip", "Julian")
colnames(traps)<-c("TransectID", "X", "Y")
head(traps)

#Read in coordinates and sex of identified elk
TiogaElk2018<-read.csv("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2018TiogaElkCoordinatesUTMSex.csv", header=TRUE)
head(TiogaElk2018)
dst<-dist(TiogaElk2018[,c(5,6)])



#Adding Unknown to Sex column for individuals with unknown Sex
levels<-levels(factor(TiogaElk2018$Sex)) ###Making NAs Unknowns
TiogaElk2018$Sex<-factor(TiogaElk2018$Sex, levels = levels)
head(TiogaElk2018)

#Changing collection date to Julian day to gather precipitation data from PRISM script for Habitat Mask
TiogaElk2018$CollectionDate<-as.Date(TiogaElk2018$CollectionDate, "%m/%d/%y")

require(lubridate)
TiogaElk2018$CollectionDate = yday(TiogaElk2018$CollectionDate)
head(TiogaElk2018)

TiogaElk2018$CollectionDate<-TiogaElk2018$CollectionDate - 89
head(TiogaElk2018)

saveRDS(TiogaElk2018, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2018TiogaElkCoordinatesUTMSexJulian.rds")

####Putting in Session, Animal ID, occassion, X, Y format
TiogaElk2018$Session<-rep("Tioga18", nrow(TiogaElk2018))
TiogaElk2018$Occasion<-rep(1, nrow(TiogaElk2018))
head(TiogaElk2018)
TiogaElk2018<-TiogaElk2018[,c(14,4,15,5,6,13,12)]
head(TiogaElk2018)
colnames(TiogaElk2018)<-c("Session", "AnimalID", "Occasion", "X", "Y", "Sex", "Collection Date")

TiogaElk2018.NoDate<-TiogaElk2018[,c(1:6)]

write.csv(TiogaElk2018, "2018TiogaElk.csv", row.names=FALSE, quote = FALSE)


#Calculate distances between individuals and nearest trap
library(sf)
a<-st_as_sf(TiogaElk2018.NoDate, coords=c("X", "Y"), crs="+proj=utm +zone=10 +datum=NAD83 
            +units=m +no_defs")
b<-st_as_sf(traps, coords=c("X", "Y"), crs="+proj=utm +zone=10 +datum=NAD83
            +units=m +no_defs")
closest<-list() #creates a blank list to store data
for(i in seq_len(nrow(a))){ #for each Elk in a
  closest[[i]]<-b[which.min( #write in each row of closest, b (Trap ID) which is the 
    st_distance(b, a[i,]) #minimum distance between all trap IDs and the current Elk, i.
  ),]
}

##Calculate distances between individuals
huzzah<-st_distance(a[-1,],a[-nrow(a),],by_element=TRUE)
write.csv(huzzah, file = "Tioga2018DistancesbetweenIndividualsmeters.csv", row.names = FALSE)
###

closestTrapID<-sapply(closest, "[[", 1) #Extracting TrapIDs from list
closestTrapID #viewing TrapIDs

#####
c<-c(0, huzzh
     d<-cbind(c, huzzah)
     test<-cbind(huzzah, closestTrapID)
###

head(TiogaElk2018)
TiogaElk2018$TrapID<-closestTrapID #creating column in Tioga Elk
head(TiogaElk2018)

Tioga2018.capthist<-TiogaElk2018[,c(1,2,3,8,6)]
names(Tioga2018.capthist)<-c("Session", "ID", "Occasion", "Detector", "Sex")
head(Tioga2018.capthist)
write.table(Tioga2018.capthist, file = "2018TiogaCH.txt", row.names = FALSE, quote = FALSE)

#Find this file, examine individuals and delete
#whenever an individual was detected at a trap more than once. Delete first column (doesn't have a header)
#Then save as a tab-deliminated txt. Open txt and add # and one space before Session.

###Scale coordinates
traps$X<-traps$X/1000
traps$Y<-traps$Y/1000

head(traps)

write.table(traps, file = "2018TiogaTraps_nosegments.txt", quote = FALSE, row.names = FALSE)

##Of the individuals that were recaptures, what percentage of those individuals were spatially recaptured?
head(Tioga2018.capthist)
recaptured.individuals<-Tioga2018.capthist[with(Tioga2018.capthist, ave(ID,ID,FUN=length))>1, ]
write.csv(recaptured.individuals, "2018TiogaRecaptured.Individuals.csv", row.names = FALSE, quote = FALSE)
##Individuals that were recaptured but not spatially:  none of 41 recaptured individuals. Therefore, 41/41 
#individuals that were recaptured were spatially recaptured (100%)