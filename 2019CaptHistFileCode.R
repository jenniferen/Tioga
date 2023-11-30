##########################################
####Build Capture History File############
##########################################


setwd("~/OSU/ElkSCR")

#Read in trap file
traps<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2019TiogaTraps_nosegments.txt")
head(traps)
#colnames(traps)<-c("TransectID", "X", "Y", "sep", "Effort", "Julian", "Precip")
colnames(traps)<-c("TransectID", "X", "Y")
head(traps)

#Read in coordinates and sex of identified elk
TiogaElk2019<-read.table("C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2019TiogaElkCoordinatesUTMSex.txt", header=TRUE)
head(TiogaElk2019)

#Adding Unknown to Sex column for individuals with unknown Sex
levels<-levels(factor(TiogaElk2019$Sex)) ###Making NAs Unknowns
levels[length(levels)+1]<-"U"
TiogaElk2019$Sex<-factor(TiogaElk2019$Sex, levels = levels)
TiogaElk2019$Sex[is.na(TiogaElk2019$Sex)]<- "U"
TiogaElk2019

#Changing collection date to Julian day to gather precipitation data from PRISM script for Habitat Mask
TiogaElk2019$Collection_Date<-as.Date(TiogaElk2019$Collection_Date, "%m/%d/%y")

require(lubridate)
TiogaElk2019$Collection_Date = yday(TiogaElk2019$Collection_Date)
head(TiogaElk2019)

TiogaElk2019$Collection_Date<-TiogaElk2019$Collection_Date - 74
head(TiogaElk2019)

saveRDS(TiogaElk2019, file = "C:/Users/nelsonj7/Desktop/MSUComputer/OSU/ElkSCR/Tioga/2019TiogaElkCoordinatesUTMSexJulian.rds")

####Putting in Session, Animal ID, occassion, X, Y format
TiogaElk2019$Session<-rep("Tioga19", nrow(TiogaElk2019))
TiogaElk2019$Occasion<-rep(1, nrow(TiogaElk2019))
head(TiogaElk2019)
TiogaElk2019<-TiogaElk2019[,c(6,4,7,1,2,5,3)]
head(TiogaElk2019)
colnames(TiogaElk2019)<-c("Session", "AnimalID", "Occasion", "X", "Y", "Sex", "Collection Date")

TiogaElk2019.NoDate<-TiogaElk2019[,c(1:6)]

write.csv(TiogaElk2019, "2019TiogaElk.csv", row.names=FALSE, quote = FALSE)


#Calculate distances between individuals and nearest trap
library(sf)
a<-st_as_sf(TiogaElk2019.NoDate, coords=c("X", "Y"), crs="+proj=utm +zone=10 +datum=NAD83 
            +units=m +no_defs")
b<-st_as_sf(traps.list[[1]], coords=c("X", "Y"), crs="+proj=utm +zone=10 +datum=NAD83
            +units=m +no_defs")
closest<-list() #creates a blank list to store data
for(i in seq_len(nrow(a))){ #for each Elk in a
  closest[[i]]<-b[which.min( #write in each row of closest, b (Trap ID) which is the 
    st_distance(b, a[i,]) #minimum distance between all trap IDs and the current Elk, i.
  ),]
}

##Calculate distances between individuals
huzzah<-st_distance(a[-1,],a[-nrow(a),],by_element=TRUE)
write.csv(huzzah, file = "Tioga2019DistancesbetweenIndividualsmeters.csv", row.names = FALSE)
###

closestTrapID<-sapply(closest, "[[", 1) #Extracting TrapIDs from list
closestTrapID #viewing TrapIDs

#####
c<-c(0, huzzh
     d<-cbind(c, huzzah)
     test<-cbind(huzzah, closestTrapID)
###

head(TiogaElk2019)
TiogaElk2019$TrapID<-closestTrapID #creating column in Tioga Elk
head(TiogaElk2019)

Tioga2019.capthist<-TiogaElk2019[,c(1,2,3,8,6)]
names(Tioga2019.capthist)<-c("Session", "ID", "Occasion", "Detector", "Sex")
head(Tioga2019.capthist)
write.table(Tioga2019.capthist, file = "2019TiogaCH.txt", row.names = FALSE, quote = FALSE)

#Find this file, examine individuals and delete
#whenever an individual was detected at a trap more than once. Delete first column (doesn't have a header)
#Then save as a tab-deliminated txt. Open txt and add # and one space before Session.

###Scale trap coordinates
traps$X<-traps$X/1000
traps$Y<-traps$Y/1000

write.table(traps, file = "2019TiogaTraps_nosegments.txt", row.names = FALSE, quote = FALSE)

##Of the individuals that were recaptures, what percentage of those individuals were spatially recaptured?
head(Tioga2019.capthist)
recaptured.individuals<-Tioga2019.capthist[with(Tioga2019.capthist, ave(TrapID,TrapID,FUN=length))>1, ]
write.csv(recaptured.individuals, "2019TiogaRecaptured.Individuals.csv", row.names = FALSE, quote = FALSE)
##Individuals that were recaptured but not spatially:  20, 33, 38, 40 = 4 individuals of 55. Therefore, 51/55 
#individuals that were recaptured were spatially recaptured (93%)