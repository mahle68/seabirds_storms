#script for investigating and running ssf on all seabird data. Peter Ryan edition. and Gremiellet (all studies that Sophie didn't do the flyingsitting assignment)
#April 6, 2021. Elham Nourani. Radolfzell am Bodensee
#flyigsitting: threshold is 3 km/h for one-hourly data (for wandering albatrosses).

#look at section 9.5 in Virgilio's book for drawing smooths for each level of a factor variable (doesnt work for binned data)

library(tidyverse)
library(lubridate)
library(move)
library(mapview)
library(sf)
library(sp)
library(move)
library(circular)
library(CircStats)
library(fitdistrplus)
library(corrr)
library(INLA)

wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")

NCEP.loxodrome.na <- function (lat1, lat2, lon1, lon2) {
  deg2rad <- pi/180
  acot <- function(x) {
    return(atan(1/x))
  }
  lat1 <- deg2rad * lat1
  lat2 <- deg2rad * lat2
  lon1 <- deg2rad * lon1
  lon2 <- deg2rad * lon2
  deltaLon <- lon2 - lon1
  pi4 <- pi/4
  Sig1 <- log(tan(pi4 + lat1/2))
  Sig2 <- log(tan(pi4 + lat2/2))
  deltaSig <- Sig2 - Sig1
  if (deltaLon == 0 && deltaSig > 0) {
    head <- 0
  }
  else if (deltaLon == 0 && deltaSig < 0) {
    head <- 180
  }
  else if (deltaSig == 0 && deltaLon > 0) {
    head <- 90
  }
  else if (deltaSig == 0 && deltaLon < 0) {
    head <- 270
  }
  else if (deltaSig < 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig < 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig > 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi
  }
  else if (deltaSig > 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 360
  }
  else {
    head <-NA}
  return(head)
}
source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")


# Peter Ryan data prep ######
#peter ryan data from Sophie
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/Peter_Ryan_data_annotated_SplitTrip.Rdata") #PR_data_split



#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(PR_data_split$TripID),timestamps = PR_data_split$date_time,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows

PR_data_split <- PR_data_split[-rows_to_delete,] 

#create move object
move_ls <- lapply(split(PR_data_split,PR_data_split$common_name),function(x){
  
  x <- x %>%
    arrange(TripID, date_time) %>% 
    as.data.frame()
  
  mv <- move(x = x$location.long,y = x$location.lat,time = x$date_time,data = x,animal = x$TripID,proj = wgs)
  mv
})

#original sampling frequency
str(lapply(move_ls, timeLag, units = "hours")) #almost hourly

#sub-sample to hourly 
(start_time <- Sys.time())
sp_ls_1hr <- lapply(move_ls, function(species){
  
  lapply(split(species), function(track){
    track_th <- track %>%
      thinTrackTime(interval = as.difftime(1, units='hours'),
                    tolerance = as.difftime(30, units='mins'))#the unselected bursts are the large gaps between the selected ones
    track_th$selected <- c(as.character(track_th@burstId),NA) #assign selected as a variable
    track_th <- as(track_th,"Move") #convert back to a move object (from a moveburst)
    track_th <- track_th[is.na(track_th$selected) | track_th$selected == "selected",]
    track_th
  }) %>% 
    reduce(rbind) #gets converted to a spatialpointsdf
  
})

Sys.time() - start_time # < 1 min


save(sp_ls_1hr, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/movels_1hr_30_PR.RData")


# Gremillet data prep ######
rtt <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/Movebank/Red-tailed tropicbirds (Phaethon rubricauda) Round Island.csv", 
                stringsAsFactors = F,fileEncoding="latin1")
cg <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/David Gremillet/DavidGremillet_CapeGannet_AlgoaBay/CapeGannet-GPS-AlgoaBay-DavidGremillet-AllYears.csv") %>% 
  mutate(common_name = "Cape gannet",
           scientific_name = "Morus capensis")
ng <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/David Gremillet/DavidGremillet_NorthernGannet_IleRouzic/Gannet-GPS-Rouzic-DavidGremillet-AllYears.csv") %>% 
  mutate(common_name = "Northern gannet",
         scientific_name = "Morus bassanus",
         DateGMT = as.character(as.Date(DateGMT, format = "%d/%m/%Y")))# %>% 
  #rowwise() %>% 
  #mutate(DateGMT = str_replace_all(DateGMT, "/", "-")) %>% 
  #ungroup()


gannets <- cg %>% 
  full_join(ng) %>%
  rowwise() %>% 
  mutate(timestamp = paste(DateGMT,TimeGMT, sep = " ")) %>% 
  ungroup() %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))


#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(gannets$TrackId),timestamps = gannets$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows

gannets <- gannets[-rows_to_delete,] 


save(gannets, file = "R_files/gannets_ls.RData")

#create one move stack

load("R_files/gannets_ls.RData")


gannets_1hr <- gannets %>% 
  mutate(hour = hour(timestamp)) %>% 
  group_by(TrackId, DateGMT, hour) %>%
  slice(1) %>% 
  ungroup()
  
move_ls <- lapply(split(gannets_1hr,gannets_1hr$common_name),function(x){
  
  x <- x %>%
    arrange(TrackId, timestamp) %>% 
    as.data.frame()
  
  mv <- move(x = x$Longitude,y = x$Latitude,time = x$timestamp,data = x,animal = x$TrackId,proj = wgs)
  mv
})


str(lapply(move_ls, timeLag, units = "mins")) #almost hourly



sample <- gannets[1:100,]

sample  %>% 
  mutate(hour = hour(timestamp)) %>% 
  group_by(TrackId, DateGMT, hour) %>%
  slice(1)

# Apr 7 -------------------------------------------------------------
mv <- gannets %>% 
  arrange(common_name, TrackId, timestamp) 
 
mv <- move(x = mv$Longitude,y = mv$Latitude,time = mv$timestamp,data = mv, animal = mv$TrackId,proj = wgs)


sp_data <- NULL

# try a for loop
for(i in split(mv)){
  
  thinned_mv <- i %>%
    thinTrackTime(interval = as.difftime(60, units='mins'),
                  tolerance = as.difftime(30, units='mins')) #the unselected bursts are the large gaps between the selected ones
  track_th$selected <- c(as.character(track_th@burstId),NA) #assign selected as a variable
  track_th <- as(track_th,"Move") #convert back to a move object (from a moveburst)
  track_th <- track_th[is.na(track_th$selected) | track_th$selected == "selected",]
  track_th
  
}


#g_sp <- lapply(split(mv), function(x){
#  track_th <- 
#}) %>% 
  reduce(rbind) #gets converted to a spatialpointsdf



move_ls <- lapply(split(gannets,gannets$common_name),function(x){
  
  x <- x %>%
    arrange(TrackId, timestamp) %>% 
    as.data.frame()
  
  mv <- move(x = x$Longitude,y = x$Latitude,time = x$timestamp,data = x,animal = x$TrackId,proj = wgs)
  mv
})

save(move_ls, file = "R_files/gannets_mv_ls.RData")

#original sampling frequency
str(lapply(move_ls, timeLag, units = "hours")) #almost hourly

#sub-sample to hourly 

load("R_files/gannets_mv_ls.RData")

mycl <- makeCluster(detectCores() - 4) 

clusterExport(mycl, "move_ls") #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(move)
  library(sp)
  library(tidyverse)
  
})
(start_time <- Sys.time())
#sp_ls_1hr <- parLapply(mycl, move_ls, function(species){
sp_ls_1hr <- lapply(move_ls, function(species){
  
  ng_sp <- lapply(split(species)[1:3], function(x){
    track_th <- x %>%
      thinTrackTime(interval = as.difftime(60, units='mins'),
                    tolerance = as.difftime(30, units='mins')) #the unselected bursts are the large gaps between the selected ones
    track_th$selected <- c(as.character(track_th@burstId),NA) #assign selected as a variable
    track_th <- as(track_th,"Move") #convert back to a move object (from a moveburst)
    track_th <- track_th[is.na(track_th$selected) | track_th$selected == "selected",]
    track_th
  }) %>% 
    reduce(rbind) #gets converted to a spatialpointsdf
  
})

Sys.time() - start_time # < 1 min


##### think the track manually

gannets_1hr <- gannets %>% 
  group_by(common_name, TrackId) %>%
  arrange(timestamp) %>% 
  rowwise() %>% 
  mutate(timediff = difftime(lag(timestamp,1),timestamp,units = "mins"))


(start_time <- Sys.time())
gannets_1hr <- gannets %>% 
  group_by(common_name, TrackId) %>%
  arrange(timestamp) %>%
  mutate(as.duration(lag(timestamp) %--% timestamp))
Sys.time() - start_time # 0.9 sec


(start_time <- Sys.time())
gannets_1hr <- gannets %>% 
  group_by(common_name, TrackId) %>%
  arrange(timestamp) %>%
  mutate(time_lag = as.duration(lag(timestamp) %--% timestamp)) %>% 
  mutate(selected = cusum())
Sys.time() - start_time # 0.9 sec

