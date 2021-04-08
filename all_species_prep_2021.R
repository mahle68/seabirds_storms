#script for investigating and running ssf on all seabird data. second batch of species
#March 11, 2021. Elham Nourani. Radolfzell am Bodensee



library(tidyverse)
library(lubridate)
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


species <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/final_datasets.csv")


# Movebank data prep ######

#movebank files from Sophie
files <- list.files("data/From_Sophie/final_list_track_split", full.names = T)
lapply(files, load,.GlobalEnv)

data_ls <- sapply(files, function(x) mget(load(x)), simplify = TRUE)

#make sure all have timestamp 
data_ls <- lapply(data_ls,function(x){
  
  if("timestamp" %in% colnames(x)){
    x$timestamp <- as.POSIXct(strptime(x$timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")
  } else {
    x$timestamp <- as.POSIXct(strptime(x$date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")
  }
  
  if(!("FlyingSitting" %in% colnames(x))) {
    x$FlyingSitting <- NA
  }
  
  x
  
})

#common column names:
cols <- Reduce(intersect, lapply(data_ls, colnames))
  

data_df <- data_ls %>% 
  map(dplyr::select, cols) %>% 
  reduce(rbind) %>% 
  rename(sci_name = individual.taxon.canonical.name,
         indID = individual.local.identifier) %>% 
  mutate(sci_name = fct_recode(sci_name, 'Fregata magnificens' = "Fregata"),
         month = month(timestamp),
         year = year(timestamp))


save(data_df, file =  "R_files/move_data_df.RData")


# one-hourly subsample

load("R_files/move_data_df.RData")
#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(data_df$TripID),timestamps = data_df$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows

data_df <- data_df[-rows_to_delete,] 

#split to one-hourly trips
#check current sampling rate
move_ls <- lapply(split(data_df,data_df$study.name),function(x){
  x <- x %>%
    arrange(TripID, timestamp) %>% 
    as.data.frame()
  mv <- move(x = x$Longitude, y = x$Latitude, time = x$timestamp, data = x, animal = x$TripID, proj = wgs)
  mv
  
})


str(lapply(move_ls, timeLag, units = "hours"))

save(move_ls, file =  "R_files/move_data_mv_ls.RData")


# spp other than nazca
load("R_files/move_data_mv_ls.RData")


#thin the tracks. one hourly itnernavls
mycl <- makeCluster(detectCores() - 4) 

clusterExport(mycl, "move_ls") #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(move)
  library(sp)
  library(tidyverse)

})

(start_time <- Sys.time())
sp_ls_1hr <- parLapply(mycl, move_ls[-7], function(species){
#sp_ls_1hr <- lapply(move_ls, function(species){
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

Sys.time() - start_time #10 min

stopCluster(mycl)

save(sp_ls_1hr, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/movels_1hr_30_no_nb.RData")



#nazca booby separately. thinTrackTime produces C stack errors. so try the quick and easy way
load("R_files/move_data_mv_ls.RData")

nb <- move_ls[[7]]  #i get an error trying to run this. so, filter out na values and sitting points
nb <- as.data.frame(move_ls[[7]])


nb_1hr <- nb %>% 
  arrange(timestamp) %>% 
  mutate(hour = hour(timestamp),
         date = as.Date(timestamp)) %>% 
  group_by(TripID, date, hour) %>% #year is already taken care of in the track ID
  slice(1) %>% 
  ungroup()
  

#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(nb_1hr$TripID),timestamps = nb_1hr$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows.

#convert to move
nb_1hr_mv <- move(x = nb_1hr$coords.x1,y = nb_1hr$coords.x2,time = nb_1hr$timestamp,data =nb_1hr,animal = nb_1hr$TripID, proj = wgs)


timeLag(nb_1hr_mv, units = "mins")

save(nb_1hr_mv, file = "R_files/nazcabooby_1hr.RData")



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
  arrange(timestamp) %>% 
  mutate(hour = hour(timestamp)) %>% 
  group_by(TrackId, DateGMT, hour) %>%
  slice(1) %>% 
  ungroup()

move_ls_gannets <- lapply(split(gannets_1hr,gannets_1hr$common_name),function(x){
  
  x <- x %>%
    arrange(TrackId, timestamp) %>% 
    as.data.frame()
  
  mv <- move(x = x$Longitude,y = x$Latitude,time = x$timestamp,data = x,animal = x$TrackId,proj = wgs)
  mv
})


str(lapply(move_ls_gannets, timeLag, units = "mins")) #almost hourly

save(gannets_1hr, file = "R_files/gannets_1hr.RData")
save(move_ls_gannets, file = "R_files/gannets_1hr_mv.RData")

# Red-tailed tropicbird ######
rtt <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/Movebank/Red-tailed tropicbirds (Phaethon rubricauda) Round Island.csv", 
                stringsAsFactors = F,fileEncoding="latin1") %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  rename(sci_name = individual.taxon.canonical.name,
         indID = individual.local.identifier,
         TripID = comments) %>% 
  mutate(month = month(timestamp),
         year = year(timestamp))

#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(rtt$TripID),timestamps = rtt$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows.

#convert to move
rtt_mv <- move(x = rtt$location.long,y = rtt$location.lat,time = rtt$timestamp,data =rtt,animal = rtt$TripID, proj = wgs)


timeLag(rtt_mv) #this is already one-hourly :D

save(rtt_mv, file = "R_files/redtailedtropicbird_1hr.RData")

