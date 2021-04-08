#script for investigating and running ssf on all seabird data. Opening them, (old version: sub-sampling to one hourly). putting everything together
#March 11, 2021. Elham Nourani. Radolfzell am Bodensee
#update: don't subsample hourly. just subset to closest minute for gannets and nazca boobies (the sub-min resolutions that cause the C stack error.)
#update: one minute is still causing errors. try rounding up to 15 min and then removing duplicated points

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

#red-tailed tropicbird
rtt <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/Movebank/Red-tailed tropicbirds (Phaethon rubricauda) Round Island.csv", 
                stringsAsFactors = F,fileEncoding="latin1") %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  rename(TripID = comments) %>% 
  mutate(month = month(timestamp),
         year = year(timestamp))

#movebank files from Sophie
files <- list.files("data/From_Sophie/final_list_track_split", full.names = T)

data_ls <- sapply(files, function(x) mget(load(x)), simplify = TRUE)

#append rtt
data_ls$red_tailed_tropicbird <- rtt


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
  
  if(!("location.lat" %in% colnames(x))){
    x <- x %>% 
      rename(location.lat = Latitude,
             location.long = Longitude)
  }
  
  
  x[x$study.name == "Foraging ecology of masked boobies (Sula dactylatra) in the world’s largest “oceanic desert”", 
    "individual.taxon.canonical.name"] <- "Sula dactylatra"
  
  x
  
})

#common column names:
cols <- Reduce(intersect, lapply(data_ls, colnames))
  

data_df <- data_ls %>% 
  map(dplyr::select, cols) %>% 
  reduce(rbind) %>% 
  rename(sci_name = individual.taxon.canonical.name,
         indID = individual.local.identifier) %>% 
  mutate(sci_name = as.character(fct_recode(sci_name, 'Fregata magnificens' = "Fregata")),
         month = month(timestamp),
         year = year(timestamp),
         indID = as.character(indID)) %>% 
  rowwise() %>% 
  mutate(TripID = paste(year, indID, as.character(TripID))) %>% 
  ungroup()


nb_15min <- data_df %>% 
  filter(sci_name == "Sula granti") %>% 
  mutate(dt_15min = round_date(timestamp, "15 minutes")) 

#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(nb_15min$TripID),timestamps = nb_15min$dt_15min,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows

nb_15min <- nb_15min[-rows_to_delete,] 


data_df <- data_df %>% 
  filter(sci_name != "Sula granti") %>% 
  full_join(nb_15min)


#make sure all track IDs are unique
data_df %>% 
  group_by(TripID) %>% 
  summarise(n = n_distinct(indID)) %>%
  filter(n > 1) #these are not unique. so, paste the ind name with the trip ID and year to make it unique
  

save(data_df, file = "R_files/mv_incl_nb15min_rtt_df.RData")

# Peter Ryan data prep ######
#peter ryan data from Sophie

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/Peter_Ryan_data_annotated_SplitTrip.Rdata") #PR_data_split
#no prep needed

# Gremillet data prep ######
cg <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/David Gremillet/DavidGremillet_CapeGannet_AlgoaBay/CapeGannet-GPS-AlgoaBay-DavidGremillet-AllYears.csv") %>% 
  mutate(common_name = "Cape gannet",
         scientific_name = "Morus capensis")
ng <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/David Gremillet/DavidGremillet_NorthernGannet_IleRouzic/Gannet-GPS-Rouzic-DavidGremillet-AllYears.csv") %>% 
  mutate(common_name = "Northern gannet",
         scientific_name = "Morus bassanus",
         DateGMT = as.character(as.Date(DateGMT, format = "%d/%m/%Y")))


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

#make sure tripIDs are unique
gannets %>% 
  group_by(TrackId) %>% 
  summarise(n = n_distinct(BirdId)) %>% 
  filter(n > 1)

#subset to minutely :p

load("R_files/gannets_ls.RData")

gannets_15min <- gannets %>%
  rename(TripID = TrackId,
         indID = BirdId) %>%
  mutate(dt_15min = round_date(timestamp, "15 minutes")) 

#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(gannets_15min$TripID),timestamps = gannets_15min$dt_15min,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows

gannets_15min <- gannets_15min[-rows_to_delete,] 



save(gannets_15min, file = "R_files/gannets_15min.RData")


# step 2: put everything together # ------------------------------------------

#open data ####
all_files <- list("R_files/mv_incl_nb15min_rtt_df.RData",
              "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/Peter_Ryan_data_annotated_SplitTrip.Rdata", #PR_data_split
              "R_files/gannets_15min.RData")


data_ls_all <- sapply(all_files, function(x) mget(load(x)), simplify = TRUE)


#make sure all have the same column names 
data_ls_all <- lapply(data_ls_all,function(x){
  
  if(!("location.lat" %in% colnames(x))){
    x <- x %>% 
      rename(location.lat = Latitude,
             location.long = Longitude)
  }
  
  if(!("indID" %in% colnames(x))) {
    x <- x %>% 
      mutate(indID = NA)
  }
  
  if(!("FlyingSitting" %in% colnames(x))) {
    x <- x %>% 
      mutate(FlyingSitting = NA)
  }
  
  if(!("sci_name" %in% colnames(x)) & "scientific_name" %in% colnames(x)) {
    x <- x %>% 
      rename(sci_name = scientific_name)
  }
  
  x <- x %>% 
    mutate_at(c("TripID", "indID", "sci_name", "FlyingSitting"), as.character)
  
  
  x
  
})



#common column names:
cols <- Reduce(intersect, lapply(data_ls_all, colnames))

#put all data together

data_df_all <- data_ls_all %>% 
  map(dplyr::select, cols) %>% 
  reduce(rbind)

save(data_df_all, file = "R_files/all_spp_df_15min.RData")
