#script for investigating and running ssf on all seabird data. Opening them, (old version: sub-sampling to one hourly). putting everything together
#March 11, 2021. Elham Nourani. Radolfzell am Bodensee
#update: don't subsample hourly. just subset to closest minute for gannets and nazca boobies (the sub-min resolutions that cause the C stack error.)
#update: one minute is still causing errors. try rounding up to 15 min and then removing duplicated points
#update: nazca booby and magnificient frigatebird are still problematic, so reduce to 25 min.

#Qs: are the red footed boobies in Europa breeding?
# are all wandering albatrosses breeding?

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


# step 1: data prep ######
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

#common column names
cols <- Reduce(intersect, lapply(data_ls, colnames))
  
colonies <- data.frame(study.name = c("Foraging habitat of white-tailed tropicbirds (data from Santos et al. 2019)", "Frigatebirds breeding at Iguana Island, Panama", "Foraging ecology of masked boobies (Sula dactylatra) in the world’s largest “oceanic desert”",
                                      "Nazca booby Sula granti Isla Espanola, Galapagos.", "Galapagos Albatrosses", "Red footed boobies (Weimerskirch)", "Great frigatebirds (Weimerskirch)", "Red-tailed tropicbirds (Phaethon rubricauda) Round Island",
                                      "MPIAB PNIC hurricane frigate tracking"),
                       colony.name = c("Fernando de Noronha", "Iguana Island", "Motu Nui", "Isla Espanola", "Isla Espanola", "Genovesa Island", "Genovesa Island", "Round Island", "Isla Contoy"),
                       colony.lat = c(-3.86, 21.06, -27.2, -1.38, -1.38, 0.32, 0.32, -19.85, 21.50),
                       colony.long = c(-32.42, -73.30, -109.4, -89.67, -89.67, -89.96, -89.96, 57.79, -86.79))


#breeding periods for species tracked during both breeding and non-breeding season
breeding <- data.frame(study.name = c("Great frigatebirds (Weimerskirch)", "Red footed boobies (Weimerskirch)", "Galapagos Albatrosses", "Foraging habitat of white-tailed tropicbirds (data from Santos et al. 2019)"), #source for frigatebird, the other frigatebird study
                       breeding_start = c(3, 2, 3, 3),
                       breeding_end = c(8, 9, 12, 4))

#for nazca boobies and mag frigatebirds, sample randomly from the individuals
nb_to_keep <- data_ls$`data/From_Sophie/final_list_track_split/NaBo_Galap_split.Rdata.NaBo_Galap_split` %>% 
  distinct(individual.local.identifier) %>% 
  sample_n(55)
  
  
data_ls$`data/From_Sophie/final_list_track_split/NaBo_Galap_split.Rdata.NaBo_Galap_split` <- data_ls$`data/From_Sophie/final_list_track_split/NaBo_Galap_split.Rdata.NaBo_Galap_split` %>% 
  filter(individual.local.identifier %in% nb_to_keep$individual.local.identifier)
           

data_df <- data_ls %>% 
  map(dplyr::select, cols) %>% 
  reduce(rbind) %>% 
  rename(sci_name = individual.taxon.canonical.name) %>%  
  full_join(colonies, by = "study.name") %>% 
  mutate(sci_name = as.character(fct_recode(sci_name, 'Fregata magnificens' = "Fregata")),
         month = month(timestamp),
         year = year(timestamp),
         indID = as.character(paste(sci_name,individual.local.identifier, sep = "_"))) %>% 
  rowwise() %>% 
  mutate(TripID = paste(year, indID, as.character(TripID))) %>% 
  ungroup()

#filter for breeding season: remove entire tracks that fall within the non-breeding period!!

data_breeding <- data_df %>% 
  filter(study.name == "Great frigatebirds (Weimerskirch)" & between(month, 3,8) | 
           study.name == "Red footed boobies (Weimerskirch)" & between(month, 2,9) |
           study.name == "Galapagos Albatrosses" & between(month, 3,12)|
           study.name == "Foraging habitat of white-tailed tropicbirds (data from Santos et al. 2019)" & between(month, 3,4)|
           !(study.name %in% c("Great frigatebirds (Weimerskirch)", "Red footed boobies (Weimerskirch)", "Galapagos Albatrosses",
                               "Foraging habitat of white-tailed tropicbirds (data from Santos et al. 2019)")))
#summarize how many tracks and individuals per study

data_breeding %>% 
  group_by(study.name) %>% 
  summarize(n_tracks = n_distinct(TripID),
            n_ind = n_distinct(indID))


#make sure all track IDs are unique
data_breeding %>% 
  group_by(TripID) %>% 
  summarise(n = n_distinct(indID)) %>%
  filter(n > 1) #these are not unique. so, paste the ind name with the trip ID and year to make it unique
  

save(data_breeding, file = "R_files/mv_nbsample_w_colony.RData") #naza booby is a sample of 55 ind from the original 700 or whatever


# Peter Ryan data prep ######
#peter ryan data from Sophie. These are all breeding

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/Peter_Ryan_data_annotated_SplitTrip.Rdata") #PR_data_split
PR_data_split <- PR_data_split %>% 
  rename(colony.name = colony_name,
         colony.lat = lat_colony,
         colony.long = lon_colony)

save(PR_data_split, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/Peter_Ryan_data_annotated_SplitTrip.Rdata")

# Gremillet data prep ######
#all breeding (incubation and chick-rearing)
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
  mutate(timestamp = paste(DateGMT,TimeGMT, sep = " "),
                colony.name = ifelse(scientific_name == "Morus capensis", "Algoa Bay", "Ile Rouzic"),
                colony.lat = ifelse(scientific_name == "Morus capensis", -33.5,48.9),
                colony.long = ifelse(scientific_name == "Morus capensis", 25.5, -3.43)) %>% 
  ungroup() %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))

##summary data for Martina's cost of transport project (22.10.2021)
summary_data <- gannets %>%  group_by(year(timestamp), colony.name, common_name) %>%  distinct(BirdId) %>% 
  rename(bird_id = BirdId, colony_name = colony.name, year = `year(timestamp)`)

write.csv(summary_data, "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/David Gremillet/gannet_summary_info.csv") 

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


# step 2: put everything together ######

#open data #
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/waal_all.RData") #waal; called waal

RFB <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/RFBO_2012_ParamWind.csv", 
                stringsAsFactors = F) %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) #red foote boobie breeding season in Europa: June-Oct. this data is for October and November...

waal_RFB <- waal %>% 
  full_join(RFB) %>% 
  mutate(sci_name = ifelse(Species == "WAAL", "Diomedea exulans", "Sula sula"),
         colony.name = ifelse(Species == "WAAL", "Crozet Island", "Europa Island"),
         colony.lat = ifelse(Species == "WAAL", -46.4,-22.3),
         colony.long = ifelse(Species == "WAAL", 51.5, 40.3)) %>% 
  rename(location.long = Longitude,
         location.lat = Latitude,
         indID = BirdID,
         timestamp = date_time)

#sf <- st_as_sf(waal_RFB, coords = c("location.long", "location.lat"), crs = wgs)

#open other files from earlier in this script
all_files <- list("R_files/mv_nbsample_w_colony.RData",
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
  reduce(rbind) %>% 
  full_join(waal_RFB[,cols])

save(data_df_all, file = "R_files/all_spp_df_colony_waal.RData")

#get summary

data_df_all %>% 
  group_by(sci_name) %>% 
  summarize(n_tracks = n_distinct(TripID),
            n_ind = n_distinct(indID))

# step 3: add flyingsitting to all ######

#and only keep the flying points.

load("R_files/all_spp_df_colony_waal.RData") #data_df_all

#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(data_df_all$TripID),timestamps = data_df_all$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows

data_df_all <- data_df_all[-rows_to_delete,] 

#convert to move objects (i need to do this to calc speed for flyingsitting assignment)

move_ls <- lapply(split(data_df_all,data_df_all$sci_name),function(x){
  x <- x %>%
    arrange(TripID, timestamp) %>% 
    as.data.frame()
  mv <- move(x = x$location.long, y = x$location.lat, time = x$timestamp, data = x, animal = x$TripID, proj = wgs)
  mv
  
})

#time lag?

move_ls_tl <- lapply(move_ls, function(x){
  x$timeLag <-  unlist(lapply(timeLag(x, units="hours"),  c, NA))
  x
})

lapply(move_ls_tl, function(x) mean(x$timeLag, na.rm = T))


#add flying sitting
mycl <- makeCluster(10) 
clusterExport(mycl, c("move_ls")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(sf)
  library(sp)
  library(tidyverse)
  library(move)
})

(b <- Sys.time())

fl_sit_mv <- parLapply(mycl, move_ls, function(species){ #each species
  
  track_flsit <- lapply(split(species), function(track){
    
    #--STEP 1: drop points where the animal is not moving (i.e. sitting)
    if("FlyingSitting" %in% colnames(track@idData)){ #if there is data fro flyingsitting, it ends up in the data, if not, it will be in iData
      
      track$speed_kmh <- c(NA, speed(track) * 3.6)
      track$FlyingSitting <- ifelse(track$speed_kmh > 2, "flying", "sitting")
      
      #take flyingsitting out of the iData
      track@idData <- track@idData[colnames(track@idData) != "FlyingSitting"]
      # track
    } else {
      track$speed_kmh <- c(NA, speed(track) * 3.6)
    }
    
    track$colony.name <- track@idData$colony.name
    track$colony.lat <- track@idData$colony.lat
    track$colony.long <- track@idData$colony.long
    track$inID <- track@idData$indID
    track$TripID <- track@idData$TripID
    track$sci_name <- track@idData$sci_name
    
    track
    
  }) %>% 
    reduce(rbind)
  track_flsit
  
})

Sys.time() - b #1.7 min

stopCluster(mycl) 


save(fl_sit_mv, file = "R_files/all_spp_sitting_flying.RData")

# step 3: prep for annotation ######

all_spp <- fl_sit_mv %>%
  reduce(rbind) %>% 
  as.data.frame() %>% 
  dplyr::select(-c("coords.x1", "coords.x2")) %>% 
  mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) %>% 
  as.data.frame()

save(all_spp, file = "R_files/all_spp_sitting_flying_df.RData")

colnames(all_spp)[c(2,3)] <- c("location-long","location-lat")

spp1 <- all_spp[1:(nrow(all_spp)/2),]
spp2 <- all_spp[((nrow(all_spp)/2) + 1): nrow(all_spp),]


write.csv(spp1, "R_files/all_spp_annotation_1.csv")
write.csv(spp2, "R_files/all_spp_annotation_2.csv")

#extract colony summaries
all_spp %>% 
  group_by(sci_name) %>% 
  summarize(n_c = n_distinct(colony.name),
            name = head(colony.name,1))

all_spp %>% 
  filter(sci_name == "Sula sula") %>% 
  group_by(colony.name) %>% 
  summarise(n())
  
all_spp %>% 
  filter(sci_name == "Fregata magnificens") %>% 
  group_by(colony.name) %>% 
  summarise(n())
  
#extract colony IDs for each track/species

trip_col <- all_spp %>% 
  group_by(sci_name, TripID) %>% 
  summarize(colony.name = head(colony.name,1),
            colony.lat = head(colony.lat,1),
            colony.long = head(colony.long,1))

save(trip_col, file = "R_files/trip_colony.RData")
