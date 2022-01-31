#script for generating random steps and running ssf for seabird data
#follows up from all_species_prep_2021.R
#Elham Nourani. Radofzell, DE. April 8. 2021
#clogit z transofmation didn't change things.
#update Aug. 19: separate the code into step generation and analysis. The final manuscript will include step generation, then permutations (null_modeling_seabirds.R)

library(tidyverse)
library(lubridate)
library(sf)
library(parallel)
library(move)
library(corrr)
library(INLA)
library(survival)
library(TwoStepCLogit)
library(CircStats)
library(circular)
library(fitdistrplus)


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
source("/home/mahle68/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")

species <- read.csv("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/data/final_datasets.csv")

# ----------- Step 1: generate alternative steps ####

#open the data. already has flyingsitting assignments
load("R_files/all_spp_sitting_flying_df.RData") #all_spp

all_spp$timestamp <- as.POSIXct(strptime(all_spp$timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")
all_spp$year <- year(all_spp$timestamp)

#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(all_spp$TripID),timestamps = all_spp$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows

data_df_all <- all_spp[-rows_to_delete,] 


#summary info

all_spp %>% 
  group_by(sci_name, colony.name) %>% 
  summarize(n_tracks = n_distinct(TripID),
            n_ind = n_distinct(inID),
            n_yr = n_distinct(year),
            min_yr = min(year),
            max_yr = max(year))


#convert to move objects (i need to do this to calc speed for flyingsitting assignment)

move_ls <- lapply(split(all_spp, paste(all_spp$sci_name, all_spp$colony.name, sep = "_")),function(x){
  x <- x %>%
    arrange(TripID, timestamp) %>% 
    as.data.frame()
  mv <- move(x = x$location.long, y = x$location.lat, time = x$timestamp, data = x, animal = x$TripID, proj = wgs)
  mv
  
})

save(move_ls, file = "20_spp_col_move_ls.RData")


#### summary info

lapply(move_ls, function(x) length(split(x)))

#### the three species that are problematic, c("Sula granti","Fregata magnificens", "Diomedea exulans"), sub-sample manually to hourly. also reduce the number of individuals.

load("R_files/all_spp_sitting_flying_df.RData") #all_spp
#make sure frigatebird is only in march and april: breeding season

 lapply(move_ls, function(x){
  unlist(lapply(timeLag(x, units="hours"),  c, NA))
  x
})

 
#for Fregata magnificens_Isla_Contoy, remove tracks with very low resolution (ie. > 3 hr)
freg_contoy <- move_ls[[5]]
#extract average time lag at each track
lapply(timeLag(freg_contoy, units = "hours"), mean, na.rm = T)

freg_contoy$timelag <- unlist(lapply(timeLag(freg_contoy, units = "hours"),c, NA)) #all are over 3 hours. don't use

unlist(lapply(timeLag(move_ls[[12]], units = "hours"),c, NA)) 

#for 
 
# spp3 <- all_spp %>% 
#   mutate(timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
#   filter(sci_name == "Fregata magnificens" & between(month(timestamp), 3,4)|
#            sci_name %in% c("Diomedea exulans", "Sula granti")) %>% 
#   group_by(sci_name, TripID, as.Date(timestamp), hour(timestamp)) %>% 
#   slice(1) %>% 
#   ungroup()
# 
# #remove dublicates
# rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(spp3$TripID),timestamps = spp3$timestamp,
#                                                         sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows
# 
# spp3 <- spp3[-rows_to_delete,] 
# 
# #convert to move objects (i need to do this to calc speed for flyingsitting assignment)
# 
# move_ls <- lapply(split(spp3,spp3$sci_name),function(x){
#   x <- x %>%
#     arrange(TripID, timestamp) %>% 
#     as.data.frame()
#   mv <- move(x = x$location.long, y = x$location.lat, time = x$timestamp, data = x, animal = x$TripID, proj = wgs)
#   mv
#   
# })
# 
# save(move_ls, file = "3_spp_move_ls.RData")

load("20_spp_col_move_ls.RData") #move_ls

int <- 60 #in minutes. for galapagos albatross, make this 90
tol <- 15 #in minutes
n_alt <- 30 #n of alternative points... 


mycl <- makeCluster(10) 
clusterExport(mycl, c("move_ls", "wgs", "meters_proj", "NCEP.loxodrome.na", "int", "tol", "n_alt")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(sf)
  library(sp)
  library(tidyverse)
  library(move)
  library(CircStats)
  library(circular)
  library(fitdistrplus)
})

(b <- Sys.time())

#used_av_ls_60_15 <- parLapply(mycl, move_ls, function(species){ #each species
  
parLapply(mycl, move_ls[13:20], function(species){ #each species
  
#lapply(move_ls[c(13:20)], function(species){ #each species
  sp_obj_ls <- lapply(split(species), function(track){
    
    #--STEP 1: drop points where the animal is not moving (i.e. sitting)
    track_flying <- track[is.na(track$FlyingSitting) | track$FlyingSitting == "flying"]
    
    
    track_th <- track_flying %>%
      thinTrackTime(interval = as.difftime(int, units='mins'),
                    tolerance = as.difftime(tol, units='mins')) #the unselected bursts are the large gaps between the selected ones
    #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst... longer gaps will divide the brusts) 
    track_th$selected <- c(as.character(track_th@burstId),NA) #assign selected as a variable
    track_th$burst_id <-c(1,rep(NA,nrow(track_th)-1)) #define value for first row
    
    if(nrow(track_th@data) == 1){
      track_th@data$burst_id <- track_th$burst_id
    } else {for(i in 2:nrow(track_th@data)){
      
      if(i== nrow(track_th@data)){
        track_th@data$burst_id[i] <- NA
      } else
        if(track_th@data[i-1,"selected"] == "selected"){
          track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"]
        } else {
          track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"] + 1
        }
    }
    }
    #convert back to a move object (from move burst)
    track_th <- as(track_th,"Move")
    
    #--STEP 3: calculate step lengths and turning angles 
    #sl_ and ta_ calculations should be done for each burst. converting to a move burst doesnt make this automatic. so just split manually
    burst_ls <- split(track_th, track_th$burst_id)
    burst_ls <- Filter(function(x) length(x) >= 3, burst_ls) #remove bursts with less than 3 observations
    
    burst_ls <- lapply(burst_ls, function(burst){
      burst$step_length <- c(distance(burst),NA) 
      burst$turning_angle <- c(NA,turnAngleGc(burst),NA)
      burst
    })
    
    #put burst_ls into one sp dataframe
    bursted_sp <- do.call(rbind, burst_ls)
    
    #reassign values
    if(length(bursted_sp) >= 1){
      bursted_sp$sci_name <- track@idData$sci_name
      bursted_sp$indID <- track@idData$indID
      bursted_sp$TripID <- track@idData$TripID
      bursted_sp$colony.name <- track@idData$colony.name
    }
    
    #bursted_sp$track<-track@idData$seg_id 
    
    bursted_sp
    
  }) %>% 
    Filter(function(x) length(x) > 1, .) #remove segments with no observation
  
  #--STEP 4: estimate step length and turning angle distributions
  #put everything in one df
  bursted_df <- sp_obj_ls %>%  
    reduce(rbind) %>% 
    as.data.frame() %>% 
    dplyr::select(-c("coords.x1","coords.x2"))
  
  #estimate von Mises parameters for turning angles
  #calculate the averages (mu).steps: 1) convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan.OR use circular::mean.circular
  mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
  sl <- bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  
  #plot
  pdf(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/2021_all_spp/1hr_15min_",species@idData$sci_name[1], "_", species@idData$colony.name[1], ".pdf"))
  
  par(mfrow=c(1,2))
  hist(sl,freq=F,main="",xlab = "Step length (km)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  
  hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
  dev.off()
  
  #--STEP 5: produce alternative steps
  used_av_track <- lapply(sp_obj_ls, function(track){ #for each trackment
    
    lapply(split(track,track$burst_id),function(burst){ #for each burst,
      
      #assign unique step id
      burst$step_id <- 1:nrow(burst)
      
      lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
        
        current_point<- burst[this_point,]
        previous_point <- burst[this_point-1,] #this is the previous point, for calculating turning angle.
        used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
        
        #calculate bearing of previous point
        #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
        prev_bearing <- NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
                                          previous_point@coords[,1], current_point@coords[,1])
        
        current_point_m <- spTransform(current_point, meters_proj) #convert to meters proj
        
        #randomly generate n alternative points
        rnd <- data.frame(turning_angle = as.vector(rvonmises(n = n_alt, mu = mu, kappa = kappa)), #randomly generate n step lengths and turning angles
                          step_length = rgamma(n = n_alt, shape = fit.gamma1$estimate[[1]], rate = fit.gamma1$estimate[[2]]) * 1000) %>% 
          #find the gepgraphic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
          mutate(lon = current_point_m@coords[,1] + step_length*cos(turning_angle),
                 lat = current_point_m@coords[,2] + step_length*sin(turning_angle))
        
        
        #covnert back to lat-lon proj
        rnd_sp <- rnd
        coordinates(rnd_sp) <- ~lon+lat
        proj4string(rnd_sp) <- meters_proj
        rnd_sp <- spTransform(rnd_sp, wgs)
        
        #put used and available points together
        df <- used_point@data %>%  
          slice(rep(row_number(), n_alt+1)) %>% #paste each row n_alt times for the used and alternative steps
          mutate(location.long = c(head(location.long,1),rnd_sp@coords[,1]), #the coordinates were called x and y in the previous version
                 location.lat = c(head(location.lat,1),rnd_sp@coords[,2]),
                 turning_angle = c(head(turning_angle,1),rnd_sp$turning_angle),
                 step_length = c(head(step_length,1),rnd_sp$step_length),
                 used = c(1,rep(0,n_alt)))  %>%
          rowwise() %>% 
          mutate(heading = NCEP.loxodrome.na(lat1 = current_point$location.lat, lat2 = location.lat, lon1 = current_point$location.long, lon2 = location.long)) %>% 
          as.data.frame()
        
        df
        
      }) %>% 
        reduce(rbind)
      
    }) %>% 
      reduce(rbind)
    
  }) %>% 
    reduce(rbind)
  #used_av_track
  
  #save the file
      save(used_av_track, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/ssf_input/1hr_15min/",species@idData$sci_name[1], "_", species@idData$colony.name[1], ".RData"))
  
  })

Sys.time() - b 

stopCluster(mycl) 



# ----------- Step 2: annotate#####

files <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/ssf_input/1hr_15min/", full.names = T)

used_av_ls_60_15 <- sapply(files, function(x) mget(load(x)), simplify = TRUE)

save(used_av_ls_60_15, file = "R_files/ssf_input_all_60_15.RData")

load("R_files/ssf_input_all_60_30_3spp.RData") #used_av_ls_60_15

#create one dataframe with movebank specs
used_av_df_60_15 <- lapply(c(1:length(used_av_ls_60_15)), function(i){
  
  data <- used_av_ls_60_15[[i]] %>% 
    mutate(date_time = timestamp,
           timestamp = paste(as.character(timestamp),"000",sep = "."),
           stratum = paste(TripID, burst_id, step_id, sep = "_"),
           year = year(timestamp)) %>% 
    
    as.data.frame()
}) %>% 
  reduce(rbind)

#rename columns
colnames(used_av_df_60_15)[c(2,3)] <- c("location-long","location-lat")

save(used_av_df_60_15, file = "R_files/ssf_input_df_60_15_50.RData")


#create multiple files

n_chunks <- ceiling(nrow(used_av_df_60_15)/1e6)
nrow_chunks <- round(nrow(used_av_df_60_15)/n_chunks)

r <- rep(1: n_chunks,  each = nrow_chunks)
r <- r[-length(r)] #remove one value to make sure the length of r is the same as nrow(used_av_df_60_50)

chunks <- split(used_av_df_60_15, r)


#save as csv
lapply(c(1:length(chunks)), function(i){
  
  write.csv(chunks[[i]], paste0("R_files/ssf_input_df_60_15_50_", i, "_.csv"))
  
})


# # summary stats
# used_av_df_60_15 %>% 
#   #group_by(species) %>% 
#   group_by(group) %>% 
#   summarise(yrs_min = min(year(date_time)),
#             yrs_max = max(year(date_time)),
#             n_ind = n_distinct(ind),
#             n_tracks = n_distinct(track))

#visual inspection
data_sf <- used_av_df_60_15 %>% 
  filter(used == 1) %>% 
  st_as_sf(coords = c(2,3), crs = wgs)

mapview(data_sf, zcol = "species")


#put all annotated data together
#load 18 species
load("R_files/ssf_input_df_60_15_50.RData") #used_av_df_60_15 #i don't know what this is. probably from an earlier version.

used_av_df_60_15 <- used_av_df_60_15 %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  dplyr::select(-"date_time")

#rename lat and long
colnames(used_av_df_60_15)[c(2,3)] <- c("location.long", "location.lat")
  

#load annotated files
files_ls <- list.files("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/annotation/10_15_all_colony_sp/", pattern = ".csv",recursive = T, full.names = T)

#append all species
ann_18 <- lapply(files_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  drop_na(ECMWF.ERA5.SL.Sea.Surface.Temperature) %>% #remove points over land
  dplyr::select(-"date_time") %>% 
  mutate(stratum = paste(TripID, burst_id, step_id, sep = "_")) %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.,
         v10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component.,
         wh = ECMWF.ERA5.SL.Significant.Wave.Height) %>%
  mutate(row_id = row_number(),
         delta_t = sst - t2m,
         wind_support= wind_support(u = u10m,v = v10m, heading = heading),
         cross_wind= cross_wind(u = u10m, v = v10m,heading = heading),
         wind_speed = sqrt(u10m^2 + v10m^2),
         abs_cross_wind = abs(cross_wind(u = u10m, v = v10m, heading = heading)),
         year = year(timestamp)) %>%
  left_join(species[,c(1:2)], by = c("sci_name" = "scientific.name")) %>% 
  rename(common_name = species) #%>%
  #full_join(used_av_df_60_15) #what is this!?


save(ann_18, file = "R_files/ann_18.RData")

#extract startum IDs for those that have less than 50 alternative points over the sea
less_than_30 <- ann_18 %>% 
  filter(used == 0) %>% 
  group_by(stratum) %>% 
  summarise(n = n()) %>% 
  filter(n < 30) 

# retain 50 alternative steps per stratum
used <- ann_18 %>%
  filter(!(stratum %in% less_than_30$stratum)) %>%
  filter(used == 1)

used_avail_30 <- ann_18 %>%
  filter(!(stratum %in% less_than_30$stratum)) %>%
  filter(used == 0) %>%
  group_by(stratum) %>%
  sample_n(30, replace = F) %>%
  ungroup() %>%
  full_join(used) #append the used levels

#make sure all strata have 51 points
no_used <- used_avail_30 %>%
  group_by(stratum) %>%
  summarise(n = n()) %>%
  filter(n < 31)

ann_30 <- used_avail_30 %>%
  filter(!(stratum %in% no_used$stratum))
  
#do we have 51 rows per stratum?
ann_30 %>% 
  group_by(stratum) %>% 
  summarise(n = n()) %>% 
  filter(n != 31) 


save(ann_30, file = "R_files/ssf_input_annotated_60_15_30alt_18spp.RData")


# ----------- Step 3: box plots ####

load("R_files/ssf_input_annotated_60_15_30alt_18spp.RData") #ann_30

ann_30 <- ann_30 %>%  # a lot of NAs in common names
  mutate(group = paste(common_name, colony.name, sep = "_")) %>% 
  as.data.frame()

#split the data into dynamic soarers and non dynamic soarers
soarers <- ann_30[ann_30$sci_name %in% species[species$flight.type  == "dynamic soaring", "scientific.name"],]
non_soarers <- ann_30[ann_30$sci_name %in% species[species$flight.type  != "dynamic soaring", "scientific.name"],]

#pdf("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/2021_boxplots/15spp_soarers_40_1hr.pdf", width = 15, height = 9)

#separate boxplots for colonies. Maybe put each colony in one figure.

labels <- c("Wind support", "cross wind","Wind speed")
variables <- c("wind_support", "cross_wind", "wind_speed")

X11(width = 15, height = 9)
par(mfrow= c(3,1), 
    oma = c(2,0,3,0), 
    las = 1)


for(i in 1:length(variables)){
  
  boxplot(soarers[,variables[i]] ~ soarers[,"common_name"], data = soarers, boxfill = NA, border = NA, main = labels[i], xlab = "", ylab = "",xaxt = "n")
  
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(soarers[soarers$used == 1, variables[i]] ~ soarers[soarers$used == 1,"common_name"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = "orange", lwd = 0.7, outpch = 20, outcex = 0.8,
          boxwex = 0.25, at = 1:length(unique(soarers[soarers$used == 1, "common_name"])) - 0.15)
  
  boxplot(soarers[soarers$used == 0, variables[i]] ~ soarers[soarers$used == 0, "common_name"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(soarers[soarers$used == 1 , "common_name"])) + 0.15)
 
  axis(side = 1, at = 1:length(unique(soarers$common_name)), line = 0, labels = FALSE, 
       tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015)
  
   text(x = 1:length(unique(soarers$common_name))+0.18,
       y = par("usr")[3],
       labels = unique(soarers$common_name), xpd = NA, srt = 35, cex = 1, adj = 1.18)
  
}

mtext("Instantaneous values at each step 1hr", side = 3, outer = T, cex = 1.3)

dev.off()

pdf("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/2021_boxplots/15spp_non_soarers_40_1hr.pdf", width = 15, height = 9)

X11(width = 7, height = 9)
par(mfrow= c(3,1), 
    oma = c(2,0,3,0), 
    las = 1)

for(i in 1:length(variables)){
  
  boxplot(non_soarers[,variables[i]] ~ non_soarers[,"common_name"], data = non_soarers, boxfill = NA, border = NA, main = labels[i], xlab = "", ylab = "",xaxt = "n")
  
  if(i == 1){
    legend("topleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(non_soarers[non_soarers$used == 1, variables[i]] ~ non_soarers[non_soarers$used == 1,"common_name"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(non_soarers[non_soarers$used == 1, "common_name"])) - 0.15)
  
  boxplot(non_soarers[non_soarers$used == 0, variables[i]] ~ non_soarers[non_soarers$used == 0, "common_name"], 
          yaxt = "n", xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(non_soarers[non_soarers$used == 1 , "common_name"])) + 0.15)
  
  axis(side = 1, at = 1:length(unique(non_soarers$common_name)), line = 0, labels = FALSE, 
       tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015)
  
  text(x = 1:length(unique(non_soarers$common_name))+0.18,
       y = par("usr")[3],
       labels = unique(non_soarers$common_name), xpd = NA, srt = 35, cex = 1, adj = 1.2)
  
}

mtext("Instantaneous values at each step 1hr", side = 3, outer = T, cex = 1.3)

dev.off()
  
# ----------- Step 4: clogit ####
load("R_files/ssf_input_annotated_60_30_40alt_18spp.RData") #ann_40

f1 <- formula(used ~  wind_speed + wind_support + common_name +
                           strata(stratum))
f2 <- formula(used ~  wind_speed * common_name + wind_support * common_name +
                strata(stratum))
m1 <- clogit(f1, data = ann_40)
m2 <- clogit(f2, data = ann_40)

#interpret the coefficients
#https://it.unt.edu/interpreting-glm-coefficients

wspd_coeffs <- data.frame(coef = summary(m2)$coefficients[,1],
                                  se = summary(m2)$coefficients[,3]) %>% 
  filter(str_detect(row.names(.), "wind_speed")) %>% 
  mutate(species = c(levels(as.factor(ann_40$common_name))[1], row.names(.)[2:nrow(.)])) %>% 
  mutate(across("species", str_replace, "wind_speed:common_name", ""))
  
wspd_coeffs$logit <- c(wspd_coeffs[1, "coef"],  wspd_coeffs[2:18, "coef"] + wspd_coeffs[1, "coef"])

wspd_coeffs <- wspd_coeffs %>% 
  rowwise() %>% 
  mutate(upper = logit + se,
         lower = logit - se)

#plot
X11(width = 9, height = 5.5)
par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(3, 6.7, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(min(wspd_coeffs$lower)-0.05,max(wspd_coeffs$upper) + -0.05), 
     ylim = c(1,19), xlab = "wind speed logit", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(wspd_coeffs$logit, factor(wspd_coeffs$species), col = "steelblue1", pch = 20, cex = 1.3)
arrows(wspd_coeffs$lower, c(1:n_distinct(wspd_coeffs$species)),
       wspd_coeffs$upper, c(1:n_distinct(wspd_coeffs$species)),
       col = "steelblue1", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line
#add axes
axis(side= 1, at= c(-4,-2,0,2), labels= c("-4","-2", "0", "2"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:n_distinct(wspd_coeffs$species)), #line=-4.8, 
     labels = levels(factor(wspd_coeffs$species)),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 

#### one model per species (consider Two step clogit to control for year and individual)

f3 <- formula(used ~  wind_speed +
                strata(stratum))

ms_wspd <- lapply(split(ann_40, ann_40$common_name),function(x){
  clogit(f3 , data = x)
}
)


wspd <- lapply(ms_wspd, function(x){
  data.frame(coef = summary(x)$coefficients[,1],
             se = summary(x)$coefficients[,3],
             p = summary(x)$coefficients[,5])
}) %>% 
  reduce(rbind) %>% 
  mutate(species = names(ms_wspd),
         lower = coef - se,
         upper = coef + se) %>% 
  mutate(sig = ifelse(p <= 0.05, "yes", "no"))


#plot
X11(width = 9, height = 5.5)
par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(3, 6.7, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(min(wspd$lower)-0.05,max(wspd$upper) + -0.05), 
     ylim = c(1,19), xlab = "wind speed logit", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(wspd$coef, factor(wspd$species), col = "steelblue1", pch = 20, cex = 1.3)
arrows(wspd$lower, c(1:n_distinct(wspd$species)),
       wspd$upper, c(1:n_distinct(wspd$species)),
       col = "steelblue1", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line
#add axes
axis(side= 1, at= c(-4,-2,0,2), labels= c("-4","-2", "0", "2"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:n_distinct(wspd$species)), #line=-4.8, 
     labels = levels(factor(wspd$species)),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 

#wind support
f4 <- formula(used ~  wind_support +
                strata(stratum))

ms_wspt <- lapply(split(ann_40, ann_40$common_name),function(x){
  clogit(f4 , data = x)
}
)


wspt <- lapply(ms_wspt, function(x){
  data.frame(coef = summary(x)$coefficients[,1],
             se = summary(x)$coefficients[,3],
             p = summary(x)$coefficients[,5])
}) %>% 
  reduce(rbind) %>% 
  mutate(species = names(ms_wspt),
         lower = coef - se,
         upper = coef + se) %>% 
  mutate(sig = ifelse(p <= 0.05, "yes", "no"))


#plot
X11(width = 9, height = 5.5)
par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(3, 6.7, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(min(wspt$lower)-0.05,max(wspt$upper) + 0.05), 
     ylim = c(1,19), xlab = "wind speed logit", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(wspt$coef, factor(wspt$species), col = "orange", pch = 20, cex = 1.3)
arrows(wspt$lower, c(1:n_distinct(wspt$species)),
       wspt$upper, c(1:n_distinct(wspt$species)),
       col = "orange", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line
#add axes
axis(side= 1, at= c(-0.1,0,0.1), labels= c("-0.1", "0", "0.1"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:n_distinct(wspt$species)), #line=-4.8, 
     labels = levels(factor(wspt$species)),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 


# ----------- Step 5: two step clogit ####

load("R_files/ssf_input_annotated_60_30_40alt_18spp.RData") #ann_30
load("R_files/input_for_smooth_inla.RData") #ann_40 with z-transformation

#empty data frame for coeffs and standard errors
coefs <- data.frame(species = names(split(ann_40, ann_40$common_name)),
                    wind_support_coef = NA,
                    wind_support_se = NA,
                    wind_speed_coef = NA,
                    wind_speed_se = NA)

for(i in 1:length(split(ann_40, ann_40$common_name))){
  
  x <- split(ann_40, ann_40$common_name)[[i]]
    
  if(length(unique(x$year)) > 1){
    m <- Ts.estim(used ~ wind_support_z + wind_speed_z +  strata(stratum) + cluster(year),
             random = ~ wind_support_z + wind_speed_z, 
             data = x, 
             D="UN(1)") #D is the structure of the output matrix
    
    coefs[i,"wind_support_coef"] <- m$beta[1]
    coefs[i,"wind_speed_coef"] <- m$beta[2]
    
    coefs[i,"wind_support_se"] <- m$se[1]
    coefs[i,"wind_speed_se"] <- m$se[2]
      
  } else {
    m <- clogit(used ~ wind_support_z + wind_speed_z +  strata(stratum), data = x)
    
    coefs[i,"wind_support_coef"] <- m$coefficients[1]
    coefs[i,"wind_speed_coef"] <- m$coefficients[2]
    
    coefs[i,"wind_support_se"] <- summary(m)$coefficients[1,3]
    coefs[i,"wind_speed_se"] <- summary(m)$coefficients[2,3]
  }
  
 # models[i] <- m
}


coefs <- coefs %>% 
  rowwise() %>% 
  mutate(upper_wsd = wind_speed_coef + wind_speed_se,
         lower_wsd = wind_speed_coef - wind_speed_se,
         upper_wst = wind_support_coef + wind_support_se,
         lower_wst = wind_support_coef - wind_support_se)

save(coefs, file = "R_files/two_step_coefs.RData")
#save(models, file = "R_files/two_step_models.RData")

#get rid of white-tailed tropicbird... the range of values is bizarre
coefs <- coefs[-18,]


#plot
X11(width = 9, height = 5.5)
par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(3, 6.7, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(min(coefs$lower_wsd)-0.05,max(coefs$upper_wsd) + -0.05), 
     ylim = c(1,19), xlab = "wind speed coef", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(coefs$wind_speed_coef, factor(coefs$species), col = "steelblue1", pch = 20, cex = 1.3)
arrows(coefs$lower_wsd, c(1:n_distinct(coefs$species)),
       coefs$upper_wsd, c(1:n_distinct(coefs$species)),
       col = "steelblue1", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line
#add axes
axis(side= 1, at= c(-4,-2,0,2, 4), labels= c("-4","-2", "0", "2", "4"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:n_distinct(coefs$species)), #line=-4.8, 
     labels = levels(factor(coefs$species)),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 



jpeg("/home/enourani/ownCloud/Work/conferences/seabirds_2021/windspt_coefs.jpg", width = 6, height = 3.5, units = "in", res = 500)

#plot wind support
X11(width = 6, height = 3.5)
par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(4.5, 6.7, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-1,1), 
     ylim = c(1,18), xlab = "Wind support coefficients", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(coefs$wind_support_coef, factor(coefs$species), col = "steelblue1", pch = 20, cex = 2.5)
arrows(coefs$lower_wst, c(1:n_distinct(coefs$species)),
       coefs$upper_wst, c(1:n_distinct(coefs$species)),
       col = "steelblue1", code = 3, length = 0.03, angle = 90, lwd = 2) #angle of 90 to make the arrow head as straight as a line
#add axes
axis(side= 1, at= c(-0.5,0,0.5), labels= c("-0.5", "0", "0.5"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:n_distinct(coefs$species)), #line=-4.8, 
     labels = levels(factor(coefs$species)),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 

dev.off()

#all species at once
M1 <- Ts.estim(used ~ wind_support_z + wind_speed_z +  strata(stratum) + cluster(common_name),
         random = ~ wind_support_z + wind_speed_z, 
         data = ann_40, 
         D="UN(1)") #D is the structure of the output matrix

save(M1, file = "R_files/two_step_one_model.RData")

# ----------- Step 6: INLA ####
load("R_files/ssf_input_annotated_60_30_40alt_18spp.RData") #ann_40
#all flying, all points over the sea, check later that all are adults.

#correlation
ann_40 %>% 
  dplyr::select(c("location.lat", "wind_speed", "wind_support", "cross_wind", "abs_cross_wind")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: none...wind speed and abs_cross_wind have a corr of 0.598

# #z-transform
# all_data <- ann_50_all %>% 
#   #group_by(species) 
#   mutate_at(c("delta_t", "wind_speed", "wind_support", "wind_support_var", "abs_cross_wind", "delta_t_var"),
#             #list(z = ~scale(.)))
#             list(z = ~as.numeric(scale(.)))) %>%
#   arrange(stratum, desc(used)) %>% 
#   group_by(stratum) %>%  
#   mutate(lat_at_used = head(zone,1)) %>%  #add a variable for latitudinal zone. This will assign the lat zone of the used point to the entire stratum
#   ungroup() 

#model
mean.beta <- 0
prec.beta <- 1e-4 

formulaM1 <- used ~ -1 + wind_speed_z + wind_support_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(year1, wind_speed, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(year2, wind_support,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


(b <- Sys.time())

lapply(split(ann_40, ann_40$common_name), function(x){
  
  x <- x %>% 
    mutate(year1 = factor(year),
           year2 = factor(year),
           year3 = factor(year),
           # ind1 = factor(indID),
           # ind2 = factor(indID),
           # ind3 = factor(indID),
           stratum = factor(stratum))
  
M1 <- inla(formula = formulaM1, family ="Poisson",  
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           #control.inla = list(force.diagonal = T),
           data = x,
           num.threads = 10,
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T, cpo = T))

})

Sys.time() - b #2.3 min


## all species at once

ann_40 <- ann_40 %>% 
  mutate(species1 = factor(common_name),
         species2 = factor(common_name),
         species3 = factor(common_name),
         stratum = factor(stratum)) %>% 
  mutate_at(c("wind_speed", "wind_support"),
            list(z = ~as.numeric(scale(.)))) %>% 
  mutate_at(c("wind_speed","wind_support"),
            list(group = ~inla.group(.,n = 50, method = "cut"))) %>%
  as.data.frame()

save(ann_40, file = "R_files/input_for_smooth_inla.RData")

formulaM2 <- used ~ -1 + wind_speed + wind_support +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
f(species1, wind_speed, model = "iid", 
  hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
f(species2, wind_support,  model = "iid",
  hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

formulaM3 <- used ~ -1 + species1 +
  f(wind_support_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T)))


M2 <- inla(formula = formulaM2, family ="Poisson",  
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           control.inla = list(force.diagonal = T),
           data = ann_40,
           num.threads = 10,
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T, cpo = T))

(b <- Sys.time())
M3 <- inla(formulaM3, family ="Poisson", 
          control.fixed = list(
            mean = mean.beta,
            prec = list(default = prec.beta)),
          data = ann_40,
          num.threads = 10,
          control.predictor = list(compute = T), #list(link = 1), #link is only relevant for NA observations. required to set the right link (i.e., the logit function) 
          #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
          control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))
Sys.time() - b 