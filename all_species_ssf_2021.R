#script for generating random steps and running ssf for seabird data
#follows up from all_species_prep_2021.R
#Elham Nourani. Radofzell, DE. April 8. 2021

library(tidyverse)
library(lubridate)
library(sf)
library(parallel)
library(move)


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



# ------------ step 1 : generate alternative steps ####

#open the data. this data is already subsampled to hourly
load("R_files/all_spp_df_25min.RData") #data_df_all

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

save(move_ls, file = "17_spp_move_ls.RData")
save(move_ls, file = "17_spp_move_ls_25.RData")


#### summary info

lapply(move_ls, function(x) length(split(x)))

####
load("17_spp_move_ls_25.RData")

#exclude nazca booby and mag frigatebird
move_no_nb_mf <- move_ls[!(names(move_ls) %in% c("Sula granti","Fregata magnificens"))] 


mycl <- makeCluster(10) 
clusterExport(mycl, c("move_ls", "wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the function

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

#used_av_ls_60_30 <- parLapply(mycl, move_ls, function(species){ #each species
  
parLapply(mycl, move_ls, function(species){ #each species
  
  sp_obj_ls <- lapply(split(species), function(track){
    
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
    
    track_flying <- track[is.na(track$FlyingSitting) | track$FlyingSitting == "flying"]
    
    
    track_th <- track_flying %>%
      thinTrackTime(interval = as.difftime(60, units='mins'),
                    tolerance = as.difftime(30, units='mins')) #the unselected bursts are the large gaps between the selected ones
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
    
    #put burst_ls into one dataframe
    bursted_sp <- do.call(rbind, burst_ls)
    
    #reassign values
    
    if(length(bursted_sp) >= 1){
      bursted_sp$sci_name<-track@idData$sci_name
      bursted_sp$indID<-track@idData$indID
      bursted_sp$TripID<-track@idData$TripID
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
  pdf(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/2021_all_spp/",species@idData$sci_name[1], ".pdf"))
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
        
        #randomly generate 50 step lengths and turning angles
        rta <- as.vector(rvonmises(n = 50, mu = mu, kappa = kappa)) #generate random turning angles with von mises distribution (in radians)
        rsl <- rgamma(n = 50, shape=fit.gamma1$estimate[[1]], rate = fit.gamma1$estimate[[2]]) * 1000  #generate random step lengths from the gamma distribution. make sure unit is meters
        
        #calculate bearing of previous point
        #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
        prev_bearing <- NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
                                          previous_point@coords[,1], current_point@coords[,1])
        
        #find the gepgraphic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
        current_point_m <- spTransform(current_point, meters_proj) #convert to meters proj
        rnd <- data.frame(lon = current_point_m@coords[,1] + rsl*cos(rta),lat = current_point_m@coords[,2] + rsl*sin(rta)) #for this to work, lat and lon should be in meters as well. boo. coordinates in meters?
        
        #covnert back to lat-lon proj
        rnd_sp <- rnd
        coordinates(rnd_sp) <- ~lon+lat
        proj4string(rnd_sp) <- meters_proj
        rnd_sp <- spTransform(rnd_sp, wgs)
        
        #put used and available points together
        df <- used_point@data %>%  
          slice(rep(row_number(),51)) %>% #paste each row 50 times for the used and alternative steps
          mutate(x = c(head(location.long,1),rnd_sp@coords[,1]),
                 y = c(head(location.lat,1),rnd_sp@coords[,2]),
                 used = c(1,rep(0,50)))  %>%
          rowwise() %>% 
          mutate(heading = NCEP.loxodrome.na(lat1 = current_point$location.lat, lat2 = y, lon1 = current_point$location.long, lon2 = x)) %>% 
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
  save(used_av_track, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/ssf_input/",species@idData$sci_name[1], ".RData"))
  
})

Sys.time() - b 

stopCluster(mycl) 



# ---------- STEP 4: annotate#####

files <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/ssf_input/", full.names = T)

used_av_ls_60_30 <- sapply(files, function(x) mget(load(x)), simplify = TRUE)

load("2021/ssf_input_all_60_30_150_updated.RData") #used_av_ls_60_30

#create one dataframe with movebank specs
used_av_df_60_30 <- lapply(c(1:length(used_av_ls_60_30)), function(i){
  
  data <- used_av_ls_60_30[[i]] %>% 
    mutate(date_time = timestamp,
           timestamp = paste(as.character(timestamp),"000",sep = ".")) %>% 
    
    as.data.frame()
}) %>% 
  reduce(rbind)

save(used_av_ls_60_30, file = "R_files/ssf_input_all_df_60_30_50_15spp.RData")

#rename columns
colnames(used_av_df_60_30)[c(14,15)] <- c("location-long","location-lat")

write.csv(used_av_df_60_30, "R_files/ssf_input_df_60_30_50_15spp.csv")


#create 3 files
df_1 <- used_av_df_60_30 %>% 
  slice(1:(nrow(used_av_df_60_30)/3))
write.csv(df_1, "R_files/ssf_input_df_60_30_50_15spp_1.csv")

df_2 <- used_av_df_60_30 %>% 
  slice(((nrow(used_av_df_60_30)/3) + 1):((nrow(used_av_df_60_30)/3)*2))
write.csv(df_2, "R_files/ssf_input_df_60_30_50_15spp_2.csv")

df_3 <- used_av_df_60_30 %>% 
  slice((((nrow(used_av_df_60_30)/3)*2) + 1):nrow(used_av_df_60_30))
write.csv(df_3, "R_files/ssf_input_df_60_30_50_15spp_3.csv")


# # summary stats
# used_av_df_60_30 %>% 
#   #group_by(species) %>% 
#   group_by(group) %>% 
#   summarise(yrs_min = min(year(date_time)),
#             yrs_max = max(year(date_time)),
#             n_ind = n_distinct(ind),
#             n_tracks = n_distinct(track))

#visual inspection
data_sf <- used_av_df_60_30 %>% 
  filter(used == 1) %>% 
  st_as_sf(coords = c(2,3), crs = wgs)

mapview(data_sf, zcol = "species")

