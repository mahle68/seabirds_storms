# scrip for plotting used vs available boxplots for wandering albatross data (not subsetted for strong wind)
# Elham Nourani. Jan 25. 2021


library(tidyverse)
library(lubridate)
library(move)
library(parallel)
library(sf)
library(raster)
library(sp)
library(circular)
library(CircStats)
library(fitdistrplus)

#
wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

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


#open data
teil1 <- read.csv("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/all_WAAL/WAAL_GPS_2010-17_ParamWind.csv")
teil2 <- read.csv("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/all_WAAL/WAAL_XGPS_2016-19_ParamWind.csv")

columns_to_keep <- union(colnames(teil1), colnames(teil2))[-13] #remove julian date. format issues

waal <- teil1 %>% 
  dplyr::select(columns_to_keep) %>% 
  full_join(teil2[,columns_to_keep]) %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  filter (FlyingSitting == "flying") %>% #filter out sitting positions
  arrange(TripID,date_time)



# STEP 1: prepare alternative steps#####

#create a move object

#remove duplicated rows
duplicated_rows <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(waal$TripID),timestamps = waal$date_time,
                                                         sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows

waal <- waal[-duplicated_rows,] 

save(waal, file = "/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/R_files/waal_complete.RData")

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/waal_complete.RData")


mv <- move(x = waal$Longitude, y = waal$Latitude, time = waal$date_time, data = waal, animal = waal$TripID, proj = wgs)

hrs <- 2 #how long should the steps be? temporally
n <- 20 #how many alternative steps?


mycl <- makeCluster(detectCores() - 2)

clusterExport(mycl, c("mv", "n", "hrs", "wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(tidyverse)
  library(sf)
  library(raster)
  library(move)
  library(sp)
  library(circular)
  library(CircStats)
  library(fitdistrplus)
  library(tidyr)
})



  (start_time <- Sys.time())
  
  sp_obj_ls <- parLapply(mycl, split(mv), function(trip){
    #sp_obj_ls <- lapply(split(mv),function(trip){
    
    #--STEP 1: thin the data to n-hourly intervals
    trip_th<-trip%>%
      thinTrackTime(interval = as.difftime(hrs, units='hours'),
                    tolerance = as.difftime(15, units='mins')) #the unselected bursts are the large gaps between the selected ones
    #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst... longer gaps will divide the brusts) 
    trip_th$selected <- c(as.character(trip_th@burstId),NA) #assign selected as a variable
    trip_th$burst_id <-c(1,rep(NA,nrow(trip_th)-1)) #define value for first row
    
    if(nrow(trip_th@data) == 1){ #assign burst_ids
      trip_th@data$burst_id <- trip_th$burst_id
    } else {for(i in 2:nrow(trip_th@data)){
      
      if(i == nrow(trip_th@data)){
        trip_th@data$burst_id[i]<-NA
      } else
        if(trip_th@data[i-1,"selected"] == "selected"){
          trip_th@data$burst_id[i]<-trip_th@data[i-1,"burst_id"]
        } else {
          trip_th@data$burst_id[i]<-trip_th@data[i-1,"burst_id"]+1
        }
    }
    }
    #convert back to a move object (from move burst)
    trip_th <- as(trip_th,"Move")
    
    #trip_th <- trip_th[trip_th$selected == "selected" | is.na(trip_th$selected),]
    
    #--STEP 3: calculate step lengths and turning angles 
    #sl_ and ta_ calculations should be done for each burst. converting to a move burst doesnt make this automatic. so just split manually
    burst_ls<-split(trip_th,trip_th$burst_id)
    burst_ls<-Filter(function(x) length(x)>=3, burst_ls) #remove bursts with less than 3 observations
    
    burst_ls<-lapply(burst_ls,function(burst){
      burst$step_length<-c(distance(burst),NA) 
      burst$turning_angle<-c(NA,turnAngleGc(burst),NA)
      burst
    })
    
    #put burst_ls into one dataframe
    bursted_sp<-do.call(rbind,burst_ls)
    
    #reassign values
    
    bursted_sp$TripID<-trip@idData$TripID 
    bursted_sp$BirdID<-trip@idData$BirdID 
    
    
    bursted_sp
  }) %>% 
    Filter(function(x) nrow(x) > 1, .) #remove segments with no observation (these have only one obs due to the assignment of segment id)
  
  Sys.time() - start_time
  
  stopCluster(mycl) #54.94844 mins
  
  save(sp_obj_ls, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/sp_obj_ls_waal_all.RData")

  load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/sp_obj_ls_waal_all.RData")
  
  #only keep spatialpoindataframes (removes the empty elements)
  sp_obj_ls <- sp_obj_ls %>% 
    Filter(function(x) class(x) == "SpatialPointsDataFrame", .) 
  
  
  
#--STEP 4: estimate step length and turning angle distributions
#put everything in one df
bursted_df <- sp_obj_ls %>%  
  reduce(rbind) %>% 
  as.data.frame() %>% 
  dplyr::select(-c("coords.x1","coords.x2"))

#estimate von Mises parameters for turning angles
#calculate the averages (mu).steps: 1)convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan.OR use circular::mean.circular
mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))

#estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
sl <- bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")

#plot
#jpeg("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/ta_sl_dist_4hr.jpeg")
jpeg(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/ta_sl_dist_WAAL_all", hrs, "hr_", n, "n",".jpeg"))
#X11();
par(mfrow=c(1,2))
hist(sl,freq=F,main="",xlab = "Step length (km)")
plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                        rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")

hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
mtext(paste0("WAAL_all Step length and turning angle distributions (", hrs, "-hr)"), side = 3, outer =T,line = -4, font = 2) 
dev.off()

#--STEP 5: generate alternative steps
mycl <- makeCluster(detectCores() - 2)

clusterExport(mycl, c("sp_obj_ls", "n", "hrs", "wgs", "mu", "kappa", "fit.gamma1", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(tidyverse)
  library(sf)
  library(raster)
  library(move)
  library(sp)
  library(circular)
  library(CircStats)
  library(fitdistrplus)
  library(tidyr)
})



(start_time <- Sys.time())
used_av_trip <- parLapply(mycl, sp_obj_ls, function(trip){ #for each trip
  
  used_av_burst <- lapply(split(trip,trip$burst_id),function(burst){ #for each burst,
    
    #assign unique step id
    burst$step_id <- 1:nrow(burst)
    
    used_av_step <- lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
      
      current_point<- burst[this_point,]
      previous_point<-burst[this_point-1,] #this is the previous point, for calculating turning angle.
      used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
      
      #randomly generate n step lengths and turning angles
      rta <- as.vector(rvonmises(n = n, mu = mu, kappa = kappa)) #generate random turning angles with von mises distribution (in radians)
      rsl<-rgamma(n= n, shape=fit.gamma1$estimate[[1]], rate= fit.gamma1$estimate[[2]]) * 1000  #generate random step lengths from the gamma distribution. make sure unit is meters
      
      #calculate bearing of previous point
      #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
      prev_bearing<-NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
                                      previous_point@coords[,1], current_point@coords[,1])
      
      #find the gepgraphic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
      current_point_m <- spTransform(current_point, meters_proj) #convert to meters proj
      rnd <- data.frame(lon = current_point_m@coords[,1] + rsl*cos(rta),lat = current_point_m@coords[,2] + rsl*sin(rta)) #for this to work, lat and lon should be in meters as well. boo. coordinates in meters?
      
      #covnert back to lat-lon proj
      rnd_sp<-rnd
      coordinates(rnd_sp)<-~lon+lat
      proj4string(rnd_sp)<-meters_proj
      rnd_sp<-spTransform(rnd_sp,wgs)
      
      #put used and available points together
      df <- used_point@data %>%  
        slice(rep(row_number(),n+1)) %>% #paste each row 20 times for the used and alternative steps
        mutate(Longitude = c(head(Longitude,1),rnd_sp@coords[,1]),
               Latitude = c(head(Latitude,1),rnd_sp@coords[,2]),
               used = c(1,rep(0,n)))  %>% #one hour after the start point of the step
        rowwise() %>% 
        mutate(heading = NCEP.loxodrome.na(lat1=current_point$Latitude,lat2=Latitude,lon1=current_point$Longitude,lon2= Longitude)) %>% 
        #dplyr::select(-c("accumulated_dist", "turningAngle_deg", "turningAngle_rad", "speed_kmh")) %>% 
        as.data.frame()
      
      df[df$used == 0, c("turning_angle", "step_length")] <- NA
      
      df
      
    }) %>% 
      reduce(rbind)
    used_av_step
  }) %>% 
    reduce(rbind)
  used_av_burst
}) %>% 
  reduce(rbind)
#used_av_trip
#})
Sys.time() - start_time


save(used_av_trip, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_all_alt_steps_", hrs, "hr_", n, "n",".RData"))

# #plot
# jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/alt_step_eg.jpeg")
# #X11()
# maps::map("world",xlim = c(49,57), ylim = c(-49,-44.9),fil = TRUE,col = "ivory") #flyway
# title("Example of used and available steps")
# legend("bottomleft", legend = c("track","previous point", "start point","alternative end points", "used end point"),
#       col = c("grey", "yellow4", "firebrick1","orange", "firebrick4"), pch = 16, bty = "n")
# points(burst,col = "grey", pch = 16, cex = 1)
# points(previous_point,col = "yellow4", pch = 16, cex = 2)
# points(current_point,col = "firebrick1", pch = 16, cex = 2)
# points(rnd_sp, col = "orange", pch = 16, cex = 1)
# points(used_point, col = "firebrick4", pch = 16, cex = 2)
# dev.off()

#have a look
X11();par(mfrow= c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0))
maps::map("world",fil = TRUE,col = "grey85", border=NA) 
points(used_av_trip[used_av_trip$used == 0,c("Longitude","Latitude")], pch = 16, cex = 0.2, col = "gray55")
points(used_av_trip[used_av_trip$used == 1,c("Longitude","Latitude")], pch = 16, cex = 0.2, col = "orange")


#prep for movebank annotation
used_av_trip <- used_av_trip %>% 
  mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% 
  as.data.frame()

colnames(used_av_trip)[c(3,2)] <- c("location-long","location-lat") #rename columns to match movebank format
write.csv(used_av_trip, paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_all_alt_steps_", hrs, "hr_", n, "n",".csv"))


## STEP 2: data exploration!#####

#open annotated data
ann <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/waal_all_2hr/WAAL_all_alt_steps_2hr_20n.csv-30396442191641561.csv", stringsAsFactors = F)

ann_cmpl <-  ann %>% 
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u10 = ECMWF.Interim.Full.Daily.SFC.Wind..10.m.above.Ground.U.Component.,
         v10 = ECMWF.Interim.Full.Daily.SFC.Wind..10.m.above.Ground.V.Component.,
         sea_s_pr = ECMWF.Interim.Full.Daily.SFC.Mean.Sea.Level.Pressure,
         wave_h = ECMWF.Interim.Full.Daily.SFC.Significant.Wave.Height,
         air_pr = ECMWF.Interim.Full.Daily.SFC.Surface.Air.Pressure) %>% 
  mutate(delta_t = sst - t2m,
         wind_support_ms = wind_support(u=u10,v=v10,heading=heading),
         cross_wind_ms = cross_wind(u=u10,v=v10,heading=heading),
         abs_cross_wind_ms = abs(cross_wind(u=u10,v=v10,heading=heading)),
         wind_speed_ms = sqrt(u10^2 + v10^2),
         stratum = paste(TripID, burst_id, step_id, sep = "_")) %>% 
  mutate(wind_support_kmh = wind_support_ms * 3.6,
         cross_wind_kmh = cross_wind_ms * 3.6,
         wind_speed_kmh = wind_speed_ms * 3.6,
         abs_cross_wind_kmh = abs_cross_wind_ms * 3.6) %>% 
  as.data.frame()

save(ann_cmpl, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_all_ssf_ann_2hr.RData")

# 
# #plots
# 
# #boxplots
# jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/Used_avail_WAAL_6hr.jpeg", width = 10, height = 4, units = "in", res = 300)
# #X11(width = 13, height = 4)
# par(mfrow= c(1,3), oma = c(0,0,3,0))
# for(i in c("wind_speed_kmh", "cross_wind_kmh", "wind_support_kmh")){ #tried with non-interpolated wind support and crosswind and point of sail, but no difference
#   boxplot(ann_cmpl[,i] ~ ann_cmpl[,"used"], 
#           xaxt = "n", boxfill = c("gray","orange"), main = i, xlab = "", ylab = "")
#   axis(1, labels = c("available", "used"), at = c(1,2), tck = 0)
# }
# mtext("Tracks containing 50 km/h winds- 6hrly steps", side = 3, outer = T, cex = 1)
# 
# dev.off()


