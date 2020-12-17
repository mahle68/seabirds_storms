#script for subsetting tracking data and retaining only segments with winds higher than x kmh
#Dec. 15, 2020. Elham Nourani. Radolfzell am Bodensee

library(tidyverse)
library(lubridate)
library(mapview)
library(sf)
library(sp)
library(move)
library(circular)
library(CircStats)
library(fitdistrplus)

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



#----------- STEP 1: open albatross data ----

#one-hourly data annotated by Sophie
data <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/DATSETS/WAAL_allGPS_2010-2020_homo_R1h_TrackParam_Wind50kmh.csv",
                 na.strings = c("NA",""), fileEncoding="latin1") %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  filter (FlyingSitting == "flying") %>% #filter out sitting positions
  arrange(TripID,date_time) %>% 
  #rowwise() %>% 
  mutate(row_id = row_number())

#----------- STEP 2: subset points with strong wind ----

windy <- data %>% 
  filter(windSpeed_kmh >= 50) %>% 
  mutate(delta = lead(row_id,1) - row_id ) 

#----------- STEP 3: assign segment IDs to consecutive points retained ----

windy$segID <- c(1,rep(NA,nrow(windy)-1))
  
for(i in 2:nrow(windy)){
  if(windy$delta[i-1] == 1){
    windy$segID[i] <- windy$segID[i-1]
  } else {
    windy$segID[i] <- windy$segID[i-1] + 1
  }
}

#----------- STEP 4: keep segments with at least n points ----

segs_to_keep <- windy %>% 
  group_by(segID) %>% 
  summarise(n_pts = length(segID)) %>% 
  filter(n_pts >= 10)

windy_sf <- windy %>% 
  filter(segID %in% segs_to_keep$segID) %>% 
  st_as_sf(coords = c("Longitude", "Latitude")) %>% 
  st_set_crs(wgs)

mapview(windy_sf,zcol = "segID")

#----------- STEP 5: sub-sample to n hours ----

hrs <- 2 #how long should the steps be? temporally

#selecting every nth row in each segment
windy_2hr <- windy_sf %>% 
  group_by(segID) %>% 
  filter(row_number() %% hrs == 1) %>% 
  ungroup() %>% 
  arrange(segID, date_time) #this is necessary for making a move object later

mapview(windy_2hr,zcol = "segID")

#----------- STEP 6: calculate step lengths and turning angles ----

#create a move object

mv <- move(x = st_coordinates(windy_2hr)[,1], y = st_coordinates(windy_2hr)[,2], time = windy_2hr$date_time, 
           data = st_drop_geometry(windy_2hr), animal = windy_2hr$segID, proj = wgs)

#sl_ and ta_ calculations should be done for each burst. converting to a move burst doesnt make this automatic. so just split manually
burst_ls <- split(mv)
burst_ls <- Filter(function(x) length(x) >= 3, burst_ls) #remove bursts with less than 3 observations

#calculate step lengths and turning angles
burst_ls <- lapply(burst_ls,function(burst){
  burst$step_length <- c(distance(burst), NA) 
  burst$turning_angle <- c(NA,turnAngleGc(burst), NA)
  burst$segID <- burst@idData$segID
  burst
})

bursted_sp <- do.call(rbind,burst_ls)

#----------- STEP 7: estimate step length and turning angle distributions ----

#create a df
bursted_df <- as.data.frame(bursted_sp) %>% 
  filter(step_length < 500000) %>% #filter out points with over 500 km step length. outliers 

#estimate von Mises parameters for turning angles
#calculate the averages (mu).steps: 1)convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan.OR use circular::mean.circular
mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))

#estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
sl <- bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")

#plot
#jpeg("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/ta_sl_dist_4hr.jpeg")
#jpeg(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/ta_sl_dist_WAAL_", hrs, "hr_", n, "n",".jpeg"))
#X11();
par(mfrow=c(1,2))
hist(sl,freq=F,main="",xlab = "Step length (km)")
plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                        rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")

hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
mtext(paste0("WAAL_strong winds_ Step length and turning angle distributions (", hrs, "-hr)"), side = 3, outer =T,line = -4, font = 2) 
dev.off()

#----------- STEP 8: generate alternative steps ----

n <- 20 #how many alternative steps?

bursted_sp <- bursted_df
coordinates(bursted_sp) <- ~ coords.x1 + coords.x2
proj4string(bursted_sp) <- wgs

used_av_burst <- lapply(split(bursted_sp,bursted_sp$segID),function(burst){ #for each burst (segID),
  
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
    df <- as.data.frame(used_point) %>%  
      slice(rep(row_number(),n+1)) %>% #paste each row 20 times for the used and alternative steps
      mutate(Longitude = c(head(coords.x1,1),rnd_sp@coords[,1]),
             Latitude = c(head(coords.x2,1),rnd_sp@coords[,2]),
             used = c(1,rep(0,n)))  %>% #one hour after the start point of the step
      rowwise() %>% 
      mutate(heading = NCEP.loxodrome.na(lat1 = current_point@coords[,2], lat2 = Latitude, lon1 = current_point@coords[,1], lon2 = Longitude)) %>% 
      #dplyr::select(-c("accumulated_dist", "turningAngle_deg", "turningAngle_rad", "speed_kmh")) %>% 
      as.data.frame()
    
    # df[df$used == 0, c("turning_angle", "step_length", "angleBirdWind", "windAzimuth", "windSpeed_kmh", "windDirection", "windSpeed_ms", "v_wind_interp", 
    #                   "u_wind_interp", "v_wind", "t", "travelled_distance_km", "DistColo")] <- NA
    
    df
    
  }) %>% 
    reduce(rbind)
  used_av_step
}) %>% 
  reduce(rbind)
used_av_burst

#have a look
X11();par(mfrow= c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0))
maps::map("world",fil = TRUE,col = "grey85", border=NA) 
points(used_av_burst[used_av_burst$used == 0,c("Longitude","Latitude")], pch = 16, cex = 0.2, col = "gray55")
points(used_av_burst[used_av_burst$used == 1,c("Longitude","Latitude")], pch = 16, cex = 0.2, col = "orange")

save(used_av_burst, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_onlyabove50kmh_alt_steps_", hrs, "hr_", n, "n",".RData"))

#----------- STEP 9: annotation ----

#prep for movebank annotation
used_av_burst <- used_av_burst %>% 
  mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% 
  as.data.frame()

colnames(used_av_burst)[c(30,31)] <- c("location-long","location-lat") #rename columns to match movebank format
write.csv(used_av_burst, paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_onlyabove50kmh_alt_steps_", hrs, "hr_", n, "n",".csv"))

#----------- STEP 10: data exploration ----



