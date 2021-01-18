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
library(corrr)
library(INLA)

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

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/Peter_Ryan_data_annotated_SplitTrip_filterWind50kmh.RData") #PR_tot
alb <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/DATSETS/WAAL_allGPS_2010-2020_homo_R1h_TrackParam_Wind50kmh.csv",
                na.strings = c("NA",""), fileEncoding="latin1") %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         scientific_name = "Diomedea exulans",
         common_name = "Wandering albatross") %>%
  rename(location.long = Longitude,
         location.lat = Latitude)
filter (FlyingSitting == "flying") %>% #filter out sitting positions
arrange(TripID,date_time)

#put together
cols <- intersect(names(alb), names(PR_tot))

data <- PR_tot[,cols] %>%
  filter(common_name %in% c("Atlantic Petrel","Soft-plumaged Petrel")) %>%  #keep the species of interest
  full_join(alb[,c(cols, "BirdID")])

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
ann <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/windy_WAAL_2hr/WAAL_onlyabove50kmh_alt_steps_2hr_20n.csv-517941447563460742/WAAL_onlyabove50kmh_alt_steps_2hr_20n.csv-517941447563460742.csv", stringsAsFactors = F)

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
         stratum = paste(TripID, step_id, sep = "_")) %>% #here, each trip is considered a burst, so there is no burst id
  mutate(wind_support_kmh = wind_support_ms * 3.6,
         cross_wind_kmh = cross_wind_ms * 3.6,
         wind_speed_kmh = wind_speed_ms * 3.6,
         abs_cross_wind_kmh = abs(cross_wind_ms * 3.6)) %>% 
  as.data.frame()

save(ann_cmpl, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/windy_WAAL_ssf_ann_2hr.RData")

#boxplots
jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/Used_avail_windy_WAAL_2hr.jpeg", width = 10, height = 4, units = "in", res = 300)
#X11(width = 13, height = 4)
par(mfrow= c(1,3), oma = c(0,0,3,0))
for(i in c("wind_speed_kmh", "abs_cross_wind_kmh", "wind_support_kmh")){ #tried with non-interpolated wind support and crosswind and point of sail, but no difference
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"used"], 
          xaxt = "n", boxfill = c("gray","orange"), main = i, xlab = "", ylab = "")
  axis(1, labels = c("available", "used"), at = c(1,2), tck = 0)
}
mtext("Tracks over 50 km/h winds- 2hrly steps", side = 3, outer = T, cex = 1)

dev.off()


#----------- STEP 11: analysis ----

#correlation
ann_cmpl %>% 
  dplyr::select(c("wind_speed_kmh","abs_cross_wind_kmh","wind_support_kmh")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: var_cw with location.lat and var_delta_t with location.lat

#z-transform and bin
all_data <- ann_cmpl %>% 
  mutate_at(c("wind_speed_kmh","cross_wind_kmh","abs_cross_wind_kmh","wind_support_kmh"),
            list(z = ~scale(.))) %>%
  #bin the data for smooth terms, otherwise I get an error that locations are too close.
  mutate_at(c("wind_speed_kmh","abs_cross_wind_kmh","wind_support_kmh"),
            list(group = ~inla.group(.,n = 50, method = "cut"))) %>%
  as.data.frame()

#repeat the individual ID column for INLA

all_data <- all_data %>% 
  mutate(BirdID1 = factor(BirdID),
         BirdID2 = factor(BirdID),
         BirdID3 = factor(BirdID),
         BirdID4 = factor(BirdID))

# set mean and precision for the priors of slope coefficients (fixed effects)
mean.beta <- 0
prec.beta <- 1e-4 #precision of 1e-4 equals a variance of 1e4 ;)

#model with smooth terms for wind support and crosswind. omit wind speed

model_title <- "Smooth terms for winds support, abs(cross wind) (binned data; no z_transformation)"

formula <- used ~ -1 +
  f(wind_support_kmh_group, model = "rw2", constr = F) + 
  f(abs_cross_wind_kmh_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T)))

(b <- Sys.time())
m4 <- inla(formula, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data,
           num.threads = 10,
           control.predictor = list(compute = T), #list(link = 1), #link is only relevant for NA observations. required to set the right link (i.e., the logit function) 
           #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
           control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))
Sys.time() - b 

summary(m4) #WAIC = 10274.53 ; MLik = -12130.70

save(m4, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/windy_WAAL_wspt_cw_50_cuts.RData") #n of 50 for binned wind

#plot
wspt <- data.frame(x = m4$summary.random$wind_support_kmh_group[, "ID"],
                   y = m4$summary.random$wind_support_kmh_group[, "mean"],
                   ll95 = m4$summary.random$wind_support_kmh_group[,"0.025quant"],
                   ul95 = m4$summary.random$wind_support_kmh_group[,"0.975quant"]
)

cw <- data.frame(x = m4$summary.random$abs_cross_wind_kmh_group[, "ID"],
                 y = m4$summary.random$abs_cross_wind_kmh_group[, "mean"],
                 ll95 = m4$summary.random$abs_cross_wind_kmh_group[,"0.025quant"],
                 ul95 = m4$summary.random$abs_cross_wind_kmh_group[,"0.975quant"]
)

#plot with raw value and credible intervals
jpeg(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/INLA_results_figures/windy_WAAL_2hrs_smooth_wspt_abs_cw.jpeg"), 
     width = 10, height = 4, units = "in", res = 300)
#X11(width = 8, height = 4)
par(mfrow= c(1,2), oma = c(0,1,3,0), bty = "l")

#plot(x = ann_cmpl$wind_speed_kmh, y = ann_cmpl$used, pch = 16,  col = adjustcolor("grey", alpha.f = 0.1), xlab = "wind speed (kmh)", ylab = "exp(y)")
#lines(wspd$x,exp(wspd$y)) 
#polygon(x = c(wspd$x, rev(wspd$x)), y = c(exp(wspd$ll95),rev(exp(wspd$ul95))), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

plot(x = ann_cmpl$wind_support_kmh, y = ann_cmpl$used, pch = 16,  col = adjustcolor("grey", alpha.f = 0.1), xlab = "wind support (kmh)", ylab = "exp(y)")
lines(wspt$x,exp(wspt$y)) 
polygon(x = c(wspt$x, rev(wspt$x)), y = c(exp(wspt$ll95),rev(exp(wspt$ul95))), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

plot(x = ann_cmpl$abs_cross_wind_kmh, y = ann_cmpl$used, pch = 16,  col = adjustcolor("grey", alpha.f = 0.1), xlab = "abs(cross wind) (kmh)", ylab = "exp(y)")
lines(cw$x,exp(cw$y)) 
polygon(x = c(cw$x, rev(cw$x)), y = c(exp(cw$ll95),rev(exp(cw$ul95))), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

mtext (paste0("INLA SSF analysis. Wandering Albatross (segments with wind >= 50 kmh)"), side = 3, outer = T, line = -0.2)
mtext (paste0("2-hrly steps. ", model_title), 
       side = 3, outer = T, line = -2)
dev.off()
