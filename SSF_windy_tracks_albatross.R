#script for exploration of seabird data to find their response to strong wind conditions
#Elham Nourani, PhD. Nov. 12. 2020. Radolfzell am Bodensee, Germany
#update Nov 26: filter out sitting positions and redo the alternative steps generation
# see previous versions for hourly data (annotated by Sophie)

library(tidyverse)
library(move)
library(sf)
library(sp)
library(circular)
library(CircStats)
library(fitdistrplus)
library(RNCEP)
library(lubridate)
library(mapview)
library(parallel)
library(tidyr)
library(corrr)
library(lme4)
library(MuMIn)
library(mgcv)
library(survival)
library(INLA)
library(ggregplot) #devtools::install_github("gfalbery/ggregplot")
library(maptools)
library(brinla)

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

#open data: prepared by Sophie

data <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/DATSETS/WAAL_allGPS_2010-2020_homo_R1h_TrackParam_Wind50kmh.csv",
                 na.strings = c("NA",""), fileEncoding="latin1") %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  filter (FlyingSitting == "flying") %>% #filter out sitting positions
  arrange(TripID,date_time)

# STEP 1: prepare alternative steps#####

#create a move object

mv <- move(x = data$Longitude,y = data$Latitude,time = data$date_time, data = data, animal = data$TripID,proj = wgs)

hrs <- 6 #how long should the steps be? temporally
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



start_time <- Sys.time()

sp_obj_ls <- parLapply(mycl, split(mv), function(trip){
#sp_obj_ls <- lapply(split(mv),function(trip){
      
      #--STEP 1: thin the data to 1-hourly intervals
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
      
      #if(length(bursted_sp) >= 1){
      #  bursted_sp$track<-trip@idData$track
      #  bursted_sp$group<-trip@idData$group
      #}
      
      bursted_sp$TripID<-trip@idData$TripID 
      bursted_sp$BirdID<-trip@idData$BirdID 
      bursted_sp$sex<-trip@idData$sex 
      
      bursted_sp
    }) %>% 
      Filter(function(x) length(x) > 1, .) #remove segments with no observation (these have only one obs due to the assignment of segment id)
  
Sys.time() - start_time
  
stopCluster(mycl)
    
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
    jpeg(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/ta_sl_dist_WAAL_", hrs, "hr_", n, "n",".jpeg"))
    #X11();
    par(mfrow=c(1,2))
    hist(sl,freq=F,main="",xlab = "Step length (km)")
    plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                            rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
    
    hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
    plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
    mtext(paste0("WAAL_Step length and turning angle distributions (", hrs, "-hr)"), side = 3, outer =T,line = -4, font = 2) 
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
    
    
    
    start_time <- Sys.time()
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
            dplyr::select(-c("accumulated_dist", "turningAngle_deg", "turningAngle_rad", "speed_kmh")) %>% 
            as.data.frame()
          
          df[df$used == 0, c("turning_angle", "step_length", "angleBirdWind", "windAzimuth", "windSpeed_kmh", "windDirection", "windSpeed_ms", "v_wind_interp", 
                             "u_wind_interp", "v_wind", "t", "travelled_distance_km", "DistColo")] <- NA
          
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
    

save(used_av_trip, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_Wind50kmh_alt_steps_", hrs, "hr_", n, "n",".RData"))

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

colnames(used_av_trip)[c(1,2)] <- c("location-long","location-lat") #rename columns to match movebank format
write.csv(used_av_trip, paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_Wind50kmh_alt_steps_", hrs, "hr_", n, "n",".csv"))


# STEP 2: data exploration!#####

#open annotated data
ann <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/WAAL_6hr/WAAL_Wind50kmh_alt_steps_6hr_20n.csv-3127490191197375873.csv", stringsAsFactors = F)

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
         wind_speed_kmh = wind_speed_ms * 3.6) %>% 
  as.data.frame()

save(ann_cmpl, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_data_ssf_ann_6hr.RData")


#plots

#boxplots
jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/Used_avail_WAAL_6hr.jpeg", width = 10, height = 4, units = "in", res = 300)
#X11(width = 13, height = 4)
par(mfrow= c(1,3), oma = c(0,0,3,0))
for(i in c("wind_speed_kmh", "cross_wind_kmh", "wind_support_kmh")){ #tried with non-interpolated wind support and crosswind and point of sail, but no difference
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"used"], 
          xaxt = "n", boxfill = c("gray","orange"), main = i, xlab = "", ylab = "")
  axis(1, labels = c("available", "used"), at = c(1,2), tck = 0)
}
mtext("Tracks containing 50 km/h winds- 6hrly steps", side = 3, outer = T, cex = 1)

dev.off()

#density plots
X11()

plot(density(ann_cmpl[ann_cmpl$used == 0,"wind_support_kmh"]), col= "grey")
lines(density(ann_cmpl[ann_cmpl$used == 1,"wind_support_kmh"]), col = "orange")

plot(density(ann_cmpl[ann_cmpl$used == 1,"cross_wind_kmh"]))
lines(density(ann_cmpl[ann_cmpl$used == 0,"cross_wind_kmh"]), col= "grey")

plot(density(ann_cmpl[ann_cmpl$used == 1,"wind_speed_kmh"]), col= "red")
lines(density(ann_cmpl[ann_cmpl$used == 0,"wind_speed_kmh"]), col= "blue")


hist(ann_cmpl[ann_cmpl$used == 0,"wind_support_kmh"], col = "yellow", freq = F)
hist(ann_cmpl[ann_cmpl$used == 1,"wind_support_kmh"], add= T, freq = F)

# STEP 3: analysis#####

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_data_ssf_ann_6hr.RData")

#correlation
ann_cmpl %>% 
  dplyr::select(c("wind_speed_kmh","cross_wind_kmh","wind_support_kmh")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: var_cw with location.lat and var_delta_t with location.lat

#z-transform
all_data <- ann_cmpl %>% 
  mutate(abs_cw = abs(cross_wind_kmh)) %>% 
  #group_by(species) # z_transform for each species separately. or not? ... huh!
  mutate_at(c("wind_speed_kmh","cross_wind_kmh","abs_cw","wind_support_kmh"),
            list(z = ~scale(.))) %>%
  as.data.frame()

#repeat the individual ID column for INLA

all_data <- all_data %>% 
  mutate(BirdID1 = factor(BirdID),
         BirdID2 = factor(BirdID),
         BirdID3 = factor(BirdID),
         BirdID4 = factor(BirdID))

#work on a sample of individuals for now
sample <- all_data[all_data$BirdID %in% c("BS8484", "BS27620"),]

#---- clogit: run a quick ssf and see what happens ----
form1 <- formula(used ~ wind_speed_kmh_z  + abs_cw_z + wind_support_kmh_z + I(wind_support_kmh_z ^ 2) +
                   strata(stratum))
m1 <- clogit(form1, data = all_data)

form2 <- formula(used ~ wind_speed_kmh_z  + cross_wind_kmh_z + wind_support_kmh_z + I(wind_support_kmh_z ^ 2) +
                   strata(stratum))
m2 <- clogit(form2, data = all_data)

d <- fit_clogit(all_data, form2, model = TRUE)

log_rss(x1 = all_data, object = d) #relative selection strength

#---- INLA: also try adding step length and maybe cos(ta_) ----

# Set mean and precision for the priors of slope coefficients (fixed effects)
mean.beta <- 0
prec.beta <- 1e-4 #precision of 1e-4 equals a variance of 1e4 ;)


#visualize the priors (stats rethinking p. 349)
# N <- 100
# a <- rnorm(N, 0,1e6) #prior of the intercept
# b <- rnorm(N, 3, 0.05)
# plot(NULL, xlim = c(-2,2), ylim = c(0,100))
# for(i in 1:N) curve(exp(a[i] + b[i]*x), add = T)


# prec <- seq(0.01, 2, by = 0.01)
# alpha <- c(0.95, 0.90, 0.80, 0.50, 0.25, 0.1)
# pc.prec <- lapply(alpha, function(a) {
#   inla.pc.dprec(prec, 1, a)
# })
# pc.prec <- do.call(c, pc.prec)
# 
# tab <- data.frame(prec = rep(prec, length(alpha)), pc.prec = pc.prec,
#                   alpha = rep(as.character(alpha), each = length(prec)))
# 

library(ggplot2)
ggplot(tab, aes(x = prec, y = pc.prec, linetype = alpha)) + geom_line() +
  xlab(expression(paste("precision ", tau))) +
  ylab("density") +
  ggtitle(expression(paste("PC prior on the precision ", tau,
                           " with P(", sigma, " > 1) = ", alpha))) +
  scale_linetype_discrete(name = expression(alpha))

#model 1: plain. just wspd, cw, and wspt
formula1 <- used ~ -1 + wind_speed_kmh_z + cross_wind_kmh_z + wind_support_kmh_z + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) #+
  #f(BirdID1, wind_speed_kmh_z, model = "iid", 
  #  hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  #f(BirdID2, cross_wind_kmh_z,  model = "iid",
  #  hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  #f(BirdID3, wind_support_kmh_z,  model = "iid",
  #  hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
m1 <- inla(formula1, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = all_data,
            num.threads = 10,
            control.predictor = list(link = 1), #required to set the right link (i.e., the logit function) 
            #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))
Sys.time() - b

summary(m1) #WAIC 90050.57; Mlik -96029.98

save(m1, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/m1_all_wind.RData")
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/m1_all_wind.RData")


#---------------model 2: only wind support as a smooth term (not centered). sample of data
#add missing predictor values to the dataset, for prediction and plotting
#xx <- seq(-90, 90, by = 1)
#new_data <- cbind(wind_support_kmh = xx, used = NA)
#new_data <- rbind(all_data, new_data)

#bin the data, otherwise I get an error that locations are too close.
sample$wspt_group <- inla.group(sample$wind_support_kmh, n = 20, method = "quantile")

formula2 <- used ~ -1 + f(wspt_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) #+ 
  #f(BirdID1, wspt_group, model = "iid", 
  #  hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
m2 <- inla(formula2, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data,
           num.threads = 10,
           control.predictor = list(compute = T), #list(link = 1), #link is only relevant for NA observations. required to set the right link (i.e., the logit function) 
           #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
           control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))
Sys.time() - b

summary(m2)


#---------------model 2a: full model 1 showed that wind speed had the highest coefficient. So, look at wind speed as a smooth term (not centered). all data
#bin the data, otherwise I get an error that locations are too close.
all_data$wspd_group <- inla.group(all_data$wind_speed_kmh, n = 50, method = "quantile") 
all_data$wspt_group <- inla.group(all_data$wind_support_kmh, n = 50, method = "quantile")
all_data$cw_group <- inla.group(all_data$cross_wind_kmh, n = 50, method = "quantile")
all_data$abs_cw_group <- inla.group(all_data$abs_cw, n = 50, method = "quantile")

formula2a <- used ~ -1 + f(wspd_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T)))

(b <- Sys.time())
m2a <- inla(formula2a, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data,
           num.threads = 10,
           control.predictor = list(compute = T), #list(link = 1), #link is only relevant for NA observations. required to set the right link (i.e., the logit function) 
           #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
           control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))
Sys.time() - b

save(m2a, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/m2a_full.RData")

summary(m2a)

#plot the smooth effect of wind support
## ----label = "lidarsmooth", fig = TRUE, echo = FALSE, fig.cap = '(ref:lidarsmooth)'----

tab.rw2 <- data.frame(x = m2a$summary.random$wspd_group[, "ID"],
                      y = m2a$summary.random$wspd_group[, "mean"]
)

plot(tab.rw2$x,tab.rw2$y) #I don't understand what the units of the y axis are. if I exp() them, they will go from 0.02 to 0.035
plot(tab.rw2$x,exp(tab.rw2$y)) 

#ggplot(aes(x = wind_support_kmh, y = used), data = sample) + 
  #geom_point(aes(x = wspt_group, y = used), colour = "grey") +
  #geom_point() +
  #ggtitle("wind support") +
  #geom_line(aes(x = x, y = y, linetype = "solid"), data = tab.rw2) #+
  #geom_line(aes(x = x, y = y, linetype = "dashed"), data = tab.rw2) +
  #scale_linetype_manual(name = "Smoothing method (grouped range)",
  #                      values = c("solid", "dashed"),
  #                      labels = c("rw2", "rw1")) +
  #theme(legend.position = "bottom")

marginals_mean <- sapply(m2$marginals.hyperpar,
                         function(x)
                           inla.emarginal(function(x) x, inla.tmarginal(function(x) 1/x, x)))

names(marginals_mean) <- sapply(as.vector(as.character(names(marginals_mean))),
                                function(y) gsub("Precision", x=y, "Mean of variance"))
marginals_mean


m = m2a$marginals.random$wspd_group
plot(inla.smarginal(m), type="l")


#---------------model 3: full model with 3 smooth terms (not centered; no random effect for individuals- n = 193)
formula3 <- used ~ -1 + f(wspd_group, model = "rw2", constr = F) + 
  f(wspt_group, model = "rw2", constr = F) + 
  f(cw_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T)))

(b <- Sys.time())
m3 <- inla(formula3, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = all_data,
            num.threads = 10,
            control.predictor = list(compute = T), #list(link = 1), #link is only relevant for NA observations. required to set the right link (i.e., the logit function) 
            #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))
Sys.time() - b #2.5253 hours

summary(m3) #WAIC 85746.35; MLik -113256.38

save(m3, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/m3_full.RData") #n of 50 for binned wind
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/m3_full.RData")


wspd <- data.frame(x = m3$summary.random$wspd_group[, "ID"],
                      y = m3$summary.random$wspd_group[, "mean"],
                   ll95 = m3$summary.random$wspd_group[,"0.025quant"],
                   ul95 = m3$summary.random$wspd_group[,"0.975quant"]
)

wspt <- data.frame(x = m3$summary.random$wspt_group[, "ID"],
                   y = m3$summary.random$wspt_group[, "mean"],
                   ll95 = m3$summary.random$wspt_group[,"0.025quant"],
                   ul95 = m3$summary.random$wspt_group[,"0.975quant"]
)

cw <- data.frame(x = m3$summary.random$cw_group[, "ID"],
                   y = m3$summary.random$cw_group[, "mean"],
                 ll95 = m3$summary.random$cw_group[,"0.025quant"],
                 ul95 = m3$summary.random$cw_group[,"0.975quant"]
)

jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/INLA_results_figures/WAAL_6hr_all_smooth.jpeg", width = 11, height = 3.5, units = "in", res = 300)
#X11(width = 11, height = 3.5)
par(mfrow= c(1,3), oma = c(0,1,3,0), bty = "l")
plot(wspd$x,exp(wspd$y), xlab = "wind speed (kmh)", ylab = "exp(y)", type = "l") 
points(x = all_data$wspd_group, y = all_data$used)
plot(wspt$x,exp(wspt$y), xlab = "wind support (kmh)", ylab = "", type = "l") 
abline(v = 0, col = "gray", lty = 2)
plot(cw$x,exp(cw$y), xlab = "cross wind (kmh)", ylab = "", type = "l") 
abline(v = 0, col = "gray", lty = 2)
mtext ("INLA SSF analysis. Wandering Albatross.", side = 3, outer = T, line = -0.2)
mtext ("6-hrly steps. Smooth terms for wind speed, wind support, and cross wind (binned data; no z-transofmration). no random effects.", 
       side = 3, outer = T, line = -2)
dev.off()

#plot with raw value and credible intervals
jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/INLA_results_figures/WAAL_6hr_all_smooth_0_1.jpeg", width = 11, height = 3.5, units = "in", res = 300)
#X11(width = 11, height = 3.5)
par(mfrow= c(1,3), oma = c(0,1,3,0), bty = "l")

plot(x = all_data$wspd_group, y = all_data$used, col = "grey90", xlab = "wind speed (kmh)", ylab = "exp(y)")
lines(wspd$x,exp(wspd$y)) 
polygon(x = c(wspd$x, rev(wspd$x)), y = c(exp(wspd$ll95),exp(wspd$ul95)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

plot(x = all_data$wspt_group, y = all_data$used, col = "grey90", xlab = "wind support (kmh)", ylab = "exp(y)")
lines(wspt$x,exp(wspt$y)) 
polygon(x = c(wspt$x, rev(wspt$x)), y = c(exp(wspt$ll95),exp(wspt$ul95)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)


plot(x = all_data$cw_group, y = all_data$used, col = "grey90", xlab = "cross wind (kmh)", ylab = "exp(y)")
lines(cw$x,exp(cw$y)) 
polygon(x = c(cw$x, rev(cw$x)), y = c(exp(cw$ll95),exp(cw$ul95)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

mtext ("INLA SSF analysis. Wandering Albatross.", side = 3, outer = T, line = -0.2)
mtext ("6-hrly steps. Smooth terms for wind speed, wind support, and cross wind (binned data; no z-transofmration). no random effects.", 
       side = 3, outer = T, line = -2)
dev.off()


#---------------model 4: full model with 3 smooth terms. absolute value of wind support 
formula4 <- used ~ -1 + f(wspd_group, model = "rw2", constr = F) + 
  f(wspt_group, model = "rw2", constr = F) + 
  f(abs_cw_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T)))

(b <- Sys.time())
m4 <- inla(formula4, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data,
           num.threads = 10,
           control.predictor = list(compute = T), #list(link = 1), #link is only relevant for NA observations. required to set the right link (i.e., the logit function) 
           #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
           control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))
Sys.time() - b 

summary(m4) #WAIC = 85827.00 ; MLik = -113325.29

save(m4, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/m4_full.RData") #n of 50 for binned wind

wspd <- data.frame(x = m4$summary.random$wspd_group[, "ID"],
                   y = m4$summary.random$wspd_group[, "mean"],
                   ll95 = m4$summary.random$wspd_group[,"0.025quant"],
                   ul95 = m4$summary.random$wspd_group[,"0.975quant"]
)

wspt <- data.frame(x = m4$summary.random$wspt_group[, "ID"],
                   y = m4$summary.random$wspt_group[, "mean"],
                   ll95 = m4$summary.random$wspt_group[,"0.025quant"],
                   ul95 = m4$summary.random$wspt_group[,"0.975quant"]
)

cw <- data.frame(x = m4$summary.random$abs_cw_group[, "ID"],
                 y = m4$summary.random$abs_cw_group[, "mean"],
                 ll95 = m4$summary.random$abs_cw_group[,"0.025quant"],
                 ul95 = m4$summary.random$abs_cw_group[,"0.975quant"]
)

#plot with raw value and credible intervals
jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/INLA_results_figures/WAAL_6hr_all_smooth_0_1_abs_cw.jpeg", width = 11, height = 3.5, units = "in", res = 300)
#X11(width = 11, height = 3.5)
par(mfrow= c(1,3), oma = c(0,1,3,0), bty = "l")

plot(x = all_data$wspd_group, y = all_data$used, col = "grey90", xlab = "wind speed (kmh)", ylab = "exp(y)")
lines(wspd$x,exp(wspd$y)) 
polygon(x = c(wspd$x, rev(wspd$x)), y = c(exp(wspd$ll95),exp(wspd$ul95)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

plot(x = all_data$wspt_group, y = all_data$used, col = "grey90", xlab = "wind support (kmh)", ylab = "exp(y)")
lines(wspt$x,exp(wspt$y)) 
polygon(x = c(wspt$x, rev(wspt$x)), y = c(exp(wspt$ll95),exp(wspt$ul95)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)


plot(x = all_data$abs_cw_group, y = all_data$used, col = "grey90", xlab = "abs(cross wind) (kmh)", ylab = "exp(y)")
lines(cw$x,exp(cw$y)) 
polygon(x = c(cw$x, rev(cw$x)), y = c(exp(cw$ll95),exp(cw$ul95)), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

mtext ("INLA SSF analysis. Wandering Albatross.", side = 3, outer = T, line = -0.2)
mtext ("6-hrly steps. Smooth terms for wind speed, wind support, and ab(cross wind) (binned data; no z-transofmration). no random effects.", 
       side = 3, outer = T, line = -2)
dev.off()


#--------------------------------------------------------------------------------------------------------------------------------------------------------

##################### only look at tracks with over 60 km/hr winds... doesnt change much. also try 70
#get a summary of number of tracks with over a certain wind speed
over_70 <- Waal_R1_windParam %>% 
  filter(wind_speed_kmh >= 70) %>% 
  summarize(ID = unique(TripID)) #149; was 219 with over 50

Waal_R1_windParam_70 <- Waal_R1_windParam[Waal_R1_windParam$TripID %in% over_70$ID,]

jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/alt_vs_used_70.jpeg", width = 15, height = 5, units = "in", res = 300)
X11(width = 13, height = 4)
par(mfrow= c(1,4), oma = c(0,0,3,0))
for(i in c("windDirection", "windSpeed_kmh", "cross_wind_kmh", "wind_support_kmh")){ #tried with non-interpolated wind support and crosswind and point of sail, but no difference
  boxplot(Waal_R1_windParam_70[,i] ~ Waal_R1_windParam_70[,"used"], 
          xaxt = "n", boxfill = c("gray","orange"), main = i, xlab = "", ylab = "")
  axis(1, labels = c("available", "used"), at = c(1,2), tck = 0)
}
mtext("Tracks containing 70 km/h winds (n = 33)", side = 3, outer = T, cex = 1)
dev.off()

#################### look at the mean wind speeds on each track. try mean over 40 (roughly the 3rd Qu.)
means_over_40 <- Waal_R1_windParam %>% 
  group_by(TripID) %>% 
  summarise(mean_ws = mean(wind_speed_kmh)) %>% 
  filter(mean_ws >= 40) #n = 18

Waal_R1_windParam_mean_40 <- Waal_R1_windParam[Waal_R1_windParam$TripID %in% means_over_40$TripID,]

jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/alt_vs_used_mean_40.jpeg", width = 15, height = 5, units = "in", res = 300)
X11(width = 13, height = 4)
par(mfrow= c(1,4), oma = c(0,0,3,0))
for(i in c("wave_h", "wind_speed_kmh", "cross_wind_kmh", "wind_support_kmh")){ #tried with non-interpolated wind support and crosswind and point of sail, but no difference
  boxplot(Waal_R1_windParam_mean_40[,i] ~ Waal_R1_windParam_mean_40[,"used"], 
          xaxt = "n", boxfill = c("gray","orange"), main = i, xlab = "", ylab = "")
  axis(1, labels = c("available", "used"), at = c(1,2), tck = 0)
}
mtext("Tracks with mean wind speed over 40 km/h (n = 18)", side = 3, outer = T, cex = 1)
dev.off()

X11(width = 15, height = 10);par(mfrow= c(3,1), oma = c(0,0,3,0))
for(i in c("wind_support_kmh", "cross_wind_kmh","wind_speed_kmh")){
  
  boxplot(Waal_R1_windParam_mean_40[,i] ~ Waal_R1_windParam_mean_40[,"FlyingSitting"], data = Waal_R1_windParam_mean_40, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")
  
  legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  
  boxplot(Waal_R1_windParam_mean_40[Waal_R1_windParam_mean_40$used == 1,i] ~ Waal_R1_windParam_mean_40[Waal_R1_windParam_mean_40$used == 1,"FlyingSitting"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(Waal_R1_windParam_mean_40$FlyingSitting)) - 0.15)
  boxplot(Waal_R1_windParam_mean_40[Waal_R1_windParam_mean_40$used == 0,i] ~ Waal_R1_windParam_mean_40[Waal_R1_windParam_mean_40$used == 0,"FlyingSitting"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(Waal_R1_windParam_mean_40$FlyingSitting)) + 0.15)
  
}

