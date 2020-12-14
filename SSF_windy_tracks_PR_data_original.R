#script for exploration of seabird data to find their response to strong wind conditions
#Elham Nourani, PhD. Nov. 23. 2020. Radolfzell am Bodensee, Germany
#tried filtering for tracks that contain wind speeds equal to or higher than the thrid quartile. turns out all tracks do. 
#try filtering for tracks with a mean wind speed equal to or greater than the third quartile for the species.

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

# STEP 1: identify windy tracks ##### 

#open annotated data
data_an <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/Peter Ryan/Peter_Ryan_data.csv-6108201444985501942.csv", stringsAsFactors = F) %>% 
  filter(device == "GPS") %>%  #only use GPS sensors for even sampling frequency
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature ,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u10m = ECMWF.Interim.Full.Daily.SFC.Wind..10.m.above.Ground.U.Component.,
         v10m = ECMWF.Interim.Full.Daily.SFC.Wind..10.m.above.Ground.V.Component.,
         press = ECMWF.Interim.Full.Daily.SFC.Mean.Sea.Level.Pressure) %>% 
  mutate(timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         wind_speed_ms = sqrt(u10m^2 + v10m^2),
         wind_speed_kmh = sqrt(u10m^2 + v10m^2)*3.6,
         delta_t = sst - t2m)

#remove duplicated timestamps
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(data_an$track_id),timestamps = data_an$timestamp,
                                                        sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows

data_an <- data_an[-rows_to_delete,] 

#calculate summaries of wind speed for each species (add to the spreadsheet)
summaries <- data_an %>% 
  group_by(scientific_name) %>% 
  summarise(min_ws = min(wind_speed_kmh, na.rm = T),
            median_ws = median(wind_speed_kmh, na.rm = T),
            quant_75_ws = quantile(wind_speed_kmh, 0.75, na.rm = T),
            max_ws = max(wind_speed_kmh, na.rm = T),)


windy_tracks <- lapply(split(data_an, data_an$scientific_name), function(x){ #take the tracks that contain wind speeds equal to or stronger than the 3rd quartile
  q3 <- summaries[summaries$scientific_name == unique(x$scientific_name), "quant_75_ws"] #extract value for third quartile
  #mx <- summaries[summaries$scientific_name == unique(x$scientific_name), "max_ws"] #extract value for third quartile
  
  #track_means <- x %>% 
  #  group_by(track_id) %>% 
  #  summarise(mean_ws = mean(wind_speed_kmh, na.rm = T))
  
  w_tracks <- x %>% 
    filter(wind_speed_kmh >= q3$quant_75_ws) %>% 
    summarise(tracks = unique(track_id))
  
  x %>% 
    filter(track_id %in% w_tracks$tracks) #only keep tracks that contain high wind speeds 
  
})

#Just use the whole data, since all tracks contain winds higher than the 3rd quartile

# STEP 2: prepare alternative steps#####

move_ls<-lapply(windy_tracks,function(x){
  
  x <- x %>%
    arrange(track_id, timestamp) %>% 
    as.data.frame()
  
  mv <- move(x = x$location.long,y = x$location.lat,time = x$timestamp,data = x,animal = x$track_id,proj = wgs)
  mv
})

#investage recording regime for each specie
str(lapply(move_ls, timeLag, units = "hours"))

hrs <- 6 #how long should the steps be? temporally
n <- 20 #how many alternative steps?

mycl <- makeCluster(detectCores() - 5) #7 cores, one for each species

clusterExport(mycl, c("move_ls", "hrs", "n",  "wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(dplyr)
  library(purrr)
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

used_av_ls <- parLapply(cl = mycl, X = move_ls,fun = function(species){ 
  
  #used_av_ls <- lapply(move_ls, function(species){ 
  sp_obj_ls <- lapply(split(species),function(track){ #sp_obj_ls will have the filtered and bursted trackments
    
    #--STEP 1: thin the data to n-hourly intervals
    track_th <- track %>%
      thinTrackTime(interval = as.difftime(hrs, units='hours'),
                    tolerance = as.difftime(15, units='mins')) #the unselected bursts are the large gaps between the selected ones
    #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst... longer gaps will divide the brusts) 
    track_th$selected <- c(as.character(track_th@burstId),NA) #assign selected as a variable
    track_th$burst_id <-c(1,rep(NA,nrow(track_th)-1)) #define value for first row
    
    if(nrow(track_th@data) == 1){
      track_th@data$burst_id <- track_th$burst_id
    } else {for(i in 2:nrow(track_th@data)){
      
      if(i== nrow(track_th@data)){
        track_th@data$burst_id[i]<-NA
      } else
        if(track_th@data[i-1,"selected"] == "selected"){
          track_th@data$burst_id[i]<-track_th@data[i-1,"burst_id"]
        } else {
          track_th@data$burst_id[i]<-track_th@data[i-1,"burst_id"]+1
        }
    }
    }
    #convert back to a move object (from move burst)
    track_th <- as(track_th,"Move")
    
    #--STEP 3: calculate step lengths and turning angles 
    #sl_ and ta_ calculations should be done for each burst. converting to a move burst doesnt make this automatic. so just split manually
    burst_ls<-split(track_th,track_th$burst_id)
    burst_ls<-Filter(function(x) length(x) >= 3, burst_ls) #remove bursts with less than 3 observations
    
    burst_ls<-lapply(burst_ls,function(burst){
      burst$step_length<-c(distance(burst),NA) #
      burst$turning_angle<-c(NA,turnAngleGc(burst),NA)
      burst
    })
    
    #put burst_ls into one dataframe
    bursted_sp <- do.call(rbind,burst_ls) 
    
    #reassign values
    
    if(length(bursted_sp) >= 1){
      bursted_sp$track<-track@idData$track_id
      bursted_sp$species<-track@idData$common_name
    }
    
    bursted_sp$track_id<-track@idData$track_id 
    bursted_sp
  }) %>% 
    Filter(function(x) length(x) > 1, .) #remove tracks with no observation (these have only one obs due to the assignment of trackment id)
  
    #--STEP 4: estimate step length and turning angle distributions
  #put everything in one df
  bursted_df <- sp_obj_ls %>%  
    reduce(rbind) %>% 
    as.data.frame() %>%
    filter(step_length < 500000) %>% #filter out points with over 500 km step length. outliers 
    dplyr::select(-c("coords.x1","coords.x2"))
  
  #estimate von Mises parameters for turning angles
  #calculate the averages (mu).steps: 1)convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan.OR use circular::mean.circular
  mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
  sl<-bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  #plot
  #jpeg(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/ta_sl_dist_2hr_", bursted_df$species[1],".jpeg"))
  jpeg(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/ta_sl_dist_",bursted_df$species[1], "_", hrs, "hr_", n, "n",".jpeg"))
  #X11()
  par(mfrow=c(1,2))
  hist(sl,freq=F,main="",xlab = "Step length (km)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  
  hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]), freq=F, main="",xlab="Turning angles (radians)")
  plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
  #mtext(paste0("Step length and turning angle distributions (2-hr) ",bursted_df$species[1]), side = 3, outer =T,line = -3)
  mtext(paste0("Step length and turning angle distributions (", hrs, "-hr) ", bursted_df$species[1]), side = 3, outer =T,line = -3, font = 2) 
  dev.off()
  
  
  #--STEP 5: produce alternative steps
  used_av_track <- lapply(sp_obj_ls, function(track){ #for each track
    
    used_av_burst <- lapply(split(track,track$burst_id),function(burst){ #for each burst,
      
      #assign unique step id
      burst$step_id <- 1:nrow(burst)
      
      used_av_step <- lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
        
        current_point<- burst[this_point,]
        previous_point<-burst[this_point-1,] #this is the previous point, for calculating turning angle.
        used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
        
        #randomly generate 20 step lengths and turning angles
        rta <- as.vector(rvonmises(n = n, mu = mu, kappa = kappa)) #generate random turning angles with von mises distribution (in radians)
        rsl <- rgamma(n = n, shape = fit.gamma1$estimate[[1]], rate = fit.gamma1$estimate[[2]])*1000  #generate random step lengths from the gamma distribution. make sure unit is meters
        
        #calculate bearing of previous point
        #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
        prev_bearing <- NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
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
          mutate(location.long = c(head(location.long,1),rnd_sp@coords[,1]),
                 location.lat = c(head(location.lat,1),rnd_sp@coords[,2]),
                 used = c(1,rep(0,n)))  %>% #one hour after the start point of the step
          rowwise() %>% 
          mutate(heading = NCEP.loxodrome.na(lat1=current_point$location.lat,lat2=location.lat,lon1=current_point$location.long,lon2= location.long)) %>% 
          dplyr::select(-c("u10m", "t2m", "press", "sst", "v10m", "X", "selected")) %>% 
          as.data.frame()
        
        df[df$used == 0, c("wind_speed_ms", "wind_speed_kmh", "delta_t", "step_length", "turning_angle")] <- NA
        df
        
      }) %>% 
        reduce(rbind)
      used_av_step
    }) %>% 
      reduce(rbind)
    used_av_burst
  }) %>% 
    reduce(rbind)
  used_av_track
})

Sys.time() - start_time #22 min

stopCluster(mycl)

save(used_av_ls, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/PR_alt_steps_", hrs, "hr_", n, "n",".RData"))

#prepare to submit to Movebank
used_av_all <- lapply(used_av_ls, function(x){
  x %>% 
    mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) %>% 
    as.data.frame()
}) %>% 
  reduce(rbind)

#have a look
X11();par(mfrow= c(1,1), mar = c(0,0,0,0), oma = c(0,0,0,0))
maps::map("world",fil = TRUE,col = "grey85", border=NA, ylim = c(-65,-10), xlim = c(-60,30)) 
points(used_av_all[used_av_all$used == 0,c("location.long","location.lat")], pch = 16, cex = 0.2, col = "gray55")
points(used_av_all[used_av_all$used == 1,c("location.long","location.lat")], pch = 16, cex = 0.2, col = "orange")


#row numbers are over a million, so do separate into two dfs for annotation
colnames(used_av_all)[c(4,3)] <- c("location-long","location-lat") #rename columns to match movebank format

write.csv(used_av_all, paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/PR_alt_steps_", hrs, "hr_", n, "n",".csv"))


# STEP 3: Exploration of annotated data#####
#open annotated data:
#open annotated data and add wind support and crosswind

ann_cmpl <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/PR_6hr/PR_alt_steps_6hr_20n.csv-6356372864490415351.csv", stringsAsFactors = F) %>% 
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
         stratum = paste(track_id, burst_id, step_id, sep = "_")) %>% 
  mutate(wind_support_kmh = wind_support_ms * 3.6,
         cross_wind_kmh = cross_wind_ms * 3.6,
         wind_speed_kmh = wind_speed_ms * 3.6) %>% 
  as.data.frame()

save(ann_cmpl, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/PR_data_ssf_ann_", hrs, "hr_", n, "n",".RData"))

load(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/PR_data_ssf_ann_", hrs, "hr_", n, "n",".RData"))

#plot
jpeg("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/Used_avail_PR_6hr.jpeg", width = 12, height = 10, units = "in", res = 300)
     
X11(width = 12, height = 10)
par(mfrow= c(3,1), oma = c(0,0,3,0))
for(i in c("wind_support_kmh", "cross_wind_kmh","wind_speed_kmh")){
  
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")
  #if(i == "wind_support_kmh"){
    legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  #}
  boxplot(ann_cmpl[ann_cmpl$used == 1,i] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0,i] ~ ann_cmpl[ann_cmpl$used == 0,"species"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  #  } 
}
mtext(paste0("Used and available wind conditions (", hrs, " hr)"), side = 3, outer = T, cex = 1.3)
dev.off()


# STEP 4: Analysis #####

#only keep data for soft-plumaged and atlantic petrels
ann_cmpl <- ann_cmpl %>% 
  filter(species %in% c("Atlantic Petrel", "Soft-plumaged Petrel"))

#correlation
ann_cmpl %>% 
  dplyr::select(c("wind_speed_kmh","cross_wind_kmh","wind_support_kmh")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: var_cw with location.lat and var_delta_t with location.lat

#model each species separately

# Set mean and precision for the priors of slope coefficients (fixed effects)
mean.beta <- 0
prec.beta <- 1e-4 #precision of 1e-4 equals a variance of 1e4 ;)

formula <- used ~ -1 + f(wind_speed_kmh_group, model = "rw2", constr = F) + 
  f(wind_support_kmh_group, model = "rw2", constr = F) + 
  f(abs_cw_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T)))

formula2 <- used ~ -1 + 
  f(wind_speed_kmh_group, model = "rw2", hyper = list(
    theta = list(prior = "pc.prec", param = c(0, 0.03))), constr = F) + 
  f(wind_support_kmh_group, model = "rw2", hyper = list(
    theta = list(prior = "pc.prec", param = c(0, 0.03))), constr = F) + 
  f(abs_cw_group, model = "rw2", hyper = list(
    theta = list(prior = "pc.prec", param = c(0, 0.03))), constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6), fixed = T)))

model_title <- "Smooth terms for wind speed, wind support, and ab(cross wind) (binned data; no z-transofmration). no random effects."

models <- lapply(split(ann_cmpl, ann_cmpl$species), function(x){
  
  #prep data
  x <- x %>% 
    mutate(abs_cw = abs(cross_wind_kmh)) %>% 
    #z_transform
    mutate_at(c("wind_speed_kmh","cross_wind_kmh","abs_cw","wind_support_kmh"),
              list(z = ~scale(.))) %>%
    #bin the data for smooth terms, otherwise I get an error that locations are too close.
    mutate_at(c("wind_speed_kmh","cross_wind_kmh","abs_cw","wind_support_kmh"),
              list(group = ~inla.group(.,n = 50, method = "cut"))) %>%
    as.data.frame()
  
  #model
  (b <- Sys.time())
  m <- inla(formula, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = x,
            num.threads = 10,
            control.predictor = list(compute = T), #list(link = 1), #link is only relevant for NA observations. required to set the right link (i.e., the logit function) 
            #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))
  Sys.time() - b 
  
  save(m, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/", x$species[1], "_m4_full.RData")) #n of 50 for binned wind
  
  #plot
  wspd <- data.frame(x = m$summary.random$wind_speed_kmh_group[, "ID"],
                     y = m$summary.random$wind_speed_kmh_group[, "mode"],
                     ll95 = m$summary.random$wind_speed_kmh_group[,"0.025quant"],
                     ul95 = m$summary.random$wind_speed_kmh_group[,"0.975quant"]
  )
  
  wspt <- data.frame(x = m$summary.random$wind_support_kmh_group[, "ID"],
                     y = m$summary.random$wind_support_kmh_group[, "mean"],
                     ll95 = m$summary.random$wind_support_kmh_group[,"0.025quant"],
                     ul95 = m$summary.random$wind_support_kmh_group[,"0.975quant"]
  )
  
  cw <- data.frame(x = m$summary.random$abs_cw_group[, "ID"],
                   y = m$summary.random$abs_cw_group[, "mean"],
                   ll95 = m$summary.random$abs_cw_group[,"0.025quant"],
                   ul95 = m$summary.random$abs_cw_group[,"0.975quant"]
  )
  
  #plot with raw value and credible intervals
  jpeg(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/INLA_results_figures/", x$species[1], "_",hrs ,"hrs_all_smooth_0_1_abs_cw.jpeg"), width = 11, height = 3.5, units = "in", res = 300)
  #X11(width = 11, height = 3.5)
  par(mfrow= c(1,3), oma = c(0,1,3,0), bty = "l")
  
  plot(x = x$wind_speed_kmh, y = x$used, pch = 16,  col = adjustcolor("grey", alpha.f = 0.1), xlab = "wind speed (kmh)", ylab = "exp(y)")
  lines(wspd$x,exp(wspd$y)) 
  polygon(x = c(wspd$x, rev(wspd$x)), y = c(exp(wspd$ll95),rev(exp(wspd$ul95))), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
  
  plot(x = x$wind_support_kmh, y = x$used, pch = 16,  col = adjustcolor("grey", alpha.f = 0.1), xlab = "wind support (kmh)", ylab = "exp(y)")
  lines(wspt$x,exp(wspt$y)) 
  polygon(x = c(wspt$x, rev(wspt$x)), y = c(exp(wspt$ll95),rev(exp(wspt$ul95))), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
  
  plot(x = x$abs_cw, y = x$used, pch = 16,  col = adjustcolor("grey", alpha.f = 0.1), xlab = "abs(cross wind) (kmh)", ylab = "exp(y)")
  lines(cw$x,exp(cw$y)) 
  polygon(x = c(cw$x, rev(cw$x)), y = c(exp(cw$ll95),rev(exp(cw$ul95))), col = adjustcolor("grey", alpha.f = 0.3), border = NA)
  
  mtext (paste0("INLA SSF analysis. ", x$species[1]), side = 3, outer = T, line = -0.2)
  mtext (paste0(hrs,"-hrly steps. ", model_title), 
         side = 3, outer = T, line = -2)
  dev.off()
})

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/Soft-plumaged Petrel_m4_full.RData")
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/Atlantic Petrel_m4_full.RData")
