# applying the step selection analysis to three species of seabirds
# We decided to do this for Peter Ryan's soft-plumaged and atlantic petrels and the albatross data (boxplots for these species showed some variation)
# Elham Nourani. Radolfzell, Germany. Dec. 10. 2020

#code partly from SSF_windy_tracks_albatrosses.R

library(tidyverse)
library(lubridate)
library(INLA)
library(move)
library(corrr)

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

# ---- STEP 1: open the data #####
# all have been filtered by Sophie to include trips with wind over 50 kmh

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/Peter_Ryan_data_annotated_SplitTrip_filterWind50kmh.RData") #PR_tot
alb <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/DATSETS/WAAL_allGPS_2010-2020_homo_R1h_TrackParam_Wind50kmh.csv",
                        na.strings = c("NA",""), fileEncoding="latin1") %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         scientific_name = "Diomedea exulans",
         common_name = "Wandering albatross") %>%
  rename(location.long = Longitude,
         location.lat = Latitude)
  #filter (FlyingSitting == "flying") %>% #filter out sitting positions
  #arrange(TripID,date_time)

#put together
cols <- intersect(names(alb), names(PR_tot))

data <- PR_tot[,cols] %>%
  filter(common_name %in% c("Atlantic Petrel", "Soft-plumaged Petrel")) %>%  #keep the species of interest
  full_join(alb[,c(cols, "BirdID")])

# ---- STEP 2: prepare alternative steps#####

move_ls<-lapply(split(data, data$common_name),function(x){
  x <- x %>%
    arrange(TripID, date_time) %>% 
    as.data.frame()
  mv <- move(x = x$location.long,y = x$location.lat,time = x$date_time,data = x,animal = x$TripID,proj = wgs)
  mv
})

#investage recording regime for each specie
str(lapply(move_ls, timeLag, units = "hours"))

hrs <- 4 #how long should the steps be? temporally
n <- 20 #how many alternative steps?

mycl <- makeCluster(detectCores() - 6) #6 cores, two for each species

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
      bursted_sp$TripID<-track@idData$TripID
      bursted_sp$species<-track@idData$common_name
    }
    
    bursted_sp$TripID<-track@idData$TripID 
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
        rnd_sp <- rnd
        coordinates(rnd_sp) <- ~lon+lat
        proj4string(rnd_sp) <- meters_proj
        rnd_sp <- spTransform(rnd_sp,wgs)
        
        #put used and available points together
        df <- used_point@data %>%  
          slice(rep(row_number(),n+1)) %>% #paste each row 20 times for the used and alternative steps
          mutate(location.long = c(head(location.long,1),rnd_sp@coords[,1]),
                 location.lat = c(head(location.lat,1),rnd_sp@coords[,2]),
                 used = c(1,rep(0,n)))  %>% #one hour after the start point of the step
          rowwise() %>% 
          mutate(heading = NCEP.loxodrome.na(lat1=current_point$location.lat,lat2=location.lat,lon1=current_point$location.long,lon2= location.long)) %>% 
          #dplyr::select(-c("u10m", "t2m", "press", "sst", "v10m", "X", "selected")) %>% 
          as.data.frame()
        
        #df[df$used == 0, c("wind_speed_ms", "wind_speed_kmh", "delta_t", "step_length", "turning_angle")] <- NA
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

save(used_av_ls, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/3spp_alt_steps_", hrs, "hr_", n, "n",".RData"))

#prepare to submit to Movebank
used_av_all <- lapply(used_av_ls, function(x){
  x %>% 
    mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% 
    as.data.frame()
}) %>% 
  reduce(rbind)
#row numbers are over a million, so do separate into two dfs for annotation
colnames(used_av_all)[c(1,2)] <- c("location-long","location-lat") #rename columns to match movebank format

write.csv(used_av_all, paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/3spp_alt_steps_", hrs, "hr_", n, "n",".csv"))

# ---- STEP 3: Exploration of annotated data#####
#open annotated data:

ann_cmpl <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/3spp_ssf/3spp_alt_steps_4hr_20n.csv-3788506526285067566.csv", stringsAsFactors = F) %>% 
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

# ---- STEP 4: Analysis #####

# ---- check for correlation
ann_cmpl %>% 
  dplyr::select(c("wind_speed_kmh","cross_wind_kmh","wind_support_kmh", "abs_cross_wind_kmh")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: var_cw with location.lat and var_delta_t with location.lat

# ---- prep data
ann_cmpl <- ann_cmpl %>%
  mutate(species1 = species, #repeat the grouping variable to use for random effects
         species2 = species) %>% 
  #z_transform
  mutate_at(c("wind_speed_kmh","cross_wind_kmh","abs_cross_wind_kmh","wind_support_kmh"),
            list(z = ~scale(.))) %>%
  #bin the data for smooth terms, otherwise I get an error that locations are too close.
  mutate_at(c("wind_speed_kmh","cross_wind_kmh","abs_cross_wind_kmh","wind_support_kmh"),
            list(group = ~inla.group(.,n = 50, method = "cut"))) %>%
  as.data.frame()

# ---- set mean and precision for the priors of slope coefficients (fixed effects)
mean.beta <- 0
prec.beta <- 1e-4 #precision of 1e-4 equals a variance of 1e4 ;)

# ---- set the formula (consider adding sl_ and ta_ later on)

#species as a random slope
formula <- used ~ -1 + 
  f(species, wind_speed_kmh_group, model = "rw2", constr = F) + 
  f(species1, wind_support_kmh_group, model = "rw2", constr = F) + 
  f(species2, abs_cross_wind_kmh_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T)))

#species as a fixed effect
formula2 <- used ~ -1 + species +
  f(wind_speed_kmh_group, model = "rw2", constr = F) + 
  f(wind_support_kmh_group, model = "rw2", constr = F) + 
  f(abs_cross_wind_kmh_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T)))

#only wind support. species as fixed effect
formula3 <- used ~ -1 + species +
  f(wind_support_kmh_group, model = "rw2", constr = F) + 
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T)))


# ---- model
  (b <- Sys.time())
  m <- inla(formula3, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = ann_cmpl,
            num.threads = 10,
            control.predictor = list(compute = T), #list(link = 1), #link is only relevant for NA observations. required to set the right link (i.e., the logit function) 
            #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T))
  Sys.time() - b #9.824472 hours for fomula; 8.905299 for formula2

save(m, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/3spp_wspt_full.RData") #n of 50 for binned wind


## ----label = "doseresp", fig = TRUE, echo = FALSE, fig.cap = '(ref:doseresp)'----
tab.rw1 <- data.frame(x = m$summary.random$dose[, "ID"],
                      y = m$summary.fitted.values[, "mean"], sex = H.virescens$sex)

ggplot(aes(x = dose, y = numdead / total), data = H.virescens) + 
  geom_point() +
  xlab(expression(paste("Dose (", mu, "g)"))) +
  ylab("Proportion of dead moths") +
  ggtitle("Dose response data (fitted values)") +
  geom_line(aes(x = x, y = y, linetype = "solid"), data = tab.rw1) +
  geom_line(aes(x = x, y = y, linetype = "dashed"), data = tab.rw2) +
  facet_grid(rows = vars(sex)) + 
  scale_linetype_manual(name = "Smoothing method",
                        values = c("solid", "dashed"),
                        labels = c("rw2", "rw1")) +
  theme(legend.position = "bottom")

# ---- plot
wspt <- data.frame(x = m$summary.random$wind_support_kmh_group[, "ID"],
                   y = m$summary.fitted.values[, "mean"],
                   ll95 = m$summary.fitted.values[,"0.025quant"],
                   ul95 = m$summary.fitted.values[,"0.975quant"],
                   species = ann_cmpl$species
)
# 
# wspd <- data.frame(x = m$summary.random$wind_speed_kmh_group[, "ID"],
#                    y = m$summary.random$wind_speed_kmh_group[, "mean"],
#                    ll95 = m$summary.random$wind_speed_kmh_group[,"0.025quant"],
#                    ul95 = m$summary.random$wind_speed_kmh_group[,"0.975quant"],
#                    species = ann_cmpl$species
# )
# 
# wspt <- data.frame(x = m$summary.random$wind_support_kmh_group[, "ID"],
#                    y = m$summary.random$wind_support_kmh_group[, "mean"],
#                    ll95 = m$summary.random$wind_support_kmh_group[,"0.025quant"],
#                    ul95 = m$summary.random$wind_support_kmh_group[,"0.975quant"]#,
#                    #species = ann_cmpl$species
# )
# 
# cw <- data.frame(x = m$summary.random$abs_cross_wind_kmh_group[, "ID"],
#                  y = m$summary.random$abs_cross_wind_kmh_group[, "mean"],
#                  ll95 = m$summary.random$abs_cross_wind_kmh_group[,"0.025quant"],
#                  ul95 = m$summary.random$abs_cross_wind_kmh_group[,"0.975quant"]#,
#                  #species = ann_cmpl$species
# )

#plot with raw value and credible intervals
jpeg(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/INLA_results_figures/", x$species[1], "_",hrs ,"hrs_all_smooth_0_1_abs_cw.jpeg"), width = 11, height = 3.5, units = "in", res = 300)
#X11(width = 11, height = 3.5)
par(mfrow= c(1,3), oma = c(0,1,3,0), bty = "l")

plot(x = ann_cmpl$wind_speed_kmh, y = ann_cmpl$used, pch = 16,  col = adjustcolor("grey", alpha.f = 0.1), xlab = "wind speed (kmh)", ylab = "exp(y)")
lines(wspd$x,exp(wspd$y)) 
polygon(x = c(wspd$x, rev(wspd$x)), y = c(exp(wspd$ll95),rev(exp(wspd$ul95))), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

plot(x = ann_cmpl$wind_support_kmh, y = ann_cmpl$used, pch = 16,  col = adjustcolor("grey", alpha.f = 0.1), xlab = "wind support (kmh)", ylab = "exp(y)")
lines(wspt$x,exp(wspt$y)) 
polygon(x = c(wspt$x, rev(wspt$x)), y = c(exp(wspt$ll95),rev(exp(wspt$ul95))), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

plot(x = ann_cmpl$abs_cross_wind_kmh, y = ann_cmpl$used, pch = 16,  col = adjustcolor("grey", alpha.f = 0.1), xlab = "abs(cross wind) (kmh)", ylab = "exp(y)")
lines(cw$x,exp(cw$y)) 
polygon(x = c(cw$x, rev(cw$x)), y = c(exp(cw$ll95),rev(exp(cw$ul95))), col = adjustcolor("grey", alpha.f = 0.3), border = NA)

mtext (paste0("INLA SSF analysis. ", x$species[1]), side = 3, outer = T, line = -0.2)
mtext (paste0(hrs,"-hrly steps. ", model_title), 
       side = 3, outer = T, line = -2)
dev.off()

