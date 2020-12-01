#investigating downloaded data for wind conditions 
#Oct 12, 2020. Elham Nourani. Radolfzell, Germany

library(tidyverse)
library(sf)
library(lubridate)
library(sp)
library(mapview)


wgs <- CRS("+proj=longlat +datum=WGS84")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

# mcp <- function(x, percentile=95){ #mcp function that takes sf as input
# 
#   centroid <- sf::st_centroid(sf::st_union(x))
#   dist <- as.numeric(sf::st_distance(x, centroid))
#   within_percentile_range <- dist <= quantile(dist, percentile/100)
#   x_filter <- st_union(x[within_percentile_range,])
#   st_convex_hull(x_filter)
# 
# } #from package sdmflow

#----------------------------------------------------------------------------------
#Peter Ryan's data (requested and downloaded through the seabird tracking database)

files <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/Peter Ryan", full.names = T)

data <- lapply(files, read.csv,stringsAsFactors = F,
               colClasses = c("integer",rep("character",4),rep("numeric",2),"character","integer",rep("character",7),rep("numeric",2),"character")) %>% 
  reduce(full_join) %>% 
  filter(argos_quality %in% c("","0","1","2","3"))%>% 
  mutate(date_time = paste(date_gmt, time_gmt, " ")) %>%
  mutate(date_time,date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) 

data_vis <- data %>% #data from 2009-2014
  st_as_sf(coords = c("longitude","latitude"),
           crs = wgs) 

mapview(data_vis,zcol = "track_id")

#summarise for n individuals and tracks
data %>% 
  group_by(common_name) %>% 
  summarise(n_track = n_distinct(track_id))

#prep for movebank annotation request

data_mv <- data %>% 
  mutate(timestamp = paste(as.character(date_time),"000",sep = "."))

colnames(data_mv)[c(18,17)] <- c("location-long","location-lat")

write.csv(data_mv, "/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/annotation/Peter_Ryan_data.csv")   

#open annotated file
data_an <- read.csv("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/annotation/Peter Ryan/Peter_Ryan_data.csv-6108201444985501942.csv", stringsAsFactors = F) %>% 
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature ,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u10m = ECMWF.Interim.Full.Daily.SFC.Wind..10.m.above.Ground.U.Component.,
         v10m = ECMWF.Interim.Full.Daily.SFC.Wind..10.m.above.Ground.V.Component.,
         press = ECMWF.Interim.Full.Daily.SFC.Mean.Sea.Level.Pressure) %>% 
  mutate(wind_speed = sqrt(u10m^2 + v10m^2),
         delta_t = sst - t2m)

#calculate summaries of wind speed for each species (add to the spreadsheet)
data_an %>% 
  group_by(scientific_name) %>% 
  summarise(min_ws = min(wind_speed, na.rm = T),
            max_ws = max(wind_speed, na.rm = T),
            median_ws = median(wind_speed, na.rm = T),
            quant_95_ws = quantile(wind_speed, 0.95, na.rm = T))


#create histograms of wind speed

boxplot(data_an$wind_speed ~ data_an$common_name)

X11(width = 7, height = 9)
pdf("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/figures/Hist_Peter_Ryan_data.pdf", width = 7, height = 10)
par(mfrow= c(4,2), oma = c(0,0,3,0))
for(i in unique(data$common_name)){
 
    hist(data_an[data_an$common_name == i, "wind_speed"], ylab = "", xlab = "" ,main = i, xlim=c(0,25))
}
mtext("Frequency of wind speed values (m/s)- Peter Ryan's data", side = 3, outer = T, cex = 1)

dev.off()

#----------------------------------------------------------------------------------
# David Gremiller's data
data_list <- list.files("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/data/David%20Gremillet/Scopoli_shearwater_GPS_20indiv_August2020")

#----------------------------------------------------------------------------------
# Movebank data (public data; owners not contacted)


files <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/Movebank",pattern = ".csv", full.names = T) # studies where birds were manipulated (magneitc field manipulations and dispacement) were removed

data <- lapply(files, read.csv,stringsAsFactors = F) 

data_flt <- lapply(data, function(x){
  if("argos.lc" %in% names(x)){
    x <- x %>%
      filter(argos.lc %in% c("0","1","2","3"))
  } else {
    x <- x
  }
  
  x <- x %>% 
    mutate_at(c("individual.taxon.canonical.name", "tag.local.identifier", "individual.local.identifier"),as.character) %>% 
    dplyr::select(c("timestamp", "location.long", "location.lat", "sensor.type", "study.name" , "utm.northing" , "utm.easting"  , "utm.zone" ,
                    "individual.taxon.canonical.name" , "individual.local.identifier", "tag.local.identifier"))
}) %>% 
  reduce(full_join) %>%
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  drop_na()

save(data_flt, file = "data_flt.RData")


#convert to move objects
mv_ls <- lapply(split(data_flt,data_flt$individual.taxon.canonical.name),function(x){
  
  x <- x %>%
    arrange(individual.local.identifier,timestamp) %>% 
    as.data.frame()
  
  mv <- move(x = x$location.long, y = x$location.lat,time = x$timestamp, data = x, animal = x$individual.local.identifier, proj = wgs)
  mv
})


#create a data frame containing colony coordinates. I looked these up on the papers cited on the Movebank study page
colonies <- vector("list",8)
names(colonies) <- unique(data_flt$study.name)
colonies$`Chick-rearing movements of Yelkouan Shearwaters 2012-2014 (BirdLife Malta)`$x <- 14.37
  colonies$`Chick-rearing movements of Yelkouan Shearwaters 2012-2014 (BirdLife Malta)`$y <- 35.99
  colonies$`Foraging habitat of white-tailed tropicbirds (data from Santos et al. 2019)`$x <- -32.42
  colonies$`Foraging habitat of white-tailed tropicbirds (data from Santos et al. 2019)`$y <- -3.86
  colonies$`Free-ranging foraging trips of Manx shearwaters`$x <- 
    colonies$`Free-ranging foraging trips of Manx shearwaters`$y <- 

#try sophie's function for the Yelkouan shearwater
source("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/function_SplitTrips.R")
dd <- SplitTrips(data = data_flt[data_flt$individual.taxon.canonical.name == "Puffinus yelkouan",], col = NA, buffer = 5000, BirdID = 10) 
      


#summary stats for number of individuals and tracks
load ("data_flt.RData")

data_flt %>% 
  group_by(individual.taxon.canonical.name) %>% 
  summarise(n_ind = n_distinct(individual.local.identifier))

data_flt %>% 
  group_by(study.name) %>% 
  summarise(n_ind = n_distinct(individual.local.identifier))

data_vis <- data_flt %>% #data from 2008-2017 #dataset too large. crashes R (in laptop)
  st_as_sf(coords = c("location.long","location.lat"),
           crs = wgs) 

mapview(data_vis,zcol = "individual.taxon.canonical.name")
mapview(data_vis[data_vis$individual.taxon.canonical.name == "Puffinus yelkouan",], zcol = "individual.local.identifier")

#prep for movebank annotation request. over one million rows, so prep two files


data_mv <- data_flt %>% 
  mutate(timestamp = paste(as.character(timestamp),"000",sep = "."))

data_mv1 <- data_mv[1:(nrow(data_mv)/2),]
data_mv2 <- data_mv[((nrow(data_mv)/2)+1):nrow(data_mv),]

colnames(data_mv1)[c(2:3)] <- c("location-long","location-lat")
colnames(data_mv2)[c(2:3)] <- c("location-long","location-lat")

write.csv(data_mv1, "/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/annotation/movebank_public/mv_public_data1.csv")  
write.csv(data_mv2, "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/movebank_public/mv_public_data2.csv")  


#-- investigations of available wind ----
#use a sample: the Mag Frigatebird
#tried calculating an MCP around the points and putting a buffer around it (didnt find a good mcp function that takes sf as input), 
#then tried converting the points to a line and putting a buffer around the line (very time-consuming),
#the, tried to just directly put the puffer around each point, then unify them for each month-year combo (pretty fast!)

load("data_flt.RData")

#just the points for vis purposes
MF_raw <- data_flt %>% 
  filter(individual.taxon.canonical.name == "Fregata magnificens") %>% 
  mutate(yr = year(timestamp),
         mn = month(timestamp)) %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs)
# 
# MF_pol <- data_flt %>% 
#   filter(individual.taxon.canonical.name == "Fregata magnificens") %>% 
#   mutate(year = year(timestamp),
#          month = month(timestamp)) %>% 
#   st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
#   st_transform(meters_proj) %>% 
#   group_by(year, month) %>%
#   filter(n() >= 2) %>%
#   summarize(do_union = F) %>% 
#   st_cast("LINESTRING") %>% 
#   group_by(year, month) %>% 
#   #rowwise() %>% 
#   group_map(~ st_buffer(., 100000)) #took 10 min

#create polygons by putting 500 km buffers around the points and splitting by year and month
MF_pols <- data_flt %>% 
  filter(individual.taxon.canonical.name == "Fregata magnificens") %>% 
  mutate(yr = year(timestamp),
         mn = month(timestamp)) %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
  st_transform(meters_proj) %>% #need meters units for the buffer
  group_by(yr, mn) %>%
  filter(n() >= 2) %>%
  group_map(~ st_buffer(., 500000)) %>% #way faster than converting to a line and calculating a buffer around the line
  map(. %>% summarise(yr = head(year(timestamp),1), #summarise acts as st_union, but enables me to keep the data
                      mn = head(month(timestamp),1),
                  species = head(individual.taxon.canonical.name,1),
                  days_in_mn = days_in_month(head(timestamp,1))) %>% #this is useful when assigining timestamps to the random points in the next step 
        mutate(area = as.numeric(st_area(.)/1e+6))) %>%  #calc area in km2
  map(~ st_transform(., wgs))


#place random points on the polygon. numbers relative to the size of the polygon
#each ECMWF Era-Interim is 6400 km2 (80 km cell size)... let's say, one point per 20000 km2?

MF_pts <- lapply(MF_pols, function(x){
  n <- round(x$area/2e+4)
  pts_sf <- st_sample(x, size = n)  #sample n points over the polygon
  pts <- data.frame(lon = st_coordinates(pts_sf)[,1],
                    lat = st_coordinates(pts_sf)[,2],
                    species = x$species,
                    used = 0,
                    #assing date_time for the year-month combo
                    ##briefly looked at hour in tracking data and it seemed to have a rather uniform distribution, so sample from all available on ECMWF-ERA inerim
                    timestamp = as.POSIXct(strptime(paste(paste(x$yr,x$mn,sample(c(1:x$days_in_mn), n, replace = T), sep = "-"), 
                                                          paste(sample(c("00","06","12","18"), n, replace = T),"00","00",sep = ":"), sep = " "),
                                                    format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))
})

#addd available point to the data




#create move object to compare with the linestring
library(move)
mv <- move(x = st_coordinates(MF)[,"X"],y = st_coordinates(MF)[,"Y"],time = MF$timestamp,
           data = as.data.frame(MF,row.names = c(1:nrow(MF))),animal = MF$individual.local.identifier,proj = wgs)
