#investigating downloaded data for wind conditions 
#Oct 12, 2020. Elham Nourani. Radolfzell, Germany

library(tidyverse)
library(sf)
library(lubridate)
library(sp)
library(mapview)
#library(adehabitatHR)
library(concaveman) #concave hull 

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

mapview(data_vis,zcol = "common_name")


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

files <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/Movebank", full.names = T)

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

#data_vis <- data_flt %>% #data from 2008-2017 #dataset too large. crashes R (in laptop)
#  st_as_sf(coords = c("location.long","location.lat"),
#           crs = wgs) 

#mapview(data_vis,zcol = "individual.taxon.canonical.name")


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

MF <- data_flt %>% 
  filter(individual.taxon.canonical.name == "Fregata magnificens") %>% 
  mutate(year = year(timestamp),
         month = month(timestamp)) %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs)

MF_pol <- data_flt %>% 
  filter(individual.taxon.canonical.name == "Fregata magnificens") %>% 
  mutate(year = year(timestamp),
         month = month(timestamp)) %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
  st_transform(meters_proj) %>% 
  group_by(year, month) %>%
  filter(n() >= 2) %>%
  summarize(do_union = F) %>% 
  st_cast("LINESTRING") %>% 
  group_by(year, month) %>% 
  #rowwise() %>% 
  group_map(~ st_buffer(., 100000)) #took 10 min

MF_pol2 <- data_flt %>% 
  filter(individual.taxon.canonical.name == "Fregata magnificens") %>% 
  mutate(year = year(timestamp),
         month = month(timestamp)) %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
  st_transform(meters_proj) %>% 
  group_by(year, month) %>%
  filter(n() >= 2) %>%
  #summarize(do_union = F) %>% 
  #st_cast("LINESTRING") %>% 
  #group_by(year, month) %>% 
  #rowwise() %>% 
  group_map(~ st_buffer(., 500000)) %>% #way faster than converting to a line and calculating a buffer around the line
  map(~ st_union(.))



#%>% 
  
  st_union()
  #group_map(mcp)
  group_map(concaveman)
  

#coordinates(MF) <-~ location.long + location.lat
#proj4string(MF) <- wgs

MF_ls <- split(MF, list(MF$year, MF$month), drop = T)
MF_ls <- MF_ls[-which(lapply(MF_ls,nrow) < 3)] #remove groups with less than 3 points 


mcps <- lapply(MF_ls, function(x){
  pol <- x %>% 
    spTransform(meters_proj) %>% 
    mcp(x) #%>%  #calculate the mcp
  
}              
)


mcps <- lapply(MF_ls, function(x){
  pol <- x %>% 
    sp_transform(meters_proj) %>% 
    mcp(x) #%>%  #calculate the mcp

}              
)

#create move object to compare with the linestring
library(move)
mv <- move(x = st_coordinates(MF)[,"X"],y = st_coordinates(MF)[,"Y"],time = MF$timestamp,
           data = as.data.frame(MF,row.names = c(1:nrow(MF))),animal = MF$individual.local.identifier,proj = wgs)
