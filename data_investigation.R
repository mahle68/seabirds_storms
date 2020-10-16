#investigating downloaded data for wind conditions
#Oct 12, 2020. Elham Nourani. Radolfzell, Germany

library(tidyverse)
library(sf)
library(lubridate)
library(sp)
library(mapview)

wgs <- CRS("+proj=longlat +datum=WGS84")

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

# David Gremiller's data
data_list <- list.files("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/data/David%20Gremillet/Scopoli_shearwater_GPS_20indiv_August2020")
