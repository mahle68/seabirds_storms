#script for finding data on cyclones, hurricanes and typhoons
#Sep 28, 2020. Elham Nourani. Konstanz, Germany

library(tidyverse)
library(lubridate)
library(sf)
library(sp)
library(mapview)

wgs <- CRS("+proj=longlat +datum=WGS84")

#read in the csv file downloaded from IBTrACS (https://www.ncdc.noaa.gov/ibtracs/index.php?name=ib-v4-access)
#meta-data available on https://www.ncdc.noaa.gov/ibtracs/pdf/IBTrACS_v04_column_documentation.pdf
storms <- read.csv("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/weather_data/ibtracs.ALL.list.v04r00.csv", stringsAsFactors = F)

#filter for category 5 hurricanes for 2017 and 2018

#extract names
  cat5_2017_18_names <- storms %>%
    filter(USA_SSHS == 5 & SEASON %in% c(2017, 2018)) %>% 
    select(NAME) %>% 
    distinct()

#extract entire hurricane tracks
cat5_2017_18 <- storms %>%
  filter(NAME %in% cat5_2017_18_names$NAME, SEASON %in% c(2017, 2018)) %>% 
  mutate_at(vars(9:10), as.numeric) %>% 
  st_as_sf(coords = c("LON", "LAT"),
           crs = wgs)

mapview(cat5_2017_18, zcol = "NAME")

#extract only parts of the track with category 3-5
cat5_2017_18 <- storms %>%
  filter(NAME %in% cat5_2017_18_names$NAME, SEASON %in% c(2017, 2018) & USA_SSHS %in% c(3,4,5)) %>% 
  st_as_sf(coords = c("LON", "LAT"),
           crs = wgs)

mapview(cat5_2017_18, zcol = "USA_SSHS")

#extract spatio-temporal extent

extents <- storms %>%
  filter(NAME %in% cat5_2017_18_names$NAME, SEASON %in% c(2017, 2018) & USA_SSHS %in% c(3,4,5)) %>% 
  mutate(timestamp = as.POSIXct(strptime(ISO_TIME,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  mutate_at(vars(9:10), as.numeric) %>% 
  group_by(NAME) %>%
  summarize(start = date(min(timestamp)),
            end = date(max(timestamp)),
            maxlat = max(LAT) + 5, #add 5-degree buffer to the extents
            minlat = min(LAT) - 5,
            maxlon = max(LON) + 5,
            minlon = min(LON) - 5)

write.csv(extents, "hurricanes_17_18.csv")
  
  

