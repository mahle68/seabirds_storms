#extract dates and extent for the tracking data, to investigate wind conditions at a daily scale
#Elham Nourani. Radolfzell, Germany. Jan 20, 2021

library(readxl)
library(tidyverse)
library(lubridate)

#Start with the red-footed booby

#open data: prepared by Sophie (2012)
data_12 <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/RFBO_2012_ParamWind.csv", 
                 stringsAsFactors = F) %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))

unique_days_12 <- data_12 %>% 
  mutate(date = date(date_time)) %>% 
  distinct(date) %>% 
  arrange(date)

#extent (visual inspection in mapview)
extent <- data.frame(xmin = 35, xmax = 44, ymin= -27, ymax = -14)
 

#2015 data
juv_1 <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/boobies_2015/2015 EUROPA/CSV_OK/vague 1/juv01_RFBO_DZ26601_N14_EU2015.csv")  %>% 
  mutate(ID = "juv01_RFBO_DZ26601_N14_EU2015",
         date_time = as.POSIXct(strptime(paste(Date,Time," "),format = "%Y/%m/%d %H:%M:%S"),tz = "UTC")) %>% 
  dplyr::select(c("ID","date_time", "Latitude", "Longitude", "Altitude"))

files <- c(list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/boobies_2015/2015 EUROPA/CSV_OK/vague 1", pattern = ".xls",full.names = T),
           list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/boobies_2015/2015 EUROPA/CSV_OK/vague2", pattern = ".xls",full.names = T))

data_15 <- lapply(files, function(x){
  
  sheet_name <- substr(unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))],1,nchar(unlist(strsplit(x, "/"))[length(unlist(strsplit(x, "/")))])-5) #extract the name of the file, which will be the name of the sheet in excel
  x <- read_excel(x, sheet_name) %>% 
    dplyr::select(c(1,5:8))
  
  colnames(x)[c(1,2)] <- c("ID", "date_time")
  x
  
}) %>%
  reduce(rbind) %>% 
  full_join(juv_1)

save(data_15, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/RFB_2015.RData")


unique_days <- data_15 %>% 
  mutate(date = date(date_time)) %>% 
  distinct(date) %>% 
  full_join(unique_days_12) %>% 
  arrange(date) %>%
  data.frame()


save(unique_days, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/unique_days_RFB.RData")
  


