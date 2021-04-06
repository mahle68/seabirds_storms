#clean up and prep little auk data

library(tidyverse)
library(sp)
library(sf)
library(mapview)


  #the 2011 data are very sparse and the foraging trips are not complete. so don't use


wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/David Gremillet/GPS_little_auk/")

files_txt <- c(list.files("data 2014/", pattern = "raw_csv_format.txt", full.names = T),
               list.files("GPS tracks 2012/", pattern = "raw_csv_format.txt", full.names = T))
files_txt <- files_txt[-grep("suite",files_txt)] #remove files with suit in their names

data_df <- data.frame()

for(i in files_txt){
  data <- read.table(i, sep = ";")
  colnames(data) <- c("deviceID","day","month","year", "hour", "minute", "raw_lat", "raw_long", "parameter")
  data <- data %>%
    mutate(TrackID = gsub(".{19}$","", strsplit(i, "//")[[1]][2]),
           date_time = as.POSIXct(strptime(paste(paste(paste0("20",year),month,day, sep = "-"), paste(hour, minute, "00", sep = ":"), sep = " "),
                                           format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) 
  data_df <- rbind(data_df,data)
}


#as sophie about the crs
write.csv(data_df, "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/David Gremillet/GPS_little_auk/little_auk_2012_2014.csv")


#set the coordinates
data_sf <- data_df %>% 
  dplyr::select(-c("day","hour","minute")) %>% 
  rowwise() %>% 
  mutate(lat = gsub(".{1}$","",raw_lat),
         lon = gsub(".{1}$","",raw_long)) %>% 
  st_as_sf(coords = c("lon","lat"), crs = 32622) %>% 
  st_transform(wgs)




b <- read.table(files_txt[[2]], sep = ";")
colnames(b) <- c("device","day","month","year", "hour", "minute", "raw_lat", "raw_long", "parameter")

d <- readOGR("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/David Gremillet/GPS_little_auk/GPS tracks 2012/LIAK12EG22.kml")


#########################
library(sp)
library(maptools) 


proj <- GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]

st_crs(sf) = 4326

file <- ton.tableau.de.données
df = read.table(file, header=TRUE, sep="\t", dec=".")
coordinates(df) = c("Longitude", "Latitude")
fichier_sortie = paste(file,"shp",sep="")
writePointsShape(df, fichier_sortie)

write(proj,file=paste(fichier_sortie,".prj",sep=""))
