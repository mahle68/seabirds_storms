# extract daily means for unique days (following extract_days_range.R)
# code from Era_5_download_prep.R in delta_t project

#Era 5
library(tidyverse)
library(reticulate)
library(ncdf4)
library(sf)
library(lubridate)
library(parallel)

##### STEP 1: connect to cdsapi server #####

#import the python library ecmwfapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()

##### STEP 2: request and download data #####
#open unique days
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/unique_days_RFB.RData")
#2012 (all Oct and first of Nov)
year <- "2012"
month <- "10"
days <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18",
          "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31")

#2015 (all of Feb and first of Mar)
year <- "2015"
month <- "02"
days <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18",
          "19", "20", "21", "22", "23", "24", "25", "26", "27", "28")

##
path <- "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/ECMWF_daily/"
species <- "RFB"

query <- r_to_py(list(
  product_type = "reanalysis",
  format = "netcdf",
  variable = c("10m_u_component_of_wind", "10m_v_component_of_wind", "significant_height_of_combined_wind_waves_and_swell", "surface_pressure"),
  year = year,
  month = month,
  day = days,
  time = c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00",
           "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"),
  area = c(-14, 35, -27, 44),
  dataset = "reanalysis-era5-single-levels"
))

server$retrieve("reanalysis-era5-single-levels",
                query,
                target = paste0(path,species, "_", month, "_", year, ".nc"))


##### STEP 3: open and process data #####

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms")
species <- "RFB"

file_list <- list.files(path = "data/ECMWF_daily/", pattern = species, full.names = T)

#calculate temperature gradient
vname <- c("u10", "v10", "swh", "sp") #last two: wave height and surface pressure

daily_summaries <- lapply(file_list,function(x) {
  nc <- nc_open(x)
  
  #extract lon and lat
  lat <- ncvar_get(nc,'latitude')
  nlat <- dim(lat) 
  lon <- ncvar_get(nc,'longitude')
  nlon <- dim(lon) 
  
  #extract the time
  t <- ncvar_get(nc, "time")
  nt <- dim(t)
  
  #time unit: hours since 1900-01-01
  #ncatt_get(nc,'time')
  
  #convert the hours into date + hour
  timestamp <- as_datetime(c(t*60*60), origin="1900-01-01")
  
  row_names <- expand.grid(lon,lat,timestamp)
  
  #extract data for variables
  vars <- lapply(vname, function(i){
    assign(paste0("data_",i), matrix(as.vector(ncvar_get(nc,i)),
                                     nrow = nlon * nlat * nt, ncol = 1)) #an array coverted to a vector and then a matrix
  }) %>% 
    reduce(cbind)
  
  var_df <- data.frame(cbind(row_names, vars))  
  
  colnames(var_df) <- c("lon","lat","date_time",vname) #set column names
  
  #for each day, aggregate over all lat and long. take the mean and max. ALSO, only keep data for water(i.e remove cells where wave height is NA)
  var_df %>%
    filter(!is.na(swh)) %>% 
    mutate(day = day(date_time),
           year = year(date_time),
           month = month(date_time),
           wspd_kmh = sqrt(u10^2 + v10^2)*3.6) %>% 
    group_by(year, month, day) %>%
    summarise_at(c(colnames(var_df)[4:ncol(var_df)], "wspd_kmh"), list(avg = mean, max = max, min = min), na.rm = T) #all vars excpet lat, lon, time
  
}) %>% 
  reduce(rbind)

save(daily_summaries, file = "R_files/RFB_daily_summaries.RData")

##### STEP 4: plot and summarise data #####
#plot

daily_summaries <- daily_summaries %>% 
  mutate(timestamp = as.POSIXct(strptime(paste(paste(year,month,day, sep = "-"), "00:00:00", sep = " "),format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  arrange(timestamp)

#first try max wspd

X11()
par(mfrow = c(3,1),
    oma = c(3,0,0,0))

#wind speed
plot(NULL, xlim = c(0,63), ylim = c(0,80), labels = F, tck = 0, ann = F)
abline(h = 60, lty = 2, col = "gray50")
points(as.factor(daily_summaries$timestamp), daily_summaries$wspd_kmh_max, pch = 20, col= "tomato1", cex = 1.3)
points(as.factor(daily_summaries$timestamp), daily_summaries$wspd_kmh_avg, pch = 20, col= "goldenrod1", cex = 1.3)

axis(side = 1, at = unique(as.factor(daily_summaries$timestamp)), line = 0, labels = unique(as.factor(daily_summaries$timestamp)), 
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, las = 2)

axis(side = 2, at = seq(10,80, by = 10), line = 0, labels = seq(10,80, by = 10),
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, 
     las = 2)

mtext("Average wind speed (km/h)", 2, line = 2.5 ,las = 0, cex = 0.6, font = 3)

legend("topleft", legend = c("maximum", "average"), col = c("tomato1","goldenrod1"),
       pch = 20, cex = 0.8, pt.cex = 1.3, bg = "white", bty = "n") 
mtext("Wind speed at sea surfacce (Red_footed Booby)", 3, line = 0.4, cex = 1.1, font = 2)

#wave height
plot(NULL, xlim = c(0,63), ylim = c(0,6), labels = F, tck = 0, ann = F)
abline(h = 5, lty = 2, col = "gray50")
points(as.factor(daily_summaries$timestamp), daily_summaries$swh_max, pch = 20, col= "springgreen4", cex = 1.3)
points(as.factor(daily_summaries$timestamp), daily_summaries$swh_avg, pch = 20, col= "yellowgreen", cex = 1.3)

axis(side = 1, at = unique(as.factor(daily_summaries$timestamp)), line = 0, labels = unique(as.factor(daily_summaries$timestamp)), 
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, las = 2)

axis(side = 2, at = seq(1,5, by = 1), line = 0, labels = seq(1,5, by = 1),
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, 
     las = 2)

mtext("Wave height (m)", 2, line = 2.5 ,las = 0, cex = 0.6, font = 3)

legend("topleft", legend = c("maximum", "average"), col = c("springgreen4","yellowgreen"),
       pch = 20, cex = 0.8, pt.cex = 1.3, bg = "white", bty = "n") 

#pressure
plot(NULL, xlim = c(0,63), ylim = c(89000,100800), labels = F, tck = 0, ann = F)
points(as.factor(daily_summaries$timestamp), daily_summaries$sp_avg, pch = 20, col= "slateblue3", cex = 1.3)
points(as.factor(daily_summaries$timestamp), daily_summaries$sp_min, pch = 20, col= "thistle1", cex = 1.3)

axis(side = 1, at = unique(as.factor(daily_summaries$timestamp)), line = 0, labels = unique(as.factor(daily_summaries$timestamp)), 
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, las = 2)

axis(side = 2, at = seq(1,5, by = 1), line = 0, labels = seq(1,5, by = 1),
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, 
     las = 2)

mtext("Surfac pressure", 2, line = 2.5 ,las = 0, cex = 0.6, font = 3)

legend("topleft", legend = c("maximum", "min"), col = c("thistle1","slateblue3"),
       pch = 20, cex = 0.8, pt.cex = 1.3, bg = "white", bty = "n") 

#----------------------------------------------------
#overall average over the years

output_df <- do.call("rbind", data_list)%>% 
  group_by(lat,lon)%>%
  summarise(avg_yr_u=mean(avg_u),
            avg_yr_v=mean(avg_v))  #average u and v over all years

#var_avg<-aggregate(var_area,by=list(c("lat","lon"))) #same as above



coordinates(output_df)<-~lon+lat

gridded(output_df)<-TRUE
avg_rstr<-stack(output_df) #create a raster stack

windows();plot(avg_rstr) #the average u and v over 2011-2014

#save average raster files

writeRaster(avg_rstr, filename=c("avg_wind_u.tif","avg_wind_v.tif"),bylayer=TRUE, overwrite=FALSE)
































#####era-interim
#library(reticulate)
library(ncdf4)
library(ecmwfr)

#download ecmwf data using ecmwfr
#https://cran.r-project.org/web/packages/ecmwfr/vignettes/webapi_vignette.html

wf_set_key(user = "mahle68@gmail.com",
           key = "0335373f68bfb0f54f60a8679840f302",
           service = "webapi")

year <- 2012
target <- paste0(year,"_.nc")
path <- "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/ECMWF_daily/"

#wf_get_key(user = "mahle68@gmail.com",
#           service= "webapi")

#copy Mars request and then from addins, choose mars to list
wind_and_pressure_request <- list(
  class   = "ei",
  dataset = "interim",
  date    = "2012-10-01/to/2012-10-31",
  expver  = "1",
  grid    = "0.75/0.75",
  levtype = "sfc",
  param   = "134.128/151.128/165.128/166.128",
  step    = "0",
  stream  = "oper",
  time    = "00:00:00/06:00:00/12:00:00/18:00:00",
  type    = "an",
  format  = "netcdf",
  target  = paste0("wind_press_",target)
)

wave_request <- list(
  class   = "ei",
  dataset = "interim",
  date    = "2012-10-01/to/2012-10-31",
  expver  = "1",
  grid    = "0.75/0.75",
  levtype = "sfc",
  param   = "229.140",
  step    = "0",
  stream  = "wave",
  time    = "00:00:00/06:00:00/12:00:00/18:00:00",
  type    = "an",
  format  = "netcdf",
  target  = paste0("wave_",target)
)


for (request in list(wind_and_pressure_request, wave_request)){
  
  wf_request(
  user = "mahle68@gmail.com",
  request = request,
  transfer = T,
  path = path
)
  
}



