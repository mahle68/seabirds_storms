#create animated wind fields for the special cases
#Oct. 27. 2021



library(tidyverse)
library(lubridate)


setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")

#load data for significant strata
load("R_files/sig_data.RData") #sig_data
#tracking data
load("R_files/all_spp_df_colony_waal.RData") #data_df_all


#extract entire trips for significant steps
trips_df <- data_df_all %>% 
  filter(sci_name %in% sig_data$sci_name & TripID %in% sig_data$TripID)

trips_sf <- data_df_all %>% 
  filter(sci_name %in% sig_data$sci_name & TripID %in% sig_data$TripID) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)
  
mapview(trips_sf, zcol = "TripID")

# 
# sig_data$stratum <-as.character(sig_data$stratum)
# 
# avoidance_hour <- sig_data %>% 
#   rowwise() %>% 
#   mutate(unique_hour = paste(yday(timestamp), hour(timestamp), sep = "_"))

#unique timestamps for each trip
times <- trips_df %>% 
  mutate(year = year(timestamp),
         month = month(timestamp),
         yday = yday(timestamp),
         day = day(timestamp)) %>% 
  group_by(sci_name, TripID, year, month, yday) %>% 
  summarize(day = head(day,1)) %>% 
  ungroup()

days_ls <- lapply(split(times,times$TripID), function(x){ #one row for each day
  
  d <- list()
  
  for(i in unique(x$month)){
    min_day <- min(x[x$month == i, "day"])
    max_day <- max(x[x$month == i, "day"])
    
    days <- list(str_pad(min_day:max_day, 2,"left","0"))
    
    d <- c(d,days)
  }
  
  names(d) <- as.character(unique(x$month))
  
  d
})

extents <- trips_df %>% 
  group_by(sci_name, TripID) %>% 
  summarise(max.lat = max(location.lat) + 7, #order: N,W,S,E
            min.lon = min(location.long) - 20,
            min.lat = min(location.lat) - 7 ,
            max.lon = max(location.long) + 20)

# manually update extent for Sula sula
extents[extents$sci_name == "Sula sula", c(3:6)] <- list( 4, -94,-3, -86)

#--------------------request data download -----

#only for month and dates corresponding with the trip

#import the python library ecmwfapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()

output_path <- "/home/enourani/Documents/ERA_5_seabirds/"


lapply(c(1:nrow(extents[-6,])), function(x){ #skip atlantic yellow-nosed albatross (already done in seabird_complete)
  
  species <- as.character(extents[x, "sci_name"])
  trip <- as.character(extents[x, "TripID"])
  datetimes <- days_ls[names(days_ls) == trip]
  year <- as.character(unique(times[times$TripID == trip,"year"]))
  area <- as.numeric(extents[x,3:6])
  
  for (i in unique(names(datetimes[[1]]))){ #for each month
    month <- i
    days <- as.character(unlist(datetimes[[1]][names(datetimes[[1]]) == i]))
    
    request <- r_to_py(list(
      product_type = "reanalysis",
      variable = c("10m_u_component_of_wind", "10m_v_component_of_wind","sea_surface_temperature"), #sea surface temp helps later for cropping land
      year = year,
      month = month,
      day = days,
      time = str_c(seq(0,23,1),"00",sep = ":") %>% str_pad(5,"left","0"),
      area = area,  # North, West, South, East.
      grid = c(1, 1), #reduce the grid size
      format = "netcdf",
      dataset_short_name = "reanalysis-era5-single-levels"
    ))
    
    server$retrieve("reanalysis-era5-single-levels",
                    request,
                    target = paste0(output_path,"wind","_", strsplit(species, " ")[[1]][1], "_",trip, "_",month, ".nc")) 
    
  }
  
})

#--------------------extract date from netcdf files -----
setwd("/home/enourani/Documents/ERA_5_seabirds/")

file_list <- list.files(pattern = ".nc", full.names = TRUE)
vname <- c("u10","v10", "sst")


mycl <- makeCluster(detectCores() - 3, setup_strategy = "sequential") 
clusterExport(mycl, c("file_list", "vname")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(dplyr)
  library(ncdf4)
  library(lubridate)
})

(b <- Sys.time())

parLapply(mycl, file_list,function(x){
  
  nc <- nc_open(x)
  
  #extract lon and lat
  lat <- ncvar_get(nc,'latitude')
  nlat <- dim(lat) 
  lon <- ncvar_get(nc,'longitude')
  nlon <- dim(lon) 
  
  #extract the time
  t <- ncvar_get(nc, "time")
  nt <- dim(t)
  
  #convert the hours into date + hour
  timestamp <- as_datetime(c(t*60*60),origin = "1900-01-01")
  
  #put everything in a large df
  row_names <- expand.grid(lon,lat,timestamp)
  
  var_df <- data.frame(cbind(
    row_names,
    matrix(as.vector(ncvar_get(nc,vname[1])), nrow = nlon * nlat * nt, ncol = 1), #array to vector to matrix
    matrix(as.vector(ncvar_get(nc,vname[2])), nrow = nlon * nlat * nt, ncol = 1),
    matrix(as.vector(ncvar_get(nc,vname[3])), nrow = nlon * nlat * nt, ncol = 1)))
  
  colnames(var_df) <- c("lon","lat","date_time",vname)   #set column names
  
  #remove points over land (NAs)
  
  df <- var_df %>%
    #na.omit() %>% #let's keep points over land as well. otherwise, points where significant wave height is NA, i.e. land, will be deleted
    mutate(yday = yday(date_time),
           hour = hour(date_time),
           year = year(date_time)) %>%
    data.frame()
  
  save(df,file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/", 
                        unlist(strsplit(unlist(strsplit(x, ".nc")), "./"))[2], ".RData"))
})

Sys.time() - b #3 seconds :p

stopCluster(mycl) 


#--------------------append each trip's files together -----

#extract trip names

for (i in trips_df %>% filter(TripID != "73500_1") %>% distinct(TripID) %>% .$TripID){ #remove atlantic yellow nosed
  file_ls <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/",i,full.names = TRUE)
  
  data_df <- sapply(file_ls, function(x) mget(load(x)), simplify = TRUE) %>%
  reduce(rbind)
  
  save(data_df, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/", i, "_raw_wind.RData"))
}


#--------------------PLOT!!! -----
#https://semba-blog.netlify.app/10/29/2018/animating-oceanographic-data-in-r-with-ggplot2-and-gganimate/

world <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/continent.shp") %>% 
  #st_crop(xmin = -62, xmax = 11, ymin = -51, ymax = -25) %>%
  st_union()


lapply(c(trips_df %>% filter(TripID != "73500_1") %>% distinct(TripID) %>% .$TripID), function(x){
  
  files <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/", x, full.names = TRUE)
  
  wind_df <- sapply(files, function(x) mget(load(x)), simplify = TRUE) %>%
    reduce(rbind) %>% 
    mutate(wind_speed = sqrt(u10^2 + v10^2), #m/s 
           unique_hour = paste(yday, hour, sep = "_"))
  

  trip <- trips_df %>% 
    filter(TripID == x) %>% 
    mutate(unique_hour = paste(yday(timestamp), hour(timestamp), sep = "_")) %>% 
    mutate(avoidance = ifelse(unique_hour %in% avoidance_hour$unique_hour, "avoided", "not_avoided")) %>% 
    group_by(unique_hour) %>% 
    slice(1) #hourly filter

  region <- world %>% 
    st_crop(xmin = min(wind_df$lon), xmax = max(wind_df$lon), ymin = min(wind_df$lat), ymax = max(wind_df$lat))
    
    
  dir.create(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/animation/", 
                            head(trip$sci_name,1),"_", x, "/"))
  
  path <- paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/animation/", 
                            head(trip$sci_name,1),"_", x, "/")
  
    
  for(i in unique(trip$unique_hour)){
    
    plot <- ggplot() +
      geom_raster(data = wind_df %>% filter(unique_hour == i), aes(x = lon, y = lat, fill = wind_speed))+
      geom_segment(data = wind_df %>% filter(unique_hour == i), 
                   aes(x = lon, xend = lon+u10/10, y = lat, 
                       yend = lat+v10/10), arrow = arrow(length = unit(0.12, "cm")), size = 0.3)+
      geom_sf(data = region, fill = "grey85", col = 1)+
      geom_point(data = trip %>%  filter(unique_hour == i), aes(x = location.long, y = location.lat), 
                            size = 1, colour = "red") +
      #geom_point(data = trip %>%  filter(unique_hour == i), aes(x = location.long, y = location.lat, colour = as.factor(avoidance)), 
      #           size = 1) +
      #scale_color_manual(values = c("avoided" = "red", "not_avoided" = "forestgreen")) +
      coord_sf(xlim = range(wind_df$lon), ylim =  range(wind_df$lat))+
      scale_fill_gradientn(colours = oce::oceColorsPalette(120), limits = c(0,23), 
                           na.value = "white", name = "Speed\n (m/s)")+
      theme_bw()+
      theme(axis.text = element_text(size = 12, colour = 1),
            legend.text = element_text(size = 10, colour = 1), 
            legend.title = element_text(size = 12, colour = 1),
            legend.position = c(0.08,0.23),
            legend.background = element_rect(colour = 1, fill = "white"))+
      labs(x = NULL, y = NULL, title = wind_df %>%  filter(unique_hour == i) %>% .$date_time %>% .[1] %>% paste(head(trip$sci_name,1), x))
    
    ggsave(plot = plot, filename = paste0(path, x,"_",i,".jpeg"), 
           height = 6, width = 12, dpi = 300)
    
  }
})





load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/seabirds_storms/atlanitc_yellow_nosed_raw_wind.RData") #data_df

data_df <- data_df %>% 
  mutate(wind_speed = sqrt(u10^2 + v10^2)) #m/s  


ayl <- PR_data_split %>% 
  filter(common_name == "Atlantic Yellow-nosed Albatross", 
         TripID == sig_data[sig_data$common_name == "Atlantic Yellow-nosed Albatross", "TripID"]) %>% 
  rowwise() %>% 
  mutate(unique_hour = paste(yday(timestamp), hour(timestamp), sep = "_")) %>% 
  mutate(avoidance = ifelse(unique_hour %in% avoidance_hour$unique_hour, "avoided", "not_avoided"))



#example map: aggregate to lower spatial resolution

wind_lres <- data_df %>% 
  mutate(lat_lres = round(lat),
         lon_lres = round(lon)) %>% 
  group_by(date_time,lat_lres,lon_lres) %>% 
  summarise(u10 = mean(u10, na.rm = T),
            v10 = mean(v10, na.rm = T),
            sp = mean(sp, na.rm = T),
            swh = mean(swh, na.rm = T),
            wind_speed = mean(wind_speed),
            yday = head(yday,1),
            hour = head(hour,1)) %>% 
  rename(lat = lat_lres,
         lon = lon_lres) %>% 
  mutate(unique_hour = paste(yday,hour, sep = "_"))

save(wind_lres, file = "/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/R_files/seabirds_storms/atlanitc_yellow_nosed_raw_wind_lres.RData")

# 
# wind_day <- data_df %>% 
#   group_by(lon, lat, yday) %>% 
#   summarise(u10 = median(u10, na.rm = TRUE),
#             v10 = median(v10, na.rm = TRUE), 
#             wind_speed = median(wind_speed, na.rm = TRUE))


for(i in unique(ayl$unique_hour)){
  
  plot <- ggplot() +
    geom_raster(data = wind_lres %>% filter(unique_hour == i), aes(x = lon, y = lat, fill = wind_speed))+
    geom_segment(data = wind_lres %>% filter(unique_hour == i), 
                 aes(x = lon, xend = lon+u10/10, y = lat, 
                     yend = lat+v10/10), arrow = arrow(length = unit(0.12, "cm")), size = 0.3)+
    geom_sf(data = region, fill = "grey85", col = 1)+
    geom_point(data = ayl %>%  filter(unique_hour == i), aes(x = location.long, y = location.lat), 
               size = 1, colour = "red") +
    coord_sf(xlim = c(-62, 11), ylim =  c(-51, -25))+
    scale_fill_gradientn(colours = oce::oceColorsPalette(120), limits = c(0,23), 
                         na.value = "white", name = "Speed\n (m/s)")+
    theme_bw()+
    theme(axis.text = element_text(size = 12, colour = 1),
          legend.text = element_text(size = 10, colour = 1), 
          legend.title = element_text(size = 12, colour = 1),
          legend.position = c(0.08,0.23),
          legend.background = element_rect(colour = 1, fill = "white"))+
    labs(x = NULL, y = NULL, title = wind_lres %>%  filter(unique_hour == i) %>% .$date_time %>% .[1])
  
  ggsave(plot = plot, filename = paste0("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/animation/wind_field_",i,".jpeg"), 
         height = 6, width = 12, dpi = 300)
  
}

