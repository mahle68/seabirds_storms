#create animated wind fields for the special cases
#Oct. 27. 2021

library(tidyverse)
library(lubridate)
library(sf)
library(raster)
library(cowplot)
library(rcartocolor)
library(mapview)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(ggimage)
library(rworldmap)
library(moveVis)

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

#load data for significant strata
load("R_files/sig_data.RData") #sig_data
#tracking data
load("R_files/all_spp_df_colony_waal.RData") #data_df_all

#extract entire trips for significant steps
trips_df <- data_df_all %>% 
  filter(sci_name %in% sig_data$sci_name & TripID %in% sig_data$TripID)

save(trips_df, file = "R_files/trips_df_for_windfields.RData")

trips_sf <- data_df_all %>% 
  filter(sci_name %in% sig_data$sci_name & TripID %in% sig_data$TripID) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)
  
mapview(trips_sf, zcol = "TripID")


#extract unique avoidance hour for each trip 
avoidance_hour <- sig_data %>% 
  rowwise() %>% 
  mutate(unique_hour = paste(yday(timestamp), hour(timestamp), sep = "_")) %>% 
  group_by(TripID) %>% 
  summarize(unique_hour = unique(unique_hour))
   
save(avoidance_hour, file = "R_files/avoidance_hour_per_trip.RData")

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
  
  save(df,file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/", 
                        unlist(strsplit(unlist(strsplit(x, ".nc")), "./"))[2], ".RData"))
})

Sys.time() - b #3 seconds :p

stopCluster(mycl) 


#--------------------append each trip's files together -----

#extract trip names

for (i in trips_df %>% filter(TripID != "73500_1") %>% distinct(TripID) %>% .$TripID){ #remove atlantic yellow nosed
  file_ls <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/",i,full.names = TRUE)
  
  data_df <- sapply(file_ls, function(x) mget(load(x)), simplify = TRUE) %>%
  reduce(rbind)
  
  save(data_df, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/", i, "_raw_wind.RData"))
}


#--------------------PLOT entire trips!!! -----
#https://semba-blog.netlify.app/10/29/2018/animating-oceanographic-data-in-r-with-ggplot2-and-gganimate/

world <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
  #st_crop(xmin = -62, xmax = 11, ymin = -51, ymax = -25) %>%
  st_union()


lapply(c(trips_df %>% filter(TripID != "73500_1") %>% distinct(TripID) %>% .$TripID), function(x){
  
  files <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/", x, full.names = TRUE)
  
  wind_df <- sapply(files, function(x) mget(load(x)), simplify = TRUE) %>%
    reduce(rbind) %>% 
    mutate(wind_speed = sqrt(u10^2 + v10^2), #m/s 
           unique_hour = paste(yday, hour, sep = "_"))
  

  trip <- trips_df %>% 
    filter(TripID == x) %>% 
    mutate(unique_hour = paste(yday(timestamp), hour(timestamp), sep = "_")) %>% 
    #mutate(avoidance = ifelse(unique_hour %in% avoidance_hour$unique_hour, "avoided", "not_avoided")) %>% 
    group_by(unique_hour) %>% 
    slice(1) #hourly filter

  #save the entire trip data
  save(trip, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/processed/", 
                            head(trip$sci_name,1),"_", x, "_whole_trip.RData"))
  
  #save wind data for entire trip
  save(wind_df, file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/processed/", 
                           head(trip$sci_name,1),"_", x, "_whole_trip_wind.RData"))
  
  region <- world %>% 
    st_crop(xmin = min(wind_df$lon), xmax = max(wind_df$lon), ymin = min(wind_df$lat), ymax = max(wind_df$lat))
    
    
  dir.create(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/animation/", 
                            head(trip$sci_name,1),"_", x, "/"))
  
  path <- paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/animation/", 
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

save(wind_lres, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/seabirds_storms/atlanitc_yellow_nosed_raw_wind_lres.RData")

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
  
  ggsave(plot = plot, filename = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/animation/wind_field_",i,".jpeg"), 
         height = 6, width = 12, dpi = 300)
  
}

#--------------------PLOT for manuscript -----

#open GIS info
world <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/World_Continents.shp")
library(rworldmap)
worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id
world.df <- world.points[,c("long","lat","group", "region")]

#use common name instead of sci name for plotting
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/species_summary_data.RData") #lm_input

#open data
load("R_files/avoidance_hour_per_trip.RData") #avoidance_hour

#determine where the globe should be for each fig
inset_locations <- data.frame(TripID = c("2010_24","77884_1", "73500_1", "XGPS2015_2016_CRO_106_3_1.csv"),
                              inset_x = c(0.695, 0.725, 0.71, 0.701),
                              inset_y = 0.691)

#keep tracks that will be plotted
hrs_to_plot <- avoidance_hour %>% 
  filter(TripID %in% c("2010_24",
                       "77884_1",
                       "73500_1",
                       "XGPS2015_2016_CRO_106_3_1.csv") & 
           !(unique_hour) %in% c("48_4","326_20","326_21")) %>%  #get rid of second point for XGPS2015... and atlantic yellow-nosed
  ungroup() %>% 
  full_join(inset_locations)




wind_ls <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/processed/", pattern = "wind", full.names = T)
           
gps_ls <- setdiff(list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/processed/", full.names = T),wind_ls)

#load albatross silhouette
#alb <- raster("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/figures/albatross-bird-shape2.png")

#create a list of extent data frames
region_ls <- list(data.frame(x = c(35, 65, 65, 35), #2010_24
                          y = c(-47, -47, -27, -27)),
               data.frame(x = c(-61, -32, -32, -61), #73500_1
                          y = c(-45, -45, -26, -26)), 
               data.frame(x = c(-31, 4, 4, -31), #77884_1
                          y = c(-55, -55, -36 , -36)),
               data.frame(x = c(30, 72, 72, 30), #XGPS2015_2016_CRO_106_3_1
                          y = c(-69, -69, -52, -52))
               )
names(region_ls) <- hrs_to_plot$TripID


#read in the files and prep for plotting
plots <- lapply(hrs_to_plot$TripID, function(x){
  
  #below i tried to plot the used and available points, but they are too small compared to the plot extent
  #extract stratum data
  #stratum <- sig_data %>% 
   # filter(TripID == x)
  
    
  if(x != "73500_1"){ #do ayn later. it is not in the list of files created above
    
    load(grep(wind_ls, pattern = x, value = T)) #wind_df
    load(grep(gps_ls, pattern = x, value = T)) #trip
    
    wind_df <- wind_df %>% 
      filter(unique_hour == hrs_to_plot %>% filter(TripID == x) %>%  pull(unique_hour))
    
    point <- trip %>% 
      filter(unique_hour == hrs_to_plot %>% filter(TripID == x) %>%  pull(unique_hour)) %>% 
      mutate(image = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/figures/albatross-bird-shape3.png")

    trip <- trip %>% 
      arrange(timestamp)
      
    } else {
    
    load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/atlanitc_yellow_nosed_raw_wind_lres.RData") #wind_lres
    wind_df <- wind_lres %>% 
      mutate(unique_hour = paste(yday,hour, sep = "_")) %>% 
      filter(unique_hour == "326_23")
    
    load("R_files/trips_df_for_windfields.RData") #trips_df
    point <- trips_df %>% 
      filter(sci_name == "Thalassarche chlororhynchos") %>% 
      mutate(unique_hour = paste(yday(timestamp),hour(timestamp), sep = "_"),
             image = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/figures/albatross-bird-shape3.png") %>%  #path to albatross clip art
      filter(unique_hour == "326_23")
    
    trip <- trips_df %>% 
      filter(sci_name == "Thalassarche chlororhynchos") %>% 
      arrange(timestamp)
  }
    
    ext <- names(region_ls) %>%  str_detect(x) %>%  keep(region_ls, .) %>% flatten_dfr() #extract region for x. flatten is basically unlist
    
    #create inset map
    inset <- ggplot() + 
      geom_polygon(data = world.df, aes(x = long, y = lat, group = group), fill = "grey80") +
      geom_polygon(data = ext, aes(x = x, y = y), fill = NA, color = "black", size = 0.3) + 
      coord_map("ortho", orientation = c(point$location.lat+37, point$location.long, 0)) +
      scale_y_continuous(breaks = (-2:2) * 30) +
      scale_x_continuous(breaks = (-4:4) * 45) +
      geom_path(data = trip, aes(x = location.long, y = location.lat), size = 0.17) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            panel.background = element_rect(fill = "white")) +
      labs(x = NULL, y = NULL) 
    
    #create main map
    
    if(x == "73500_1"){ #add legend only to one plot
    main_map <- ggplot() +
      geom_tile(data = wind_df, aes(x = lon, y = lat, fill = wind_speed))+
      geom_segment(data = wind_df, 
                   aes(x = lon, xend = lon+u10/10, y = lat, 
                       yend = lat+v10/10), arrow = arrow(length = unit(0.12, "cm")), size = 0.3)+
      geom_sf(data = world, fill = "grey85")+
      #geom_path(data = trip, aes(x = location.long, y = location.lat), 
      #          size = 2, colour = "red") +
      geom_image(data = point,aes( x = location.long, y = location.lat, image = image),
                 size = 0.7, by = "width", nudge_y = 0.5, nudge_x = 0.3) +
      coord_sf(xlim = range(ext$x), ylim =  range(ext$y))+
      scale_fill_gradientn(colours = oce::oceColorsPalette(120), limits = c(0,23), 
                           na.value = "white", name = expression(paste("Wind speed (m s"^-1*")"))) +
      theme_bw()+
      theme(axis.text = element_text(size = 12, colour = 1),
            legend.text = element_text(size = 9, colour = 1), 
            legend.title = element_text(size = 10, colour = 1),
            legend.position = c(0.13,0.83),
            legend.spacing.y = unit(0, 'cm'),
            legend.background = element_blank(),
            plot.title = element_text(face = "italic"))+
      labs(x = NULL, y = NULL, title = paste(lm_input[lm_input$sci_name == point$sci_name, "species"], 
                                             paste(point$timestamp, "UTC", sep = " "), sep = "   "))
    
    } else {
      main_map <- ggplot() +
        geom_tile(data = wind_df, aes(x = lon, y = lat, fill = wind_speed), show.legend = F)+
        geom_segment(data = wind_df, 
                     aes(x = lon, xend = lon+u10/10, y = lat, 
                         yend = lat+v10/10), arrow = arrow(length = unit(0.12, "cm")), size = 0.3)+
        geom_sf(data = world, fill = "grey85")+
        #geom_path(data = trip, aes(x = location.long, y = location.lat), 
        #           size = 2, colour = "red") +
        geom_image(data = point,aes( x = location.long, y = location.lat, image = image),
                   size = 0.7, by = "width", nudge_y = 0.5, nudge_x = 0.3) +
        #geom_point(data = stratum, aes( x = location.long, y = location.lat, color = used), size = 1)+
        coord_sf(xlim = range(ext$x), ylim =  range(ext$y))+
        scale_fill_gradientn(colours = oce::oceColorsPalette(120), limits = c(0,23), 
                             na.value = "white", name = expression(paste("Wind speed (m s"^-1*")"))) +
        theme_bw()+
        theme(axis.text = element_text(size = 12, colour = 1),
              plot.title = element_text(face = "italic"))+
        labs(x = NULL, y = NULL, title = paste(lm_input[lm_input$sci_name == point$sci_name, "species"]
                                               , paste(point$timestamp, "UTC", sep = " "), sep = "   "))
    }

    
    #X11(height = 5, width = 7)
    final_plot <- ggdraw() +
      draw_plot(main_map) +
      draw_plot(inset, x = as.numeric(hrs_to_plot[hrs_to_plot$TripID == x, "inset_x"]),
                y = as.numeric(hrs_to_plot[hrs_to_plot$TripID == x, "inset_y"]), width = 0.25, height = 0.25) #set this manually for each panel
    
    #inset locations
    #for "2010_24": x = 0.695, y = 0.691
    #for "73500_1": x = 0.695, y = 0.691
    #for ""77884_1": x = 0.725, y = 0.691
    #for ""XGPS2015_2016_CRO_106_3_1.csv": 
    final_plot
    
    png(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/", x, ".png"), 
        height = 5, width = 7, units = "in", res = 300)
    print(final_plot)
    dev.off()
 
})

#stitch together in a very messy way

library(abind)
library(png)
img1 <- readPNG("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/2010_24.png")
img2 <- readPNG("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/73500_1.png")
img3 <- readPNG("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/77884_1.png")
img4 <- readPNG("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/XGPS2015_2016_CRO_106_3_1.csv.png")
two1 <- abind::abind(img1, img4, along = 1)
two2 <- abind::abind(img2, img3, along = 1)
#all <- abind::abind(two1, two2, along = 2)
png::writePNG(two1,"/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/two1.png")
png::writePNG(two2,"/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/two2.png")

#crop two1 and two2 in gimp
one <- readPNG("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/two1_narrow.png")
two <- readPNG("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/two2_narrow.png")

all <- abind::abind(one, two, along = 2)
png::writePNG(all,"/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22/wind_stills.png")

#--------------------ANIMATION for manuscript -----

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/species_summary_data.RData") #lm_input

#open world map
world <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/World_Continents.shp")
worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id
world.df <- world.points[,c("long","lat","group", "region")]

#open wind data. add column for unique hour
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/atlanitc_yellow_nosed_raw_wind_lres.RData") #wind_lres
wind_df <- wind_lres %>% 
  mutate(hr = str_pad(hour,2, side = "left", pad = "0")) %>% #make sure hour has two digits
  mutate(unique_hour = paste(yday,hr, sep = "_"))

#open tracking data for the trip of interest
load("R_files/trips_df_for_windfields.RData") #trips_df
point <- trips_df %>% 
  filter(sci_name == "Thalassarche chlororhynchos") %>% 
  mutate(hr = str_pad(hour(timestamp), 2, side = "left", pad = "0")) %>% #make sure hour has two digits
  mutate(unique_hour = paste(yday(timestamp),hr, sep = "_"),
         image = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/figures/albatross-bird-shape3.png") %>% #path to albatross clip art
  filter(unique_hour != "320_23") #the wind data doesnt have this unique hour, so remove it

trip <- trips_df %>% 
  filter(sci_name == "Thalassarche chlororhynchos") %>% 
  arrange(timestamp)

#extract the extent for the plot
ext <- wind_df %>% 
  st_as_sf(coords = c("lon", "lat"), crs = wgs) %>% 
  st_bbox()

ext_df <- data.frame(x = c(-58, 8, 8, -58),
                     y = c(-50, -50, -26, -26))


#create inset map
inset <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group), fill = "grey80") +
  geom_polygon(data = ext_df, aes(x = x, y = y), fill = NA, color = "black", size = 0.3) + 
  coord_map("ortho", orientation = c(point$location.lat[1]+37, point$location.long[1], 0)) +
  geom_path(data = trip, aes(x = location.long, y = location.lat), size = 0.17) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.background = element_rect(fill = "white")) +
  labs(x = NULL, y = NULL) 

#create main map

#loop over unique days
map_ls <- lapply(split(point, point$unique_hour), function(hr_point){

  unique_hr <- hr_point$unique_hour
  wind <-  wind_df[wind_df$unique_hour == unique_hr,]
  
  main_map <- ggplot() +
    geom_tile(data = wind, aes(x = lon, y = lat, fill = wind_speed))+
    geom_segment(data = wind, 
                 aes(x = lon, xend = lon+u10/10, y = lat, 
                     yend = lat+v10/10), arrow = arrow(length = unit(0.12, "cm")), size = 0.3)+
    geom_sf(data = world, fill = "grey85")+
    geom_image(data = hr_point, aes( x = location.long, y = location.lat, image = image), size = 0.4, asp = 1.3, by = "width", nudge_y = 0.5, nudge_x = 0.4) +
    #geom_point(data = hr_point, aes( x = location.long, y = location.lat, color = "red"), size = 3) + #plot the point to check the placement of albatross image
    coord_sf(xlim = range(ext_df$x), ylim = range(ext_df$y))+
    scale_fill_gradientn(colours = oce::oceColorsPalette(120), limits = c(0,23), 
                         na.value = "white", name = expression(paste("Wind speed (m s"^-1*")"))) +
    theme_bw()+
    theme(axis.text = element_text(size = 11, colour = 1),
          legend.key.width=unit(0.3,"cm"),
          legend.text = element_text(size = 8, colour = 1), 
          legend.title = element_text(size = 8, colour = 1),
          legend.position = c(0.06,0.835),
          legend.spacing.y = unit(0, 'cm'),
          legend.background = element_blank(),
          plot.title = element_text(face = "italic"))+
    labs(x = NULL, y = NULL, title = paste(lm_input[lm_input$sci_name == hr_point$sci_name, "species"], 
                                           paste(hr_point$timestamp, "UTC", sep = " "), sep = "   "))
  

#X11(height = 5, width = 10)
final_plot <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset, x = 0.805, y = 0.68, width = 0.25, height = 0.25) 

png(paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/animation_Jan23/frames/", unique_hr, ".png"), 
    height = 5, width = 10, units = "in", res = 500)
print(final_plot)
dev.off()

final_plot

})

#Animate in pop os terminal:
#ffmpeg -framerate 25 -pattern_type glob -i "*.png" output.mp4