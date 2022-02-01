#script for analysis and figure prep for the seabird manuscript
#Sep 30. 2021. Elham Nourani, PhD. Radolfzell, DE.
#sources: null_modeling_seabirds.R and PCoA_seabirds.R

library(tidyverse)
library(parallel)
library(hrbrthemes)
library(ggridges)
library(corrr)
library(sp)
library(sf)
library(mapview)
library(lubridate)
library(reticulate)
library(ncdf4)
library(oce)
library(gridExtra)
library(ape) #for phylogeny
#library(gganimate)
#install.packages('devtools')
#devtools::install_github('mpio-be/windR')
#library(windR)

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")
source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

### STEP 1: Open and prep SSF data ------------------------------------ #####
#open dataset with alternative steps for hourly steps. prepped in random_steps.R
load("R_files/ssf_input_annotated_60_15_30alt_18spp.RData") #ann_30
  

ann_30 <- ann_30 %>%  
  mutate(group = paste(common_name, colony.name, sep = "_")) %>% 
  as.data.frame()

### STEP 2: Calculate within-stratum variances ------------------------------------ #####

#calculate and plot within stratum variances
data_var <- ann_30 %>%
  group_by(stratum) %>%
  summarise(wspd_var = var(wind_speed),
            u_var = var(u10m),
            v_var = var(v10m),
            wspt_var = var(wind_support),
            species = head(common_name, 1),
            group = head(group,1),
            year = head(year, 1),
            wspd_cov = (sd(wind_speed)/mean(wind_speed))*100) %>% #coef of variation
  mutate(log_wspd_cov = log(wspd_cov)) %>% 
  ungroup() %>% 
  as.data.frame()

save(data_var, file = "R_files/data_var.RData")


#plot
ggplot(data_var, aes(x = wspd_cov, y = group)) + 
  geom_density_ridges(scale = 3, alpha = 0.5) + 
  scale_x_continuous(limits = c(-1, 30)) +
  labs(y = "", x = "Coefficient of variation (%)") +
  theme_minimal()

data_var[data_var$species == "Galapagos albatross (waved albatross)","species"] <- "Galapagos albatross"

X11(width = 5, height = 5)

CoV_bar <- ggplot(data_var, aes(x = wspd_cov, y = reorder(as.factor(species), desc(as.factor(species))), height = stat(density))) + 
  geom_density_ridges(
    stat = "binline", bins = 20, scale = 0.98, alpha = 0.3,
    draw_baseline = FALSE
  ) +
  scale_x_continuous(limits = c(-4, 140)) +
  labs(y = "", x = "Coefficient of variation (%)") +
  theme_minimal()


png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/Cov_bar.png", width = 5, height = 5, units = "in", res = 300)
print(CoV_bar)
dev.off()


### STEP 3: Permutation test ------------------------------------ #####
#In each stratum, is the difference between selected and max available wind speed higher/lower than expected by chance?

#for each stratum, calculate the difference between observed wind speed and max wind speed (incl. observed)

observed_stat <- ann_30 %>% 
  group_by(group, year, stratum) %>% #group by species-colony instead of just species
  arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached  detach("package:plyr", unload=TRUE)
  summarize(max_minus_obs = max(wind_speed) - head(wind_speed,1))

save(observed_stat, file = "R_files/observed_stats.RData")

#randomize and calculate the same statistic
#shuffle all wind speed values within each year, then recalc the statistic within each stratum

permutations <- 1000

#prep cluster
mycl <- makeCluster(detectCores() - 2, setup_strategy = "sequential")
clusterExport(mycl, c("permutations", "ann_30")) 

clusterEvalQ(mycl, {
  library(dplyr)
})


a <- Sys.time()

rnd_stat <- parLapply(cl = mycl, X = c(1:permutations), fun = function(x){ 
  
  #rnd_stat <- lapply(1:permutations, function(x){
  ann_30 %>% 
    group_by(group,year) %>% 
    mutate(wind_speed = sample(wind_speed, replace = F)) %>% 
    group_by(group,year,stratum) %>% 
    arrange(desc(used), .by_group = TRUE) %>%
    summarize(max_minus_obs = max(wind_speed) - head(wind_speed,1)) %>% 
    mutate(perm = x)
  
}) %>% 
  reduce(rbind) %>% 
  as.data.frame()

Sys.time() - a # 2.61326 hours (1000 perms)

stopCluster(mycl)

#save(rnd_stat, file = "R_files/rnd_stats_1000_perm_df.RData")

load("R_files/observed_stats.RData")
load("R_files/rnd_stats_1000_perm_df.RData")

#extract observed and random values for each stratum

#prep cluster
mycl <- makeCluster(detectCores() - 7, setup_strategy = "sequential")
clusterExport(mycl, c("permutations", "observed_stat", "rnd_stat")) 

clusterEvalQ(mycl, {
  library(dplyr)
})


a <- Sys.time()

p_vals <- parLapply(mycl, unique(observed_stat$stratum), function(x){
  obs <- observed_stat[observed_stat$stratum == x,]
  rnd <- rnd_stat[rnd_stat$stratum == x,]
  
  obs$p_less <- sum(rnd$max_minus_obs <= obs$max_minus_obs)/permutations
  obs$p_more <- sum(rnd$max_minus_obs >= obs$max_minus_obs)/permutations
  
  obs
}) %>% 
  reduce(rbind)
  

Sys.time() - a # 9.972578 hours

stopCluster(mycl)


save(p_vals, file = "R_files/p_vals_1000_perm_df.RData")


### STEP 4: Plot permutation results ------------------------------------ #####
 
load("R_files/p_vals_1000_perm_df.RData") #p_vals

p_vals <- p_vals %>% 
  mutate(year_f = as.character(year),
         species = str_split(group, "\\_", simplify = TRUE)[, 1]) %>% 
  ungroup() %>% 
  as.data.frame()

p_vals[p_vals$species == "Galapagos albatross (waved albatross)", "species"] <- "Galapagos albatross"


X11(width = 5, height = 5)

perm_sig <- ggplot(p_vals, aes(x = p_more, y = reorder(as.factor(species), desc(as.factor(species))), height = stat(density))) + 
  #geom_density_ridges(scale = 2, alpha = 0.5,
  #                    jittered_points = TRUE,
  #                    position = position_points_jitter(width = 0.05, height = 0),
  #                    point_shape = '|', point_size = 2, point_alpha = 1) + 
  geom_density_ridges(
   stat = "binline", bins = 20, scale = 0.98, alpha = 0.3,
    draw_baseline = FALSE
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "", x = "Significance") +
  theme_minimal()

png("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/Perm_sig.png", width = 5, height = 5, units = "in", res = 300)
print(perm_sig)
dev.off()


#extract strata with significant avoidance of strong winds

# #extract strata with p-value less than 0.05
sig <- p_vals %>%
  filter(p_more <= 0.05)

sig_data <- ann_30 %>%
  filter(stratum %in% sig$stratum) %>% 
  mutate(species = fct_relevel(as.factor(common_name), levels = "Atlantic Yellow-nosed Albatross", "Wandering albatross", 
                                 "Sooty Albatross", "Red-footed booby")) %>% 
  mutate(col = as.character(fct_relevel(species), levels = "corn flower blue", "rosy brown", "yellow green", "pale violet red")) %>% 
  as.data.frame()

sig_data$common_name <- reorder(sig_data$common_name, sig_data$species)
sig_data$stratum <- reorder(as.factor(sig_data$stratum),desc(sig_data$species))

save(sig_data, file = "R_files/sig_data.RData")

#plot
used_wind <- sig_data %>% 
  group_by(stratum) %>% 
  filter(used ==1) %>% 
  ungroup() %>% 
  dplyr::select(c("stratum", "common_name", "wind_speed")) %>% 
  as.data.frame()


X11(width = 7, height = 4)
sig_plots <- ggplot(sig_data, aes(x = wind_speed, y = stratum)) + 
  geom_density_ridges(scale = 3, alpha = 0.4,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 2, point_alpha = 1, aes(fill = species)) + 
  scale_fill_manual(values = c("Atlantic Yellow-nosed Albatross" = "corn flower blue", "Wandering albatross" = "yellowgreen", 
                               "Sooty Albatross" = "lightcoral", "Red-footed booby" = "goldenrod")) +
  scale_x_continuous(limits = c(0, 25)) +
  geom_segment(data = used_wind, aes(x = wind_speed, xend = wind_speed, y = as.numeric(as.factor(stratum)),
                                     yend = as.numeric(as.factor(stratum)) + .9, linetype = "Selected wind speed"), color = "red") +
  
  #scale_linetype_manual("",values = c("Selected wind speed" = 1)) +
  guides(fill = guide_legend(order = 1),
         line = guide_legend(order = 2)) +
  labs(y = "Density", x = "Wind speed (m/s)") +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank())


 

### STEP 5: Plot wind fields for wind avoidance steps ------------------- #####

load("R_files/sig_data.RData") #sig_data

#focus on Atlantic yellow-nosed albatross first. For the other species, see wind_fields.R
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/Peter_Ryan_data_annotated_SplitTrip.Rdata") #PR_data_split

ayl <- PR_data_split %>% 
  filter(common_name == "Atlantic Yellow-nosed Albatross", TripID == sig_data[sig_data$common_name == "Atlantic Yellow-nosed Albatross", "TripID"]) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)

sig_data$stratum <-as.character(sig_data$stratum)

avoidance_hour <- sig_data %>% 
  filter(species == "Atlantic Yellow-nosed Albatross") %>% 
  rowwise() %>% 
  mutate(unique_hour = paste(yday(timestamp), hour(timestamp), sep = "_"))


#extract the whole trips (not only steps that were retained in the step selection dataset) that contain rare event steps and extract extents
# rare_trips <- ann_30 %>%
#   filter(used == 1 & TripID %in% sig_data$TripID) %>% 
#   st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)
 
#--------------------request data download -----

#only for month and dates corresponding with the trip

#import the python library ecmwfapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()

output_path <- "/home/enourani/Documents/ERA_5_seabirds/ayl/"

datetimes <- list("11" = str_pad(17:30,2,"left","0"),
                 "12" = str_pad(1:3,2,"left","0"))
  
lapply(c(1:length(datetimes)), function(x){
  
  month <- names(datetimes)[[x]]
  days <- datetimes[[x]]
    
    request <- r_to_py(list(
      product_type = "reanalysis",
      variable = c("10m_u_component_of_wind", "10m_v_component_of_wind", "significant_height_of_combined_wind_waves_and_swell", "surface_pressure"),
      year = "2014",
      month = month,
      day = days,
      time = str_c(seq(0,23,1),"00",sep=":") %>% str_pad(5,"left","0"),
      area = c(-25, -62, -51, 11),
      format = "netcdf",
      dataset_short_name = "reanalysis-era5-single-levels"
    ))
  
    server$retrieve("reanalysis-era5-single-levels",
                    request,
                    target = paste0(output_path,"wind_2014_", month, ".nc")) 
})

#--------------------extract date from netcdf files -----
setwd("/home/enourani/Documents/ERA_5_seabirds/ayl/")

file_list <- list.files(pattern = ".nc",full.names = TRUE)
vname <- c("u10","v10", "swh", "sp")

data_list <- lapply(file_list,function(x){
  
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
    matrix(as.vector(ncvar_get(nc,vname[3])), nrow = nlon * nlat * nt, ncol = 1),
    matrix(as.vector(ncvar_get(nc,vname[4])), nrow = nlon * nlat * nt, ncol = 1)))
  
  colnames(var_df) <- c("lon","lat","date_time",vname)   #set column names
  
  #remove points over land (NAs)
  
  df <- var_df %>%
    #na.omit() %>% #let's keep points over land as well. otherwise, points where significant wave height is NA, i.e. land, will be deleted
    mutate(yday = yday(date_time),
           hour = hour(date_time),
           year = year(date_time)) %>%
    data.frame()
  
  save(df,file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/atlantic_yellow_nosed_albatross_",head(month(df$date_time),1), ".RData"))
})

#--------------------append the two files -----

file_ls <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/",".RData",full.names = TRUE)

data_df <- sapply(file_ls, function(x) mget(load(x)), simplify = TRUE) %>%
  reduce(rbind)

save(data_df, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/seabirds_storms/atlanitc_yellow_nosed_raw_wind.RData")

#--------------------PLOT!!! -----
#https://semba-blog.netlify.app/10/29/2018/animating-oceanographic-data-in-r-with-ggplot2-and-gganimate/

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/atlanitc_yellow_nosed_raw_wind.RData") #data_df

data_df <- data_df %>% 
  mutate(wind_speed = sqrt(u10^2 + v10^2)) #m/s  
  

ayl <- PR_data_split %>% 
  filter(common_name == "Atlantic Yellow-nosed Albatross", 
         TripID == sig_data[sig_data$common_name == "Atlantic Yellow-nosed Albatross", "TripID"]) %>% 
  rowwise() %>% 
  mutate(unique_hour = paste(yday(timestamp), hour(timestamp), sep = "_")) %>% 
  mutate(avoidance = ifelse(unique_hour %in% avoidance_hour$unique_hour, "avoided", "not_avoided"))

region <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -62, xmax = 11, ymin = -51, ymax = -25) %>%
  st_union()

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
  
save(wind_lres, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/atlanitc_yellow_nosed_raw_wind_lres.RData")

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


#animate
wind_fields =  ggplot() +
  geom_raster(data = wind_lres, aes(x = lon, y = lat, fill = wind_speed))+
  geom_segment(data = wind_lres, aes(x = lon, xend = lon+u10/10, y = lat, 
                   yend = lat+v10/10), arrow = arrow(length = unit(0.12, "cm")), size = 0.3)+
  geom_sf(data = region, fill = "grey85", col = 1)+
  #eom_point(data = ayl %>%  filter(unique_hour == i), aes(x = location.long, y = location.lat), 
  #           size = 1, colour = "red") +
  coord_sf(xlim = c(-62, 11), ylim =  c(-51, -25))+
  scale_fill_gradientn(colours = oce::oceColorsPalette(120), limits = c(0,23), 
                       na.value = "white", name = "Speed\n (m/s)")+
  theme_bw()+
  theme(axis.text = element_text(size = 14, colour = 1),
        legend.text = element_text(size = 14, colour = 1), 
        legend.title = element_text(size = 14, colour = 1),
        legend.position = c(.12,.17),
        legend.background = element_rect(colour = 1, fill = "white"))+
  labs(x = NULL, y = NULL, title = "Date : {frame_time}") +
  transition_time(date_time) +
  ease_aes("linear")



wind.vector = ggplot() +
  geom_raster(data = data_df, aes(x = lon, y = lat, fill = wind_speed))+
  geom_segment(data = data_df, 
               aes(x = lon, xend = lon+u10/60, y = lat, 
                   yend = lat+v10/60), arrow = arrow(length = unit(0.1, "cm")))+
  geom_sf(data = region, fill = "grey85", col = 1)+
  coord_sf(xlim = c(-62, 11), ylim =  c(-51, -25))+
  scale_fill_gradientn(colours = oce::oceColorsPalette(120), limits = c(0,12), 
                       na.value = "white", name = "Speed\n (m/s)")+
  scale_x_continuous(breaks = c(38.8,40))+
  theme_bw()+
  theme(axis.text = element_text(size = 14, colour = 1),
        legend.text = element_text(size = 14, colour = 1), 
        legend.title = element_text(size = 14, colour = 1),
        legend.position = c(.12,.17),
        legend.background = element_rect(colour = 1, fill = "white"))+
  labs(x = NULL, y = NULL, title = "Month of : {frame_time}")+
  transition_time(hour) +
  ease_aes("linear")

animate(wind.vector)






### STEP 6: Raw wind plots ------------------------------------ #####

#open annotated data and species data
load("R_files/lm_input_20spp_col.RData") #lm_input 
load("R_files/ann_18spp.RData") #ann (from PCoA_seabirds.R)
#species <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/final_datasets.csv") %>% 
#  filter(scientific.name %in% c( "Ardenna gravis", "Diomedea dabbenena", "Diomedea exulans", "Fregata magnificens", "Fregata minor", "Morus bassanus", "Morus capensis", "Phaethon lepturus", "Phaethon rubricauda", 
#                                 "Phoebastria irrorata", "Phoebetria fusca", "Procellaria cinerea", "Pterodroma incerta", "Pterodroma mollis",  "Sula dactylatra", "Sula granti",  "Sula sula", "Thalassarche chlororhynchos"))

#add wing loading for great shearwater (decided to use the sooty shearwater value)
lm_input[lm_input$species == "Great Shearwater", "wing.loading..Nm.2."] <- 88
save(lm_input, file = "R_files/lm_input_20spp_col.RData")


ann <- ann %>%  
  left_join(lm_input[,c(1,2,9:13)], by = "sci_name") %>% 
  mutate(group = paste(species, colony.name, sep = "_")) %>% 
  mutate(group_f = as.factor(reorder(group, desc(wing.loading..Nm.2.)))) %>% 
  as.data.frame()

summ <- ann %>% 
  group_by(group_f) %>% 
  summarize(min_wind = min(wind_speed_ms, na.rm = T),
            max_wind = max(wind_speed_ms, na.rm = T)) %>% 
  as.data.frame()
  
  
#plot
#https://www.datanovia.com/en/blog/elegant-visualization-of-density-distribution-in-r-using-ridgeline/

cols <- oce::oceColorsPalette(10)

#with quantiles
X11(width = 8, height = 7)
  raw_wind <- ggplot(ann, aes(x = wind_speed_ms, y = group_f, fill = factor(stat(quantile)))) + 
    stat_density_ridges(jittered_points = TRUE, rel_min_height = .01,
                        point_shape = "|", point_size = 1, point_alpha = 0.7, size = 0.25,
      geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F,
      quantiles = 10, quantile_lines = F, scale = 3) +
    scale_fill_viridis_d(name = "Quantiles", alpha = 0.6) +
    scale_x_continuous(limits = c(-0.5, 25)) +
    #scale_fill_manual(values = oce::oceColorsPalette(10)) + #this works, but colors are representing quantiles, not wind speeds
    labs(y = "", x = "Wind speed (m/s)") +
    theme_minimal()
  
png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/raw_wind_jitter_desc_all18.png", width = 8, height = 8, units = "in", res = 300)
print(raw_wind)
dev.off()


### STEP 7: LM: wind strength ------------------------------------ #####

load("R_files/lm_input_20spp_col.RData") #lm_input (from PCoA_seabirds) #was sent to Emily too

lm_input <- lm_input %>% 
  mutate_at(c("colony.long","colony.lat", "wing.loading..Nm.2.", "wing.span..m.", "wing.area..m2.","aspect.ratio", "PC1", "PC2"),
            list(z = ~as.numeric(scale(.)))) %>% 
  mutate(breeding_ch = as.character(median_breeding_m),
         max_wind_ms = max_wind/3.6) %>% 
  arrange(median_breeding_yday) #order temporally, for temporal autocorrelation analysis.

lm_input %>% 
  dplyr::select(c("colony.long","colony.lat","median_breeding_m","median_breeding_yday", "wing.loading..Nm.2.",
                  "aspect.ratio","wing.area..m2.", "wing.span..m.", "PC1", "PC2")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #PC1 and PC2 are correlated; median month and yday are correlated. aspect ratio correlated with wing loading and wing span

## add data quantity
load("R_files/ann_18spp.RData") #ann; from PCoA_seabirds.R

#extract number of data points and number of trips
summary_info <- ann %>% 
  group_by(sci_name, colony.name) %>% 
  summarize(n_rows = n(),
            n_trips = n_distinct(TripID)) %>% 
  full_join(lm_input, by = c("sci_name", "colony.name"))

#modeling

#plot the relationship
ggplot(lm_input,aes(wing.loading..Nm.2., max_wind_ms)) +
  geom_point() +
  #stat_summary(fun.data = mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x, se = T) +
  theme_minimal()

ggplot(lm_input,aes(wing.area..m2., max_wind_ms)) +
  geom_point() +
  #stat_summary(fun.data = mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x, se = T) +
  theme_minimal()

#model
morph <- lm(max_wind_ms ~ wing.loading..Nm.2., data = lm_input) #adR = 0.2809 , AIC =  120.0459

#investigate residuals
par(mfrow = c(2,2))
plot(morph) #residuals are fine. 3rd plot: variance is higher at values between 50-60. perhaps because of large sample size there.


#plot residuals against other variables
plot(resid(morph) ~ wing.loading..Nm.2., data = lm_input)
plot(resid(morph) ~ wing.area..m2., data = lm_input)
plot(resid(morph) ~ colony.lat, data = lm_input)
plot(resid(morph) ~ colony.long, data = lm_input)
plot(resid(morph) ~ as.factor(breeding_ch), data = lm_input) #a little problematic. medians and variances are not equal
plot(resid(morph) ~ n_trips, data = summary_info)

#temporal correlation. result: no temporal autocorrelation!
acf(resid(morph))
acf(resid(morph),type = "p")

#spatial autocorrelation. bubble plot result: no autocorrelation
spdata <- data.frame(resid = resid(morph), x = lm_input$colony.long, y = lm_input$colony.lat)
coordinates(spdata)<-~ x + y
bubble(spdata, "resid", col = c("blue","orange"))

#check for phylogeny

#100 trees downloaded from https://birdtree.org/subsets/ for the 18 species
trees <- read.nexus("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/phylogeny_tree/tree-pruner-5edeac95-d8ed-410a-8162-982b1438f4e9/output.nex")
#use the first tree
plot(trees[[1]])
axisPhylo()
w <- 1/cophenetic(trees[[1]])
diag(w) <- 0 #set the diagonal to 0

#extract wind speed
#change scie name for great shearwater to match the sci name in the tress
lm_input[lm_input$sci_name == "Ardenna gravis", "sci_name"] <- "Puffinus gravis"
max_wspd <- lm_input$max_wind_ms
names(max_wspd) <- lm_input$sci_name

#estimate Moran's I
Moran.I(max_wspd, w)

#To Do:
#pool the two colonies for RFB and MF. add great shearwater

### STEP 8: LM: wind variability ------------------------------------ #####

load("R_files/lm_input_20spp_col.RData") #lm_input (from PCoA_seabirds) #was sent to Emily too
load("R_files/data_var.RData") #data_var

summary_info <- summary_info %>% 
  mutate(group = paste(species, colony.name, sep = "_"))

str_var <- data_var %>% 
  group_by(group) %>% 
  summarize(max_str_cov = max(wspd_cov)) %>% 
  full_join(summary_info, by = "group") %>% 
  as.data.frame()

save(str_var, file = "R_files/str_var.RData")
  
  

#lm_var <- lm(max_str_cov ~ wing.loading..Nm.2. + wing.area..m2., data = str_var)

m1 <- lm(max_str_cov ~ wing.loading..Nm.2. + scale(n_trips), data = str_var) # 0.2758  #AIC =  181.7002

# m1 <- lm(max_str_cov ~ wing.loading..Nm.2., data = str_var) #0.3199 

# m2 <- lm(max_str_cov ~  wing.area..m2., data = str_var) #0.1668 
# 
# m3 <- lm(max_str_cov ~  aspect.ratio, data = str_var) #0.2779 
# 
# m4 <- lm(max_str_cov ~  PC1, data = str_var) #0.375  
# 
# m5 <- lm(max_str_cov ~ body.mass..kg., data = str_var) # 0.2462 

ggplot(str_var,aes(wing.area..m2., max_str_cov)) +
  geom_point() +
  #stat_summary(fun.data = mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x, se = T) +
  theme_minimal()

### STEP 9: Plot linear models in one device ------------------------------------ #####

load("R_files/str_var.RData")

str_var <- str_var %>% 
  mutate(max_wind_ms = max_wind/3.6)

# m1 <- lm(str_var$max_str_cov ~ str_var$wing.loading..Nm.2.) #0.3199 
# #https://stackoverflow.com/questions/46459620/plotting-a-95-confidence-interval-for-a-lm-object
# 
# #create new data for plotting
# newx <- seq(min(str_var$wing.loading..Nm.2., na.rm = T), max(str_var$wing.loading..Nm.2., na.rm = T), by = 0.05)
# conf_interval <- predict(m1, newdata = data.frame(wing.loading..Nm.2.= newx), interval = "confidence",
#                          level = 0.95)
# 
# plot(x = c(30,140), y = c(2,140), type = "n", xlab = "", ylab = "")
# 
# abline(m1,col = max_col)
# matlines(conf_interval[-6,1], conf_interval[-6,2:3], col = "blue", lty = 2) #row 6 has NAs

max_col <- "corn flower blue" #"#69b3a2"
var_col <- "goldenrod"  # "#c7909d"

#

m1 <- lm(max_wind_ms ~ wing.loading..Nm.2., data = str_var) #0.3199 

ggplotRegression(m1)

ggplot(str_var, aes(wing.loading..Nm.2., max_str_cov))+
  geom_point() +
  # Add the line using the fortified fit data, plotting the x vs. the fitted values
  geom_line(data = fortify(m1), aes(x = wing.loading..Nm.2., y = .fitted))


#####
X11(width = 12, height = 5)
lm_maxwind <- ggplot(str_var, aes(x = wing.loading..Nm.2.)) +
  geom_smooth(aes(y = max_wind_ms), method = "lm", color = max_col, alpha = .2, fill = max_col) +
  geom_point(aes(y = max_wind_ms), color = max_col, alpha = .6) +
  labs(x = "Wind loading") +
  scale_y_continuous(
    name = "Maximum wind speed (m/s)") +
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(size=13))

lm_covwind <- ggplot(str_var, aes(x = wing.loading..Nm.2.)) +
  geom_smooth(aes(y = max_str_cov), method = "lm", color = var_col, alpha = .2, fill = var_col) +
  geom_point(aes(y = max_str_cov), color = var_col, alpha = .6) +
  labs(x = "Wind loading") +
  scale_y_continuous(
    name = "Variation in wind speed (%)") +# Features of the first axis
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(size=13))

 grid.arrange(lm_maxwind, lm_covwind, nrow = 1)



png("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/lm_output_two_panels.png", 
    width = 12, height = 5, units = "in", res = 300)
grid.arrange(lm_maxwind, lm_covwind, nrow = 1)
dev.off()


### STEP 10: Plot relationship between wind speed and latitude ------------------------------------ #####


#open dataset with alternative steps for hourly steps. prepped in random_steps.R
load("R_files/ssf_input_annotated_60_15_30alt_18spp.RData") #ann_30

observed <- ann_30 %>%  
  filter(used == 1) %>% 
  mutate(group = paste(common_name, colony.name, sep = "_")) %>% 
  as.data.frame()

#extract latitudinal ranges for each group
lat_ranges <- observed %>% 
  group_by(group) %>% 
  summarize(min_lat = min(location.lat),
            max_lat = max(location.lat),
            colony.name = head(colony.name,1),
            species = head(common_name,1)) #timestamp? 

Pal <- colorRampPalette(c( "darkslategray1", "darkgoldenrod1","lightpink1")) 
Cols <- paste0(Pal(20), "80")
#Cols <- palette(rainbow(20))

#plot
plot(observed$location.lat, observed$wind_speed, pch = 20, col = "gray")

for (i in 1:nrow(lat_ranges)){
  rect(xleft = lat_ranges[i,"min_lat"], ybottom = min(observed$wind_speed), xright = lat_ranges[i,"max_lat"], ytop = max(observed$wind_speed), 
       col = Cols[i], border = NA)
}


plot(observed$wind_speed, observed$location.lat, pch = 20, col = "gray")

for (i in 1:nrow(lat_ranges)){
  rect(xleft =  min(observed$wind_speed), ybottom =lat_ranges[i,"min_lat"], xright = max(observed$wind_speed), ytop =lat_ranges[i,"max_lat"] , 
       col = Cols[i], border = NA)
}

plot(observed$wind_speed, observed$location.lat, pch = 20, cex = 0.1, col = "orange")

for (i in 1:nrow(lat_ranges)){
  rect(xleft =  min(observed$wind_speed), ybottom =lat_ranges[i,"min_lat"], xright = max(observed$wind_speed), ytop =lat_ranges[i,"max_lat"] , 
       col = alpha("grey",0.1), border = NA)
}

points(observed$wind_speed, observed$location.lat, pch = 20, cex = 0.1, col = "orange")





### STEP 11: Correlation between data quantity and max wind speeds ------------------------------------ #####

#load annotated data (one hourly subset)
load("R_files/ann_18spp.RData") #ann; from PCoA_seabirds.R
load("R_files/lm_input_20spp_col.RData") #lm_input 

#extract number of data points and number of trips
summary_info <- ann %>% 
  group_by(sci_name, colony.name) %>% 
  summarize(n_rows = n(),
            n_trips = n_distinct(TripID)) %>% 
  full_join(lm_input, by = c("sci_name", "colony.name")) %>% 
  as.data.frame()

#are nrows and n trips correlated?
cor(summary_info$n_trips,summary_info$n_rows) #yes! 0.88

#are max wind speeds correlated with nrows?
#USED IN MS
cor.test(summary_info$n_trips,summary_info$max_wind_ms) #0.15
cor.test(summary_info[-3, "n_trips"], summary_info[-3,"max_wind_ms"]) #-0.3

cor.test(str_var$n_trips, str_var$max_str_cov) #0.32
cor.test(str_var[-18,"n_trips"], str_var[-18,"max_str_cov"]) # -0.15


#cor(summary_info$n_rows,summary_info$max_wind) #0.3


#plot without wandering albatross
plot(max_wind ~ n_trips, data = summary_info[-3,])
plot(max_wind ~ n_trips, data = summary_info)

plot(max_str_cov ~ n_trips, data = str_var[-18,])
plot(max_str_cov ~ n_trips, data = str_var)

## use n_rows instead
cor.test(summary_info$n_rows,summary_info$max_wind_ms) #0.3
cor.test(summary_info[-3, "n_rows"], summary_info[-3,"max_wind_ms"]) #-0.4

cor.test(str_var$n_rows, str_var$max_str_cov) #0.36
cor.test(str_var[-18,"n_rows"], str_var[-18,"max_str_cov"]) # -0.4

plot(max_wind ~ n_rows, data = summary_info[-3,])
plot(max_wind ~ n_rows, data = summary_info)

plot(max_str_cov ~ n_rows, data = str_var[-18,])
plot(max_str_cov ~ n_rows, data = str_var)