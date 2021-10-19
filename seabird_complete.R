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
#install.packages('devtools')
#devtools::install_github('mpio-be/windR')
library(windR)

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

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

#extract date and time of the rare events
rare_times <- ann_30 %>%
  filter(used == 1 & TripID %in% sig_data$TripID) %>% 
  group_by(TripID) %>% 
  summarize(yr = head(year,1),
            min_t = min(timestamp),
            max_t = max(timestamp)) %>% 
  ungroup() %>% 
  as.data.frame()

#extract the whole trips that contain rare event steps and extract extents
rare_trips <- ann_30 %>%
  filter(used == 1 & TripID %in% sig_data$TripID) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)
 
extents <- lapply(split(rare_trips, rare_trips$TripID), function(x) st_bbox(x))

#download era-5 data for the time periods (of entire tracks)
#import the python library ecmwfapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()


output_path <- "/home/enourani/Documents/ERA_5_seabirds/"

lapply(c(2:length(extents)), function(zone){
  
  x <- extents[[zone]]+1
  yr <- as.character(rare_times[zone,"yr"])
  if(month(rare_times[zone, "min_t"]) == month(rare_times[zone, "max_t"])){
    mn <- as.character(month(rare_times[zone, "min_t"]))
  } else {
    mn <- c(as.character(month(rare_times[zone, "min_t"])), as.character(month(rare_times[zone, "max_t"])))
  }
  
  for(mn in mn){ #for each month
  
  query <- r_to_py(list(
    product_type = "reanalysis",
    area =  c(x[4], x[1], x[2], x[3]),     # North, West, South, East.
    #grid = c(0.25, 0.25), 
    format = "netcdf",
    variable = c("10m_u_component_of_wind", "10m_v_component_of_wind", "surface_pressure"),
    year = yr,
    month = mn, #"01", #c(str_pad(seq(1:9),2,"left","0"),"10","11","12"),
    day = str_pad(1:31,2,"left","0"),
    time = str_c(seq(0,23,1),"00",sep=":") %>% str_pad(5,"left","0"), #every hour
    dataset = "reanalysis-era5-single-levels"
  ))
  
  server$retrieve("reanalysis-era5-single-levels",
                  query,
                  target = paste0(output_path,"wind_", names(extents)[[zone]], "_", yr, "_", mn, ".nc")) 
  }
  
})


#extract date from netcdf files



### STEP 6: Raw wind plots ------------------------------------ #####

#open annotated data and species data
load("R_files/ann_18spp.RData") #ann (from PCoA_seabirds.R)
species <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/final_datasets.csv") %>% 
  filter(scientific.name %in% c( "Ardenna gravis", "Diomedea dabbenena", "Diomedea exulans", "Fregata magnificens", "Fregata minor", "Morus bassanus", "Morus capensis", "Phaethon lepturus", "Phaethon rubricauda", 
                                 "Phoebastria irrorata", "Phoebetria fusca", "Procellaria cinerea", "Pterodroma incerta", "Pterodroma mollis",  "Sula dactylatra", "Sula granti",  "Sula sula", "Thalassarche chlororhynchos"))

ann <- ann %>%  
  left_join(species[,c(1,2)], by = c("sci_name" = "scientific.name")) %>% 
  mutate(group = paste(species, colony.name, sep = "_")) %>% 
  as.data.frame()

#plot

# ggplot(ann, aes(x = wind_speed_ms, y = group)) + 
#   geom_density_ridges(scale = 3, alpha = 0.4) +
#   scale_x_continuous(limits = c(0, 25)) +
# 
#   labs(y = "Density", x = "Wind speed (m/s)") +
#   theme_minimal() +
#   theme(legend.position = "none") 

#with mean
# ggplot(ann, aes(x = wind_speed_ms, y = group)) + 
#   stat_density_ridges(
#     geom = "density_ridges_gradient", calc_ecdf = TRUE, 
#     scale = 3, alpha = 0.4,
#     quantiles = 0.5, quantile_lines = TRUE) +
#   scale_x_continuous(limits = c(0, 25)) +
#   labs(y = "", x = "Wind speed (m/s)") +
#   theme_minimal() +
#   theme(legend.position = "none") 

#with max
X11(width = 8, height = 7)
raw_wind <- ggplot(ann, aes(x = wind_speed_ms, y = group)) + 
  stat_density_ridges(quantile_lines = TRUE, quantiles = 1, alpha = 0.7) +
  #scale_fill_grey(name = "Quantiles", alpha = 0.6) +
  scale_fill_viridis_d(name = "Quantiles", alpha = 0.6) +
  scale_x_continuous(limits = c(0, 50)) +
  labs(y = "", x = "Wind speed (m/s)") +
  theme_minimal() #+
theme(legend.position = "bottom")

#with quartiles
X11(width = 8, height = 7)
raw_wind <- ggplot(ann, aes(x = wind_speed_ms, y = group, fill = factor(stat(quantile)))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 10, quantile_lines = F, scale = 3,
  ) +
  #scale_fill_grey(name = "Quantiles", alpha = 0.6) +
  scale_fill_viridis_d(name = "Quantiles", alpha = 0.6) +
  scale_x_continuous(limits = c(0, 25)) +
  labs(y = "", x = "Wind speed (m/s)") +
  theme_minimal() #+
theme(legend.position = "bottom") 

png("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/raw_wind.png", width = 8, height = 7, units = "in", res = 300)
print(raw_wind)
dev.off()

### STEP 7: LM: wind strength ------------------------------------ #####

load("R_files/lm_input_20spp_col.RData") #lm_input (from PCoA_seabirds) #was sent to Emily too

lm_input <- lm_input %>% 
  mutate_at(c("colony.long","colony.lat", "wing.loading..Nm.2.", "wing.span..m.", "wing.area..m2.","aspect.ratio", "PC1", "PC2"),
            list(z = ~as.numeric(scale(.)))) %>% 
  mutate(breeding_ch = as.character(median_breeding_m)) %>% 
  arrange(median_breeding_yday) #order temporally, for temporal authocorrelation analysis.

lm_input %>% 
  dplyr::select(c("colony.long","colony.lat","median_breeding_m","median_breeding_yday", "wing.loading..Nm.2.",
                  "aspect.ratio","wing.area..m2.", "wing.span..m.", "PC1", "PC2")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #PC1 and PC2 are correlated; median month and yday are correlated. aspect ratio correlated with wing loading and wing span

#modeling

#plot the relationship
ggplot(lm_input,aes(wing.loading..Nm.2., max_wind)) +
  geom_point() +
  #stat_summary(fun.data = mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x, se = T) +
  theme_minimal()

ggplot(lm_input,aes(wing.area..m2., max_wind)) +
  geom_point() +
  #stat_summary(fun.data = mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x, se = T) +
  theme_minimal()

#model
morph <- lm(max_wind ~ wing.loading..Nm.2., data = lm_input) #adR = 0.2898, AIC =  161.2964
#morph2 <- lm(max_wind ~ wing.loading..Nm.2._z + wing.area..m2._z, data = lm_input) #0.3096 , AIC = 161.6079

#investigate residuals
par(mfrow = c(2,2))
plot(morph) #residuals are fine. 3rd plot: variance is higher at values between 50-60. perhaps because of large sample size there.

#X11()
#par(mfrow = c(2,2))
#plot(morph2) #residual plots are worse than morph

#plot residuals agains other variables
plot(resid(morph) ~ wing.loading..Nm.2.[-1], data = lm_input)
plot(resid(morph) ~ wing.area..m2.[-1], data = lm_input)
plot(resid(morph) ~ colony.lat[-1], data = lm_input)
plot(resid(morph) ~ colony.long[-1], data = lm_input)
plot(resid(morph) ~ as.factor(breeding_ch[-1]), data = lm_input) #a little problematic. medians and variances are not equal

#temporal correlation. result: no temporal autocorrelation!
acf(resid(morph))
acf(resid(morph),type = "p")

#spatial autocorrelation. bubble plot result: no autocorrelation
spdata <- data.frame(resid = resid(morph), x = lm_input$colony.long[-1], y = lm_input$colony.lat[-1])
coordinates(spdata)<-~ x + y
bubble(spdata, "resid", col = c("blue","orange"))
#try a semivariogram?

#To Do:
#pool the two colonies for RFB and MF. add great shearwater

### STEP 8: LM: wind variability ------------------------------------ #####

load("R_files/lm_input_20spp_col.RData") #lm_input (from PCoA_seabirds) #was sent to Emily too
load("R_files/data_var.RData") #data_var

lm_input <- lm_input %>% 
  mutate(group = paste(species, colony.name, sep = "_"))

str_var <- data_var %>% 
  group_by(group) %>% 
  summarize(max_str_cov = max(wspd_cov)) %>% 
  full_join(lm_input, by = "group") %>% 
  as.data.frame()

save(str_var, file = "R_files/str_var.RData")
  
  


lm_var <- lm(max_str_cov ~ wing.loading..Nm.2. + wing.area..m2., data = str_var)

m1 <- lm(max_str_cov ~ wing.loading..Nm.2., data = str_var) #0.3199 

m2 <- lm(max_str_cov ~  wing.area..m2., data = str_var) #0.1668 

m3 <- lm(max_str_cov ~  aspect.ratio, data = str_var) #0.2779 

m4 <- lm(max_str_cov ~  PC1, data = str_var) #0.375  

m5 <- lm(max_str_cov ~ body.mass..kg., data = str_var) # 0.2462 

ggplot(str_var,aes(wing.area..m2., max_str_cov)) +
  geom_point() +
  #stat_summary(fun.data = mean_cl_normal) + 
  geom_smooth(method='lm', formula= y~x, se = T) +
  theme_minimal()

### STEP 9: Plot linear models in one device ------------------------------------ #####

load("R_files/str_var.RData")

m1 <- lm(str_var$max_str_cov ~ str_var$wing.loading..Nm.2.) #0.3199 
#https://stackoverflow.com/questions/46459620/plotting-a-95-confidence-interval-for-a-lm-object

#create new data for plotting
newx <- seq(min(str_var$wing.loading..Nm.2., na.rm = T), max(str_var$wing.loading..Nm.2., na.rm = T), by = 0.05)
conf_interval <- predict(m1, newdata = data.frame(wing.loading..Nm.2.= newx), interval = "confidence",
                         level = 0.95)

plot(x = c(30,140), y = c(2,140), type = "n", xlab = "", ylab = "")

abline(m1,col = max_col)
matlines(conf_interval[-6,1], conf_interval[-6,2:3], col = "blue", lty = 2) #row 6 has NAs

max_col <- "corn flower blue" #"#69b3a2"
var_col <- "goldenrod"  # "#c7909d"


X11(width = 6, height = 5)
lm_result <- ggplot(str_var, aes(x = wing.loading..Nm.2.)) +
  geom_smooth(aes(y = max_wind), method = "lm", color = max_col, alpha = .2, fill = max_col) +
  geom_smooth(aes(y = max_str_cov), method = "lm", color = var_col, alpha = .2, fill = var_col) +
  geom_point(aes(y = max_wind), color = max_col, alpha = .6) +
  geom_point(aes(y = max_str_cov), color = var_col, alpha = .6) +
  labs(x = "Wind loading") +
  scale_y_continuous(
    name = "Maximum wind speed (m/s)",# Features of the first axis
    sec.axis = sec_axis(~ . + 10,name = "Variation in wind speed (%)")) + # Add a second axis and specify its features
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(color = max_col, size=13),
        axis.title.y.right = element_text(color = var_col, size=13)) 


png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/lm_output.png", width = 6, height = 5, units = "in", res = 300)
print(lm_result)
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



