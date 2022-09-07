#script for reproducing the results and figures of the manuscript "Maximum tolerable wind speed varies with morphology in seabirds"
#April 2022. Elham Nourani, PhD. Konstanz, Germany.
#update Aug 15: add wind conditions at distribution range for each species. problems: 1) wind is for 2017-2022 while gps is for 2005-2019; 2) only one data file per speceis,
#hence ignoring the 

#load libraries
library(tidyverse)
library(ggridges)
library(parallel)
library(hrbrthemes)
library(corrr)
library(sp)
library(sf)
library(patchwork)
library(lubridate)
library(oce)
library(gridExtra)
library(ape)
library(stargazer)
library(corrr)
library(introdataviz) #devtools::install_github("psyteachr/introdataviz")

### STEP 1: Open all input data ------------------------------------ #####

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/used_alt_annotated.RData") #ann_30; this file contains used and alternative steps (hourly). This will be used for permutation tests.
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/species_summary_data.RData") #lm_input; this file contains summary information for each species including morphology and wind conditions experienced
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/hrly_data_ann.RData") #ann; this file contains sub-sampled hourly data annotated with atmospheric information
wind_distr <- readRDS("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/wind_at_distr_17_22.rds")

#make sure species names, etc match between my data and Manos'
sp_flight <- lm_input %>% 
  mutate(flight_style_F = as.factor(flight_style)) %>% 
  dplyr::select(c("sci_name", "species", "flight_style_F")) %>% 
  group_by(species) %>% 
  slice(1)

range_data <- wind_distr %>% 
  dplyr::select(-c("Group", "species")) %>%  
  mutate(wind_data = "range") %>% 
  rename(wind_speed_ms = WindSpeed, 
         sci_name = Sci_name,
         year = Year) %>% 
  left_join(sp_flight, by = "sci_name") %>%  #add species name and flight type
  mutate()


### add range wind info to lm_input

range_summ <- range_data %>% 
  group_by(species) %>% 
  summarize(range_mean = mean(wind_speed_ms),
            range_median = median(wind_speed_ms),
            range_max = max(wind_speed_ms))
  
lm_input <- lm_input %>% 
  full_join(range_summ, by = "species")
  
saveRDS(lm_input, "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/species_summary_data_with_ranges.RDS")

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

### Plot Fig. S3 -------------------------------------------------------------------------

X11(width = 5, height = 5)
CoV_bar <- ggplot(data_var, aes(x = wspd_cov, y = reorder(as.factor(species), desc(as.factor(species))), height = stat(density))) + 
  geom_density_ridges(
    stat = "binline", bins = 20, scale = 0.98, alpha = 0.3,
    draw_baseline = FALSE
  ) +
  scale_x_continuous(limits = c(-4, 140)) +
  labs(y = "", x = "Coefficient of variation (%)") +
  theme_minimal()
 
 
### STEP 3: Permutation test ------------------------------------ #####
#In each stratum, is the difference between selected and max available wind speed higher than expected by chance?
#for each stratum, calculate the difference between observed wind speed and max wind speed (incl. observed)

observed_stat <- ann_30 %>% 
  group_by(group, year, stratum) %>% #group by species-colony (group) instead of just species
  arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached detach("package:plyr", unload=TRUE)
  summarize(max_minus_obs = max(wind_speed) - head(wind_speed,1))

#randomize and calculate the same statistic
#shuffle all wind speed values within each year, then recalculate the statistic within each stratum

permutations <- 1000

#prep cluster
mycl <- makeCluster(detectCores() - 2, setup_strategy = "sequential")
clusterExport(mycl, c("permutations", "ann_30")) 

clusterEvalQ(mycl, {
  library(dplyr)
})

rnd_stat <- parLapply(cl = mycl, X = c(1:permutations), fun = function(x){ #if not using multiple cores, call lapply instead of ParLapply
  
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


stopCluster(mycl)


#extract observed and random values for each stratum

#prep cluster
mycl <- makeCluster(detectCores() - 7, setup_strategy = "sequential")
clusterExport(mycl, c("permutations", "observed_stat", "rnd_stat")) 

clusterEvalQ(mycl, {
  library(dplyr)
})

p_vals <- parLapply(mycl, unique(observed_stat$stratum), function(x){
  obs <- observed_stat[observed_stat$stratum == x,]
  rnd <- rnd_stat[rnd_stat$stratum == x,]
  
  obs$p_less <- sum(rnd$max_minus_obs <= obs$max_minus_obs)/permutations
  obs$p_more <- sum(rnd$max_minus_obs >= obs$max_minus_obs)/permutations
  
  obs
}) %>% 
  reduce(rbind)

stopCluster(mycl)


# P_vals is available as "p_vals_1000_perm_df.RData"


### STEP 4: Plot permutation results ------------------------------------ #####
 
load("p_vals_1000_perm_df.RData") #p_vals

### Plot Fig. S4 -------------------------------------------------------------------------
X11(width = 5, height = 5)

perm_sig <- ggplot(p_vals, aes(x = p_more, y = reorder(as.factor(species), desc(as.factor(species))), height = stat(density))) + 
  geom_density_ridges(
   stat = "binline", bins = 20, scale = 0.98, alpha = 0.3,
    draw_baseline = FALSE
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "", x = "Significance") +
  theme_minimal()


### Plot Fig. S5 -------------------------------------------------------------------------
#extract strata with significant avoidance of strong winds
#extract strata with p-value less than 0.05
sig <- p_vals %>%
  filter(p_more <= 0.05)

sig_data <- ann_30 %>%
  filter(stratum %in% sig$stratum) %>% 
  mutate(species = fct_relevel(as.factor(common_name), levels = "Atlantic yellow-nosed albatross", "Wandering albatross", 
                                 "Sooty albatross", "Red-footed booby")) %>% 
  mutate(col = as.character(fct_relevel(species), levels = "corn flower blue", "rosy brown", "yellow green", "pale violet red")) %>% 
  as.data.frame()

sig_data$common_name <- reorder(sig_data$common_name, sig_data$species)
sig_data$stratum <- reorder(as.factor(sig_data$stratum),desc(sig_data$species))

#plot
used_wind <- sig_data %>% 
  group_by(stratum) %>% 
  filter(used ==1) %>% 
  ungroup() %>% 
  dplyr::select(c("stratum", "common_name", "wind_speed")) %>% 
  as.data.frame()


X11(width = 7, height = 4)
sig_plots<- ggplot(sig_data, aes(x = wind_speed, y = stratum)) + 
  geom_density_ridges(scale = 3, alpha = 0.4,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 2, point_alpha = 1, aes(fill = species)) + 
  scale_fill_manual(values = c("Atlantic yellow-nosed albatross" = "corn flower blue", "Wandering albatross" = "yellowgreen", 
                               "Sooty albatross" = "lightcoral", "Red-footed booby" = "goldenrod")) +
  scale_x_continuous(limits = c(0, 25)) +
  geom_segment(data = used_wind, aes(x = wind_speed, xend = wind_speed, y = as.numeric(as.factor(stratum)),
                                     yend = as.numeric(as.factor(stratum)) + .9, linetype = "Selected wind speed"), color = "red") +
  guides(fill = guide_legend(order = 1),
         line = guide_legend(order = 2)) +
  labs(y = "Density", x = expression("Wind speed (m s"^-1*")")) +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank())

### STEP 5: Linear Model: wind strength ------------------------------------ #####

#try one model with more stuff
lm_input <- readRDS("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/species_summary_data_with_ranges.RDS")

lm_input %>% 
  dplyr::select(c("max_wind_ms", "wing.loading..Nm.2.", "range_max", "range_median")) %>% 
  correlate()  #median_range and wing loading are identical. 97% correlated
  
range_M <- lm(max_wind_ms ~ range_max + range_median, data = lm_input)
morph_M <- lm(max_wind_ms ~ wing.loading..Nm.2. + range_max , data = lm_input)

#calculate correlation between wing loading and aspect ratio
cor.test(lm_input$wing.loading..Nm.2., lm_input$aspect.ratio) #0.60

#model
morph <- lm(max_wind_ms ~ wing.loading..Nm.2., data = lm_input) #adR = 0.31 , AIC =  119.2178


#get latex output (Table S2)
stargazer(morph)

#investigate residuals
par(mfrow = c(2,2))
plot(morph) 

#temporal correlation. result: no temporal autocorrelation
acf(resid(morph))
acf(resid(morph),type = "p")

#spatial autocorrelation. bubble plot result: no autocorrelation
spdata <- data.frame(resid = resid(morph), x = lm_input$colony.long, y = lm_input$colony.lat)
coordinates(spdata)<-~ x + y
x11()
bubble(spdata, "resid", col = c("blue","orange"))

#check for phylogeny
#Phyl tree downloaded from https://birdtree.org/subsets/ for the 18 species
trees <- read.nexus("phylogeny_tree/output.nex")
#use the first tree
plot(trees[[1]])
axisPhylo()
w <- 1/cophenetic(trees[[1]]) #Computes the cophenetic distances for a hierarchical clustering.
diag(w) <- 0 #set the diagonal to 0

### Plot Fig. S6 -------------------------------------------------------------------------
#save the tree as supplementary material
X11(width = 6, height = 7)
plot(trees[[1]])

#extract wind speed
#change sci name for great shearwater to match the sci name in the phyl. tree
lm_input[lm_input$sci_name == "Ardenna gravis", "sci_name"] <- "Puffinus gravis"
#take average of the two colonies for species that have multiple rows
max_wspd_unique <- lm_input %>% 
  group_by(sci_name) %>% 
  summarize(max_wind_ms = mean(max_wind_ms)) %>% 
  ungroup() 
max_wspd <- max_wspd_unique$max_wind_ms
names(max_wspd) <- max_wspd_unique$sci_name

#estimate Moran's I for maximum wind speed
Moran.I(max_wspd, w)

### STEP 6: Linear Model: wind variability ------------------------------------ #####

#use data_var created in STEP 2

summary_info <- ann %>% 
  group_by(sci_name, colony.name) %>% 
  summarize(n_rows = n(),
            n_trips = n_distinct(TripID)) %>% 
  full_join(lm_input, by = c("sci_name", "colony.name")) %>% 
  mutate(max_wind_ms = max_wind/3.6,
         group = paste(species, colony.name, sep = "_")) %>% 
  as.data.frame()

str_var <- data_var %>% 
  group_by(group) %>% 
  summarize(max_str_cov = max(wspd_cov)) %>% 
  full_join(summary_info, by = "group") %>% 
  as.data.frame()

m1 <- lm(max_str_cov ~ wing.loading..Nm.2., data = str_var) #AIC = 188.3322; adjRsq = 0.3534  

#get latex output (Table S3)
stargazer(m1)

#estimate moran's I

#change sci name for great shearwater to match the sci name in the tress
str_var[str_var$sci_name == "Ardenna gravis", "sci_name"] <- "Puffinus gravis"
#take average of the two colonies for species that have multiple rows
var_wspd_unique <- str_var %>% 
  group_by(sci_name) %>% 
  summarize(max_str_cov = mean(max_str_cov, na.rm = T)) %>% 
  ungroup() 
var_wspd <- var_wspd_unique$max_str_cov
names(var_wspd) <- var_wspd_unique$sci_name

#estimate Moran's I for wind variability
Moran.I(var_wspd, w)

### Plot Fig. 3 -------------------------------------------------------------------------

#extract a color from the oce palette for cohesion
clr <- oce::oceColorsPalette(120)[14]

str_var$flight_style_F <- factor(str_var$flight_style)

X11(width = 11, height = 5)
#plot for linear model predicting wind speed
lm_maxwind <- ggplot(str_var, aes(x = wing.loading..Nm.2.)) +
  geom_smooth(aes(y = max_wind_ms), method = "lm", color = clr, alpha = .1, fill = clr) +
  geom_point(aes(y = max_wind_ms, shape = flight_style_F), size = 2, stroke = 0.8, color = clr) +
  labs(x = expression("Wing loading (Nm"^-2*")")) +
  scale_shape_manual(values = c(4,0,2,1)) + #filled points: c(15, 17, 19)
  scale_y_continuous(
    name = expression("Maximum wind speed (m s"^-1*")")) +# Features of the first axis
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(size = 13)) +
  guides(shape = guide_legend("Flight style:"))

#plot for linear model predicting wind covariance
lm_covwind <- ggplot(str_var, aes(x = wing.loading..Nm.2.)) +
  geom_smooth(aes(y = max_str_cov), method = "lm", color = clr, alpha = .1, fill = clr) + #confidence intervals of the 0.95% by default
  geom_point(aes(y = max_str_cov, shape = flight_style_F), size = 2, stroke = 0.8,  color = clr) +
  labs(x = expression("Wing loading (Nm"^-2*")")) +
  scale_shape_manual(values = c(4,0,2,1)) + #filled points: c(15, 17, 19)
  scale_y_continuous(
    name = "Variation in wind speed (%)") +# Features of the first axis
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(size = 13))+
  guides(shape = guide_legend("Flight style:"))

#combine plots
combined <- lm_maxwind + lm_covwind & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")

#update Jun 17: plot body mass to respond to reviewer
ggplot(str_var, aes(x = body.mass..kg.)) +
  geom_smooth(aes(y = max_wind_ms), method = "lm", color = clr, alpha = .1, fill = clr) +
  geom_point(aes(y = max_wind_ms, shape = flight_style_F), size = 2, stroke = 0.8, color = clr) +
  labs(x = "Body mass (kg)") +
  scale_shape_manual(values = c(4,0,2,1)) + #filled points: c(15, 17, 19)
  scale_y_continuous(
    name = expression("Maximum wind speed (m s"^-1*")")) +# Features of the first axis
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(size = 13)) +
  guides(shape = guide_legend("Flight style:"))

#remove heavy albatrosses
ggplot(str_var %>% filter(!(species %in% c("Wandering albatross", "Tristan albatross"))), aes(x = body.mass..kg.)) +
  geom_smooth(aes(y = max_wind_ms), method = "lm", color = clr, alpha = .1, fill = clr) +
  geom_point(aes(y = max_wind_ms, shape = flight_style_F), size = 2, stroke = 0.8, color = clr) +
  labs(x = "Body mass (kg)") +
  scale_shape_manual(values = c(4,0,2,1)) + #filled points: c(15, 17, 19)
  scale_y_continuous(
    name = expression("Maximum wind speed (m s"^-1*")")) +# Features of the first axis
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(size = 13)) +
  guides(shape = guide_legend("Flight style:"))


### STEP 7: Correlation between data quantity and max wind speeds ------------------------------------ #####

#are nrows and n trips correlated?
cor(summary_info$n_trips,summary_info$n_rows) #yes! 0.88

#are max wind speeds correlated with nrows?
cor.test(summary_info$n_trips, summary_info$max_wind_ms) #0.15
cor.test(summary_info[, "n_trips"], summary_info[,"max_wind_ms"]) #-0.3

cor.test(str_var$n_trips, str_var$max_str_cov) #0.32
cor.test(str_var[,"n_trips"], str_var[,"max_str_cov"]) # -0.15


### STEP 8: Raw wind plots ------------------------------------ #####

### Plot Fig. 1 -------------------------------------------------------------------------

ann_ft <- ann %>% 
  left_join(str_var[,c("sci_name","flight_style_F")])

#merge the data from Manos (i.e. wind at distribution range) and tracking data annotated
#with wind.

sp_flight <- lm_input %>% 
  mutate(flight_style_F = as.factor(flight_style)) %>% 
  dplyr::select(c("sci_name", "species", "flight_style_F")) %>% 
  group_by(species) %>% 
  slice(1)

range_data <- wind_distr %>% 
  dplyr::select(-c("Group", "species")) %>%  
  mutate(wind_data = "range") %>% 
  rename(wind_speed_ms = WindSpeed, 
         sci_name = Sci_name,
         year = Year) %>% 
  left_join(sp_flight, by = "sci_name") %>%  #add species name and flight type
  mutate()

wind_raw <- ann_ft %>% 
  mutate(year = year(timestamp),
         wind_data = "gps_pts") %>% 
  dplyr::select(c("wind_speed_ms", "sci_name", "year", "species", "wind_data",  "flight_style_F")) %>% 
  bind_rows(range_data)
  
save(wind_raw, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/wind_raw.rdata")

cols <- oce::oceColorsPalette(10)
#extract a color from the oce palette for cohesion
clr <- oce::oceColorsPalette(120)[14]

X11(width = 8, height = 7.5)
raw_wind <- ggplot(wind_raw, aes(x = wind_speed_ms, y = species, fill = "wind_data", color = "wind_data")) + 
    stat_density_ridges(jittered_points = TRUE, rel_min_height = .01,
                        point_shape = "|", point_size = 0.8, point_alpha = 0.5, size = 0.25) + #,
                        #geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, #only the median line
                        #quantiles = 0.5, quantile_lines = T, scale = 3) +
    geom_point(data = wind_raw[wind_raw$wind_data == "gps_pts",], aes(x = -0.8, y = species, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
    scale_x_continuous(limits = c(-0.8, 32)) +
    #scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
    #                     na.value = "white", guide = 'none') +
    scale_shape_manual(values = c(4,0,2,1)) +
    labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
    theme_minimal() +
    guides(shape = guide_legend("Flight style:")) +
    theme(legend.position = "bottom",legend.title = element_text(size = 10), 
          legend.text=element_text(size = 7))



raw_windj <- ggplot(wind_raw, aes(x = wind_speed_ms, y = species, fill = wind_data)) + 
  geom_density_ridges(alpha = 0.5, jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7, color = wind_data)

png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/wind_density.png", width = 8, height = 8, units = "in", res = 300)
ggplot(wind_raw, aes(x = wind_speed_ms, y = species, fill = wind_data)) + 
  geom_density_ridges(alpha = 0.5) +
  theme_bw()
dev.off()

#test plot with jitters
png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/test.png", width = 8, height = 8, units = "in", res = 300)
ggplot(sample, aes(x = wind_speed_ms, y = species, fill = "wind_data", color = "wind_data")) + 
  stat_density_ridges(jittered_points = TRUE, rel_min_height = .01,
                      point_shape = "|", point_size = 0.8, point_alpha = 0.5, size = 0.25) +
  geom_point(data = wind_raw[wind_raw$wind_data == "gps_pts",], aes(x = -0.8, y = species, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
  scale_x_continuous(limits = c(-0.8, 31.5)) +
  #scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
  #                     na.value = "white", guide = 'none') +
  scale_shape_manual(values = c(4,0,2,1)) +
  labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
  theme_minimal() +
  guides(shape = guide_legend("Flight style:")) +
  theme(legend.position = "bottom",legend.title = element_text(size = 10), 
        legend.text=element_text(size = 7))

dev.off()

#boxplots
png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/wind_box.png", width = 8, height = 8, units = "in", res = 300)
ggplot(wind_raw, aes(x = wind_speed_ms, y = species, fill = wind_data)) + 
  geom_boxplot(outlier.size = 0.5) +
  theme_bw()
dev.off()

#try split-violin plots
#https://psyteachr.github.io/introdataviz/advanced-plots.html

#work on a sample
sample <- wind_raw %>% 
  filter(species %in% c("Wandering albatross" , "Nazca booby", "Red-tailed tropicbird" ))
sample$species_f <- as.factor(sample$species)
sample$wind_speed_ms <- round(sample$wind_speed_ms, digits = 2)

ggplot(sample, aes(y = wind_speed_ms, x = species_f, fill = wind_data)) +
  introdataviz::geom_split_violin(alpha = .4, trim = T) +
  #geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
  #stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
  #             position = position_dodge(.175)) +
  scale_x_discrete(name = "Condition", labels = c("Non-word", "Word")) +
  scale_y_continuous(name = "Reaction time (ms)",
                     breaks = seq(200, 800, 100), 
                     limits = c(200, 800)) +
  scale_fill_brewer(palette = "Dark2", name = "Language group") +
  theme_minimal()

#############
#compare max wind speeds in gps vs data pts
maxx <- wind_raw %>% group_by(species, wind_data) %>% summarize(max_wind = max(wind_speed_ms))

plot(as.factor(maxx$species), maxx$max_wind, col = as.factor(maxx$wind_data))

# X11(width = 8, height = 7.5)
# (raw_wind <- ggplot(wind_raw, aes(x = wind_speed_ms, y = species, fill = wind_data)) + 
#   stat_density_ridges(jittered_points = TRUE, rel_min_height = .01,
#                       point_shape = "|", point_size = 0.8, point_alpha = 0.5, size = 0.25,
#                       geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, #only the median line
#                       quantiles = 0.5, quantile_lines = T, scale = 3) +
#   geom_point(data = wind_raw[wind_raw$wind_data == "gps_pts",], aes(x = -0.8, y = species, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
#   scale_x_continuous(limits = c(-0.8, 30)) +
#   #scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
#   #                     na.value = "white", guide = 'none') +
#   scale_shape_manual(values = c(4,0,2,1)) +
#   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
#   theme_minimal() +
#   guides(shape = guide_legend("Flight style:")) +
#   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
#         legend.text=element_text(size = 7)))

### Fig 2 -------------------------------------------------------------------------
#based on https://semba-blog.netlify.app/10/29/2018/animating-oceanographic-data-in-r-with-ggplot2-and-gganimate/