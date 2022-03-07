#script for reproducing the results and figures of the manuscript "Maximum tolerable wind speed varies with morphology in seabirds"
#February 2022. Elham Nourani, PhD. Konstanz, Germany.

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

#define functions and coordinate reference system
#setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")
#source("/home/mahle68/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")

### STEP 1: Open and prep SSF data ------------------------------------ #####
#open dataset with alternative steps for hourly steps. prepped in random_steps.R
load("R_files/ssf_input_annotated_60_15_30alt_18spp.RData") #ann_30
  
ann_30 <- ann_30 %>%  
  mutate(group = paste(common_name, colony.name, sep = "_")) %>% 
  as.data.frame()

#fix the common name capitalization inconsistencies
# ann_30[ann_30$common_name == "Nazca boobies", "common_name"] <- "Nazca booby"
# ann_30[ann_30$common_name == "Tristan Albatross", "common_name"] <- "Tristan albatross"
# ann_30[ann_30$common_name == "Sooty Albatross", "common_name"] <- "Sooty albatross"
# ann_30[ann_30$common_name == "Great Shearwater", "common_name"] <- "Great shearwater"
# ann_30[ann_30$common_name == "Grey Petrel", "common_name"] <- "Grey petrel"
# ann_30[ann_30$common_name == "Galapagos albatross (waved albatross)", "common_name"] <- "Galapagos albatross"
# ann_30[ann_30$common_name == "Soft-plumaged Petrel", "common_name"] <- "Soft-plumaged petrel"
# ann_30[ann_30$common_name == "Atlantic Petrel", "common_name"] <- "Atlantic petrel"
# ann_30[ann_30$common_name == "Atlantic Yellow-nosed Albatross", "common_name"] <- "Atlantic yellow-nosed albatross"
# ann_30[ann_30$common_name == "masked booby", "common_name"] <- "Masked booby"

#write to public repo
save(ann_30, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/seabirds_storms_public/used_alt_annotated.RData")


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
 
 png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/Cov_bar.png", width = 5, height = 5, units = "in", res = 300)
 print(CoV_bar)
 dev.off()
 
### STEP 3: Permutation test ------------------------------------ #####
#In each stratum, is the difference between selected and max available wind speed higher than expected by chance?

#for each stratum, calculate the difference between observed wind speed and max wind speed (incl. observed)

observed_stat <- ann_30 %>% 
  group_by(group, year, stratum) %>% #group by species-colony instead of just species
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


#save(p_vals, file = "R_files/p_vals_1000_perm_df.RData")


### STEP 4: Plot permutation results ------------------------------------ #####
 
#load("R_files/p_vals_1000_perm_df.RData") #p_vals

p_vals <- p_vals %>% 
  mutate(year_f = as.character(year),
         species = str_split(group, "\\_", simplify = TRUE)[, 1]) %>% 
  ungroup() %>% 
  as.data.frame()

p_vals[p_vals$species == "Galapagos albatross (waved albatross)", "species"] <- "Galapagos albatross"

#fix the common name capitalization inconsistencies
p_vals[p_vals$species == "Nazca boobies", "species"] <- "Nazca booby"
p_vals[p_vals$species == "Tristan Albatross", "species"] <- "Tristan albatross"
p_vals[p_vals$species == "Sooty Albatross", "species"] <- "Sooty albatross"
p_vals[p_vals$species == "Great Shearwater", "species"] <- "Great shearwater"
p_vals[p_vals$species == "Grey Petrel", "species"] <- "Grey petrel"
p_vals[p_vals$species == "Galapagos albatross (waved albatross)", "species"] <- "Galapagos albatross"
p_vals[p_vals$species == "Soft-plumaged Petrel", "species"] <- "Soft-plumaged petrel"
p_vals[p_vals$species == "Atlantic Petrel", "species"] <- "Atlantic petrel"
p_vals[p_vals$species == "Atlantic Yellow-nosed Albatross", "species"] <- "Atlantic yellow-nosed albatross"
p_vals[p_vals$species == "masked booby", "species"] <- "Masked booby"



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

png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/Perm_sig.png", width = 5, height = 5, units = "in", res = 300)
print(perm_sig)
dev.off()



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
  
  #scale_linetype_manual("",values = c("Selected wind speed" = 1)) +
  guides(fill = guide_legend(order = 1),
         line = guide_legend(order = 2)) +
  labs(y = "Density", x = expression("Wind speed (m s"^-1*")")) +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank())

png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/sig_plots_unit_fixed.png", width = 7, height = 4, units = "in", res = 300)
print(sig_plots)
dev.off()



### STEP 5: Raw wind plots ------------------------------------ #####

#open annotated data and species data
load("R_files/lm_input_20spp_col.RData") #lm_input 
load("R_files/ann_18spp.RData") #ann (from PCoA_seabirds.R)

#save in public repo
save(ann, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/seabirds_storms_public/hrly_data_ann.RData")


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
  

### Plot Fig. 1 -------------------------------------------------------------------------

cols <- oce::oceColorsPalette(10)

X11(width = 8, height = 7)
raw_wind <- ggplot(ann, aes(x = wind_speed_ms, y = group_f, fill = stat(x))) + 
  stat_density_ridges(jittered_points = TRUE, rel_min_height = .01,
                      point_shape = "|", point_size = 0.8, point_alpha = 0.5, size = 0.25,
                      geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, #only the median line
                      quantiles = 0.5, quantile_lines = T, scale = 3) +
  scale_x_continuous(limits = c(-0.5, 25)) +
  scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                       na.value = "white") +
  labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
  theme_minimal() +
  theme(legend.position = "none")

png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/raw_wind_jitter_desc_all18_nbfixed.png", width = 8, height = 8, units = "in", res = 300)
print(raw_wind)
dev.off()

### STEP 6: Linear Model: wind strength ------------------------------------ #####

load("R_files/lm_input_20spp_col.RData") #lm_input #update 7.2.2022: wing loading and aspect ratio for nazca booby were corrected usign Jennifer McKee's paper

# lm_input[lm_input$species == "Nazca booby", "wing.loading..Nm.2."] <- 73.64
# lm_input[lm_input$species == "Nazca booby", "aspect.ratio"] <- 7.10
# 
# #also fix the species names. all sentence case
# lm_input[lm_input$species == "Nazca boobies", "species"] <- "Nazca booby"
# lm_input[lm_input$species == "Tristan Albatross", "species"] <- "Tristan albatross"
# lm_input[lm_input$species == "Sooty Albatross", "species"] <- "Sooty albatross"
# lm_input[lm_input$species == "Great Shearwater", "species"] <- "Great shearwater"
# lm_input[lm_input$species == "Grey Petrel", "species"] <- "Grey petrel"
# lm_input[lm_input$species == "Galapagos albatross (waved albatross)", "species"] <- "Galapagos albatross"
# lm_input[lm_input$species == "Soft-plumaged Petrel", "species"] <- "Soft-plumaged petrel"
# lm_input[lm_input$species == "Atlantic Petrel", "species"] <- "Atlantic petrel"
# lm_input[lm_input$species == "Atlantic Yellow-nosed Albatross", "species"] <- "Atlantic yellow-nosed albatross"
# lm_input[lm_input$species == "masked booby", "species"] <- "Masked booby"


#save(lm_input, file = "R_files/lm_input_20spp_col.RData")

lm_input <- lm_input %>% 
  #mutate_at(c("colony.long","colony.lat", "wing.loading..Nm.2.", "wing.span..m.", "wing.area..m2.","aspect.ratio", "PC1", "PC2"),
  #          list(z = ~as.numeric(scale(.)))) %>% 
  mutate(breeding_ch = as.character(median_breeding_m),
         max_wind_ms = max_wind/3.6) %>% 
  arrange(median_breeding_yday) #order temporally, for temporal autocorrelation analysis.
  

#save in public repo
lm_input <- lm_input %>% 
  dplyr::select(-c("PC1", "PC2")) %>% 
  as.data.frame()
save(lm_input, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/seabirds_storms_public/species_summary_data.R") #updated 7.3.2022

#calculate correlation between wing loading and aspect ratio
cor.test(lm_input$wing.loading..Nm.2., lm_input$aspect.ratio) #0.60

#model
morph <- lm(max_wind_ms ~ wing.loading..Nm.2., data = lm_input) #adR = 0.31 , AIC =  119.2178
save(morph, file = "max_wspd_model.RData")

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
bubble(spdata, "resid", col = c("blue","orange"))

#check for phylogeny
#Phyl tree downloaded from https://birdtree.org/subsets/ for the 18 species
trees <- read.nexus("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/phylogeny_tree/tree-pruner-5edeac95-d8ed-410a-8162-982b1438f4e9/output.nex")
#use the first tree
plot(trees[[1]])
axisPhylo()
w <- 1/cophenetic(trees[[1]]) #Computes the cophenetic distances for a hierarchical clustering.
diag(w) <- 0 #set the diagonal to 0

### Plot Fig S6 -------------------------------------------------------------------------
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

### STEP 7: Linear Model: wind variability ------------------------------------ #####

#use data_var created in STEP 2 and summary_info created in STEP 6
summary_info <- summary_info %>%  
  mutate(group = paste(species, colony.name, sep = "_"))

str_var <- data_var %>% 
  group_by(group) %>% 
  summarize(max_str_cov = max(wspd_cov)) %>% 
  full_join(summary_info, by = "group") %>% 
  as.data.frame()

m1 <- lm(max_str_cov ~ wing.loading..Nm.2., data = str_var) #AIC = 189.2436; adjRsq = 0.3216  

#get latex output (Table S3)
stargazer(morph)

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

### Plot Fig 3 -------------------------------------------------------------------------

#extract a color from the oce palette for cohesion
clr <- oce::oceColorsPalette(120)[14]

#str_var was created in STEP 7
str_var[str_var$species == "Red-tailed tropicbird","flight.type"] <- "flap-gliding"
str_var[str_var$species == "Great Shearwater","flight.type"] <- "dynamic soaring"
str_var$flight.type_F <- factor(str_var$flight.type)

X11(width = 11, height = 5)
#plot for linear model predicting wind speed
lm_maxwind <- ggplot(str_var, aes(x = wing.loading..Nm.2.)) +
  geom_smooth(aes(y = max_wind_ms), method = "lm", color = clr, alpha = .1, fill = clr) +
  geom_point(aes(y = max_wind_ms, shape = flight.type_F), size = 2, stroke = 0.8, color = clr) +
  labs(x = expression("Wing loading (Nm"^-2*")")) +
  scale_shape_manual(values = c(0,2,1)) + #filled points: c(15, 17, 19)
  scale_y_continuous(
    name = expression("Maximum wind speed (m s"^-1*")")) +# Features of the first axis
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(size = 13)) +
  guides(shape = guide_legend("Flight type:"))

#plot for linear model predicting wind covariance
lm_covwind <- ggplot(str_var, aes(x = wing.loading..Nm.2.)) +
  geom_smooth(aes(y = max_str_cov), method = "lm", color = clr, alpha = .1, fill = clr) +
  geom_point(aes(y = max_str_cov, shape = flight.type_F), size = 2, stroke = 0.8,  color = clr) +
  labs(x = expression("Wing loading (Nm"^-2*")")) +
  scale_shape_manual(values = c(0,2,1)) + #filled points: c(15, 17, 19)
  scale_y_continuous(
    name = "Variation in wind speed (%)") +# Features of the first axis
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(size = 13))+
  guides(shape = guide_legend("Flight type:"))

#combine plots
combined <- lm_maxwind + lm_covwind & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")

### STEP 8: Correlation between data quantity and max wind speeds ------------------------------------ #####

#load annotated data (one hourly subset)
load("R_files/ann_18spp.RData") #ann; from PCoA_seabirds.R
load("R_files/lm_input_20spp_col.RData") #lm_input 

#extract number of data points and number of trips
summary_info <- ann %>% 
  group_by(sci_name, colony.name) %>% 
  summarize(n_rows = n(),
            n_trips = n_distinct(TripID)) %>% 
  full_join(lm_input, by = c("sci_name", "colony.name")) %>% 
  mutate(max_wind_ms = max_wind/3.6) %>% 
  as.data.frame()

#are nrows and n trips correlated?
cor(summary_info$n_trips,summary_info$n_rows) #yes! 0.88

#are max wind speeds correlated with nrows?
cor.test(summary_info$n_trips, summary_info$max_wind_ms) #0.15
cor.test(summary_info[, "n_trips"], summary_info[,"max_wind_ms"]) #-0.3

cor.test(str_var$n_trips, str_var$max_str_cov) #0.32
cor.test(str_var[,"n_trips"], str_var[,"max_str_cov"]) # -0.15

### Fig 2 -------------------------------------------------------------------------
#based on https://semba-blog.netlify.app/10/29/2018/animating-oceanographic-data-in-r-with-ggplot2-and-gganimate/