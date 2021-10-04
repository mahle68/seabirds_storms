#script for analysis and figure prep for the seabird manuscript
#Sep 30. 2021. Elham Nourani, PhD. Radolfzell, DE.
#sources: null_modeling_seabirds.R and PCoA_seabirds.R

library(tidyverse)
library(parallel)
library(hrbrthemes)
library(ggridges)
library(corrr)
library(sp)

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")

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


png("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/Cov_bar.png", width = 5, height = 5, units = "in", res = 300)
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


png("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/sig_plots.png", width = 7, height = 4, units = "in", res = 300)
print(sig_plots)
dev.off()

### STEP 5: Raw wind plots ------------------------------------ #####

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

#with quartiles
X11(width = 8, height = 7)
raw_wind <- ggplot(ann, aes(x = wind_speed_ms, y = group, fill = factor(stat(quantile)))) + 
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = F, scale = 3,
  ) +
  scale_fill_viridis_d(name = "Quartiles", alpha = 0.6) +
  scale_x_continuous(limits = c(0, 25)) +
  labs(y = "", x = "Wind speed (m/s)") +
  theme_minimal() #+
theme(legend.position = "bottom") 

png("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/raw_wind.png", width = 8, height = 7, units = "in", res = 300)
print(raw_wind)
dev.off()

### STEP 6: LM: wind strength ------------------------------------ #####

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