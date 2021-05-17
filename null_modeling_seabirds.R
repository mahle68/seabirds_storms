#script for null modeling of seabird data
#May 11. 2021. Elham Nourani, PhD. Radolfzell, DE.
#https://www.jwilber.me/permutationtest/#:~:text=To%20calculate%20the%20p%2Dvalue,of%20test%2Dstatistics%20we%20calculated


library(tidyverse)
library(cowplot) #raincloud plot
library(parallel)
detach(package:plyr)
library(sm) #for one density plot per factor level

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")


source("/home/enourani/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/R_rainclouds.R")
source("/home/enourani/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/summarySE.R")
source("/home/enourani/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/simulateData.R")

#open dataset with alternative steps for hourly steps. prepped in all_species_ssf_2021.R
load("R_files/ssf_input_annotated_60_30_40alt_18spp.RData") #ann_40

#summary details
ann_40 %>% 
  filter(used == 1) %>% 
  group_by(common_name) %>% 
  summarize(n_tracks = n_distinct(TripID),
            n_ind = n_distinct(indID),
            n_yr = n_distinct(year),
            min_yr = min(year),
            max_yr = max(year))

# Step 1: calc within-stratum variances ####

#calculate and plot within stratum variances
data_var <- ann_40 %>%
  group_by(stratum) %>%
  summarise(wspd_var = var(wind_speed),
            u_var = var(u10m),
            v_var = var(v10m),
            wspt_var = var(wind_support),
            species = head(common_name, 1),
            year = head(year, 1))


#rain cloud plots for variances
#restructure the dataframe
wspd_var <- data_var %>% 
  dplyr::select(-c(u_var,v_var,wspt_var)) %>% 
  mutate(variable = "wind speed") %>% 
  dplyr::rename(score = wspd_var)

wspt_var <- data_var %>% 
  dplyr::select(-c(u_var,v_var,wspd_var)) %>% 
  mutate(variable = "wind support") %>% 
  dplyr::rename(score = wspt_var)

w_u_var <- data_var %>% 
  dplyr::select(-c(v_var,wspd_var, wspt_var)) %>% 
  mutate(variable = "u wind") %>% 
  dplyr::rename(score = u_var)

w_v_var <- data_var %>% 
  dplyr::select(-c(u_var,wspd_var, wspt_var)) %>% 
  mutate(variable = "v wind") %>% 
  dplyr::rename(score = v_var)


new_data_var <- rbind(wspd_var, wspt_var, w_u_var, w_v_var)
new_data_var <- rbind(w_u_var, w_v_var)

X11()
plot_variances <- ggplot(new_data_var, aes(x = variable, y = score, fill = species)) +
  ylim(0,20) +
  geom_flat_violin(aes(fill = species),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = species),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = species),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) #+
  #facet_grid(season~.)

ggsave("rain_cloud_plot_variances.tiff",plot = plot_variances, dpi = 500,
       path = "/home/enourani/ownCloud/Work/safi_lab_meeting/presentation_jan17")


# Step 2: permutation test- In each stratum, is the difference between selected and available wind speed higher/lower than expected by chance? ####

#for each stratum, calculate the difference between observed wind speed and max wind speed (incl. observed)

ann_40 <- ann_40 %>% 
dplyr::select(c(30,2,17,14,3:13,15,16,18:29,31:ncol(ann_40))) 

observed_stat <- ann_40 %>% 
  group_by(common_name, year, stratum) %>% 
  arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached
  summarize(max_minus_obs = max(wind_speed) - head(wind_speed,1))


#randomize and calculate the same statistic
#shuffle all wind speed values within each year, then recalc the statistic within each stratum

permutations <- 100

#prep cluster
  mycl <- makeCluster(detectCores() - 3)
  clusterExport(mycl, c("permutations", "ann_40")) 
  
  clusterEvalQ(mycl, {
    library(dplyr)
  })
  
  
  a <- Sys.time()
  
  rnd_stat <- parLapply(cl = mycl, X = c(1:permutations), fun = function(x){ 
    
    #rnd_stat <- lapply(1:permutations, function(x){
    
    ann_40 %>% 
      group_by(common_name,year) %>% 
      mutate(wind_speed = sample(wind_speed, replace = F)) %>% 
      group_by(common_name,year,stratum) %>% 
      arrange(desc(used), .by_group = TRUE) %>%
      summarize(max_minus_obs = max(wind_speed) - head(wind_speed,1)) %>% 
      mutate(perm = x)
    
  }) %>% 
    reduce(rbind) %>% 
    as.data.frame()
  
  Sys.time() - a # 6.626293 mins
  
  stopCluster(mycl)

save(rnd_stat, file = "R_files/rnd_stats_100_perm_df.RData")


#extract observed and random values for each stratum

#prep cluster
mycl <- makeCluster(detectCores() - 3)
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
  

Sys.time() - a # 56.95044 mins

stopCluster(mycl)


save(p_vals, file = "R_files/p_vals_100_perm_df.RData")




#Step 3: plot and conclusions #####
load("R_files/p_vals_100_perm_df.RData") #p_vals

p_vals$year_f <- as.character(p_vals$year)

#one plot per species, one curve per year
ggplot(p_vals) +
  stat_density(aes(x = p_more, group = year_f, color = year_f), adjust = 1.3,
               geom = "line", position = "identity", size = 0.55, alpha = 0.7) +
  xlab("P-value") + 
  ylab("Density") +
  theme_ipsum() +
  facet_wrap(~ common_name, ncol = 3) +
  guides(color = guide_legend(title="Year"))

#all species in one plot
ggplot(p_vals) +
  stat_density(aes(x = p_more, group = common_name, color = common_name), adjust = 1.3,
               geom = "line", position = "identity", size = 0.55, alpha = 0.8) +
  xlab("P-value") + 
  ylab("Density") +
  theme_ipsum() +
  guides(color = guide_legend(title="Year"))


#extract strata with p-value less than 0.05
sig <- p_vals %>% 
  filter(p_more <= 0.05)

sig_data <- ann_40 %>% 
  filter(stratum %in% sig$stratum)

#plot raw winds
par(mfrow = c(3,3))
for(i in unique(sig_data$stratum)){
  
  data <- sig_data[sig_data$stratum == i,]
  
  plot(density(data$wind_speed), main = data$common_name[1], ylab = "", xlab = "wind speed (m/s)")
  abline(v = data[data$used == 1, "wind_speed"], col = "red")
  
  if(i == "2009 2009-SULA12 2009-SULA12_3_1_22" ){
    legend("topright",legend = "used", col = "red", lty = 1,bty = "n")
  }
}

#plot permutation results
p_sig <- p_vals %>% 
  filter(stratum %in% sig$stratum)

#open random stats
load("R_files/rnd_stats_100_perm_df.RData") #rnd_stat

par(mfrow = c(3,3))
for(i in unique(p_sig$stratum)){
  
  data <- p_sig[p_sig$stratum == i,]
  rnds <- rnd_stat[rnd_stat$stratum == i,] 
    
  plot(density(rnds$max_minus_obs), main = data$common_name[1], ylab = "", xlab = "max_minus_obs_wind")
  abline(v = data$max_minus_obs, col = "red")
  
  if(i == "73500_1_7_15" ){
    legend("topleft",legend = "observed", col = "red", lty = 1,bty = "n")
  }
}



## ggplot?
ggplot(sig_data) +
  stat_density(aes(x = wind_speed), adjust = 1.3, geom = "line") +
  xlab("Wind speed (m/s)") + 
  ylab("Density") +
  theme_ipsum() +
  facet_wrap(~ stratum, ncol = 3) +
  geom_vline(xintercept = sig_data$wind_speed[sig_data$used == 1, "wind_speed"], linetype = "dotted", 
             color = "blue", size = 1.5) 

