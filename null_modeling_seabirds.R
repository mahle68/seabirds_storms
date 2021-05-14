#script for null modeling of seabird data
#May 11. 2021. Elham Nourani, PhD. Radolfzell, DE.
#https://www.jwilber.me/permutationtest/#:~:text=To%20calculate%20the%20p%2Dvalue,of%20test%2Dstatistics%20we%20calculated


library(tidyverse)
library(cowplot) #raincloud plot
library(parallel)


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
  arrange(desc(used), .by_group = TRUE) %>%
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


#plot the random and observed values
par(mfrow = c(3,2))
for(i in c("tradewind","temperate")){
  zone <- i
  for(j in c("delta_t", "u", "v"))
    plot(1:permutations, rnd_d_rsd[rnd_d_rsd$zone == i, grep(j,colnames(rnd_d_rsd))], type = "l", main = paste(i, "delta rsd in", j, sep = " "))
  abline(h = obs_d_rsd[obs_d_rsd$zone == i, grep(j,colnames(obs_d_rsd))], col = "red")
}

par(mfrow = c(3,2))
for(i in c("tradewind","temperate")){
  for(j in c("delta_t", "u", "v"))
    hist(rnd_d_rsd[rnd_d_rsd$zone == i, grep(j,colnames(rnd_d_rsd))], breaks = 50, col = "lightgrey", 
         xlim = c(0, obs_d_rsd[obs_d_rsd$zone == i, grep(j,colnames(obs_d_rsd))] + 0.5),main = paste(i, "delta rsd in", j, sep = " "))
  abline(v = obs_d_rsd[obs_d_rsd$zone == i, grep(j,colnames(obs_d_rsd))], col = "red")
}

par(mfrow = c(1,2))

##### STEP 5: calculate p-values 
p_values <- data.frame(NULL)

for(i in c("tradewind","temperate")){
  for(j in c("delta_t", "u", "v")){
    p <- sum(obs_d_rsd[obs_d_rsd$zone == i, grep(j,colnames(obs_d_rsd))] <= rnd_d_rsd[rnd_d_rsd$zone == i, grep(j,colnames(rnd_d_rsd))]) / permutations
    p_values[i, j] <- p
  }
}


