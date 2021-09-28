#script for null modeling of seabird data
#Sep. 7. 2021. Elham Nourani, PhD. Konstanz, DE.
#update: Jun 7. 2021: make sure to separate the colonies for species from multiple colonies.
#https://www.jwilber.me/permutationtest/#:~:text=To%20calculate%20the%20p%2Dvalue,of%20test%2Dstatistics%20we%20calculated
#using z-scores as the test statistic instead of the difference between max wind and used wind.
#it is possible to calculate z scores for non-normally distributed data. 
#https://www.daylight.com/meetings/emug97/Bradshaw/Significant_Similarity/Z-scores.html
#https://www.statisticshowto.com/probability-and-statistics/z-score/
#Sep 23: change from calculating z scores to wilcoxon's signed rank test. skip null modeling.
#http://www.sthda.com/english/wiki/one-sample-wilcoxon-signed-rank-test-in-r#:~:text=The%20one%2Dsample%20Wilcoxon%20signed,value%20(i.e.%20theoretical%20value).
#additional step: remove strata with very low variation


library(tidyverse)
#library(cowplot) #raincloud plot
library(parallel)
#detach(package:plyr)
library(hrbrthemes)
#library(sm) #for one density plot per factor level

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")

#open dataset with alternative steps for hourly steps. prepped in random_steps.R
load("R_files/ssf_input_annotated_60_15_30alt_18spp.RData") #ann_30
  
ann_30 <- ann_30 %>%  
  mutate(group = paste(common_name, colony.name, sep = "_")) %>% 
  as.data.frame()


# Step 1: calc within-stratum variances ####

#calculate and plot within stratum variances
data_var <- ann_30 %>%
  group_by(stratum) %>%
  summarise(wspd_var = var(wind_speed),
            u_var = var(u10m),
            v_var = var(v10m),
            wspt_var = var(wind_support),
            species = head(common_name, 1),
            year = head(year, 1),
            wspd_cov = (sd(wind_speed)/mean(wind_speed))*100) #coef of variation


#hist(data_var$wspd_cov) #most values are between 0 and 10%. Let's look at those over 20%?

over20 <- data_var %>% 
  filter(wspd_cov > 20) #justification for 20%???!!



# Step 2: in each stratum, perform a wilcoxon signed-rank test (one sample) #####
# H0: sample median is greater than or equal to the theoretical value

wilcox_p <- ann_30 %>% 
  filter(stratum %in% over20$stratum) %>% 
  group_by(group, year, stratum) %>% #group by species-colony instead of just species
  arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached  detach("package:plyr", unload=TRUE)
  mutate(p = wilcox.test(tail(wind_speed,29), mu = head(wind_speed,1), alternative = "two.sided")$p.value, #try a two sided version
         observed_wspd = head(wind_speed,1), median_wspd = median(tail(wind_speed,29))) %>% #exclude observed wind from the group
  ungroup() %>% 
  as.data.frame()


#plot

wilcox_p$year_f <- as.character(wilcox_p$year)

#one plot per species, one curve per year
ggplot(wilcox_p) +
  stat_density(aes(x = p, group = year_f, color = year_f), adjust = 1.3,
               geom = "line", position = "identity", size = 0.55, alpha = 0.7) +
  xlab("P-value") + 
  ylab("Density") +
  theme_ipsum() +
  facet_wrap(~ group, ncol = 3) +
  guides(color = guide_legend(title="Year"))


#all species in one plot
ggplot(wilcox_p) +
  stat_density(aes(x = p, group = group, color = group), adjust = 1.3,
               geom = "line", position = "identity", size = 0.55, alpha = 0.8) +
  xlab("P-value") + 
  ylab("Density") +
  theme_ipsum() +
  guides(color = guide_legend(title="Year"))



##----------------------------------------------------

# Step 2: In each stratum, perform a wilcoxon signed-rank test (one sample) #####
# H0: sample median is greater than or equal to the theoretical value




summaries <- ann_30 %>% 
  group_by(group, year, stratum) %>% #group by species-colony instead of just species
  arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached  detach("package:plyr", unload=TRUE)
  #summarize(mean_wspd = mean(wind_speed), sd_wspd = sd(wind_speed), max_wspd = max(wind_speed), min_wspd = min(wind_speed), 
  #          obs_wspd = head(wind_speed,1)) %>% 
  mutate(p = wilcox.test(tail(wind_speed,29), mu = head(wind_speed,1), alternative = "two.sided")$p.value) %>% #exclude observed wind from the group
  ungroup() %>% 
  as.data.frame()
  


summaries <- ann_30 %>% 
  group_by(group, year, stratum) %>% #group by species-colony instead of just species
  arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached  detach("package:plyr", unload=TRUE)
  #summarize(mean_wspd = mean(wind_speed), sd_wspd = sd(wind_speed), max_wspd = max(wind_speed), min_wspd = min(wind_speed), 
  #          obs_wspd = head(wind_speed,1)) %>% 
  summarize(p = wilcox.test(tail(wind_speed,29), mu = head(wind_speed,1), alternative = "greater")$p.value,
            observed_wspd = head(wind_speed,1), median_wspd = median(tail(wind_speed,29))) %>% #exclude observed wind from the group
  ungroup() %>% 
  as.data.frame()


Ps <- ann_30 %>% 
  group_by(group, year, stratum) %>% #group by species-colony instead of just species
  arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached  detach("package:plyr", unload=TRUE)
  summarize(p = wilcox.test(wind_speed, mu = head(wind_speed,1), alternative = "greater")$p.value,
            observed_wspd = head(wind_speed,1), median_wspd = median(wind_speed)) %>% #include observed wind from the group
  ungroup() %>% 
  as.data.frame()




save(observed_zs, file = "R_files/observed_zs.RData")

#randomize and calculate the same statistic
#shuffle all wind speed values within each year, then recalc the statistic within each stratum

permutations <- 100

#prep cluster
mycl <- makeCluster(detectCores() - 2, setup_strategy = "sequential")
clusterExport(mycl, c("permutations", "ann_30")) 

clusterEvalQ(mycl, {
  library(dplyr)
})


a <- Sys.time()

rnd_zs <- parLapply(cl = mycl, X = c(1:permutations), fun = function(x){ 
  
  #rnd_stat <- lapply(1:permutations, function(x){
  ann_30 %>% 
    group_by(group,year) %>% 
    mutate(wind_speed = sample(wind_speed, replace = F)) %>% #shuffle the wind speed values
    group_by(group,year,stratum) %>% 
    arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached  detach("package:plyr", unload=TRUE)
    summarize(mean_wspd = mean(wind_speed), sd_wspd = sd(wind_speed), max_wspd = max(wind_speed), min_wspd = min(wind_speed), 
              obs_wspd = head(wind_speed,1)) %>% 
    mutate(z_rnd = (obs_wspd - mean_wspd)/sd_wspd) %>% 
    ungroup() %>% 
    as.data.frame()
  
}) %>% 
  reduce(rbind) %>% 
  as.data.frame()

Sys.time() - a # 5.7 min (100 perms)

stopCluster(mycl)

save(rnd_zs, file = "R_files/rnd_zs_100_perm_df.RData")

load("R_files/observed_zs.RData")
load("R_files/rnd_zs_100_perm_df.RData")

#extract observed and random values for each stratum

#prep cluster
mycl <- makeCluster(detectCores() - 3, setup_strategy = "sequential")
clusterExport(mycl, c("permutations", "observed_zs", "rnd_zs")) 

clusterEvalQ(mycl, {
  library(dplyr)
})


(a <- Sys.time())

p_vals <- parLapply(mycl, unique(observed_stat$stratum), function(x){
  obs <- observed_zs[observed_zs$stratum == x,]
  rnd <- rnd_zs[rnd_zs$stratum == x,]
  
  obs$p <- sum(rnd$z_rnd <= obs$z_obs)/permutations
  
  obs
  
}) %>% 
  reduce(rbind)
  

Sys.time() - a #1.340673

stopCluster(mycl)


save(p_vals, file = "R_files/p_vals_z_100_perm_df.RData")

#Step 2.1: plot and conclusions #####
load("R_files/p_vals_z_100_perm_df.RData") #p_vals

summaries$year_f <- as.character(summaries$year)
Ps$year_f <- as.character(Ps$year)

#one plot per species, one curve per year
ggplot(Ps) +
  stat_density(aes(x = p, group = year_f, color = year_f), adjust = 1.3,
               geom = "line", position = "identity", size = 0.55, alpha = 0.7) +
  xlab("P-value") + 
  ylab("Density") +
  theme_ipsum() +
  facet_wrap(~ group, ncol = 3) +
  guides(color = guide_legend(title="Year"))

#all species in one plot
ggplot(p_vals) +
  stat_density(aes(x = p, group = group, color = group), adjust = 1.3,
               geom = "line", position = "identity", size = 0.55, alpha = 0.8) +
  xlab("P-value") + 
  ylab("Density") +
  theme_ipsum() +
  guides(color = guide_legend(title="Year"))


#extract strata with p-value less than 0.05
sig <- p_vals %>% 
  filter(p_more <= 0.05)

sig_data <- ann_30 %>% 
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
load("R_files/rnd_stats_1000_perm_df.RData") #rnd_stat

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

