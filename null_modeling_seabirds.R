#script for null modeling of seabird data
#May 11. 2021. Elham Nourani, PhD. Radolfzell, DE.
#update: Jun 7. 2021: make sure to separate the colonies for species from multiple colonies.
#https://www.jwilber.me/permutationtest/#:~:text=To%20calculate%20the%20p%2Dvalue,of%20test%2Dstatistics%20we%20calculated


library(tidyverse)
#library(cowplot) #raincloud plot
library(parallel)
#detach(package:plyr)
library(hrbrthemes)
#library(sm) #for one density plot per factor level
library(ggridges)
#https://www.datanovia.com/en/blog/elegant-visualization-of-density-distribution-in-r-using-ridgeline/


setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")


#source("/home/mahle68/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/R_rainclouds.R")
#source("/home/mahle68/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/summarySE.R")
#source("/home/mahle68/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/simulateData.R")

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
            group = head(group,1),
            year = head(year, 1),
            wspd_cov = (sd(wind_speed)/mean(wind_speed))*100) %>% #coef of variation
  mutate(log_wspd_cov = log(wspd_cov)) %>% 
  ungroup() %>% 
  as.data.frame()

#plot
# logdensity <- function (x, bw = "SJ") 
# {
#   y <- log(x)
#   g <- density(y, bw = bw, n = 1001)
#   xgrid <- exp(g$x)
#   g$y <- c(0, g$y/xgrid)
#   g$x <- c(0, xgrid)
#   return(g)
# }
# 
# 
# 
# #dens <- apply(data_var, 2, density)
# 
# plot(x = "", xlim = c(0,10), ylim = c(0,1), ylab = "Density", xlab = "Log of coeffient of variation", bty = "l")
# #mapply(lines, group, col=1:length(group))
# 
# 
# for(i in levels(as.factor(data_var$group))){
#   
#   fit <- logdensity(data_var[data_var$group == i, "log_wspd_cov"][data_var[data_var$group == i, "log_wspd_cov"] > 0]) # Only take density of positive part
#   lines(fit$x,fit$y*mean(data_var[data_var$group == i, "log_wspd_cov"] > 0), col = "red") # Scale density by proportion positive
#   
# #lines(density(data_var[data_var$group == i, "wspd_cov"], from = 0), ylab = "", xlab = "", add = T)
# }
#   
# plot(density(data$wind_speed), main = data$common_name[1], ylab = "", xlab = "wind speed (m/s)")
# 
# 
# p2 <- ggplot(data=data_var, aes(x=wspd_cov, group=group, fill=group)) +
#   geom_density(adjust=1.5, alpha=.4) +
#   theme_ipsum()





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


#ggsave("Cov_bar.png", CoV_bar, width = 11, height = 12, units = "cm", path = "/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs", device = "png")

png("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/Cov_bar.png", width = 5, height = 5, units = "in", res = 300)
print(CoV_bar)
dev.off()


# Step 2: permutation test- In each stratum, is the difference between selected and max available wind speed higher/lower than expected by chance? ####

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



#Step 3: plot and conclusions #####
load("R_files/p_vals_1000_perm_df.RData") #p_vals

p_vals <- p_vals %>% 
  mutate(year_f = as.character(year),
         species = str_split(group, "\\_", simplify = TRUE)[, 1]) %>% 
  ungroup() %>% 
  as.data.frame()

p_vals[p_vals$species == "Galapagos albatross (waved albatross)", "species"] <- "Galapagos albatross"

#plot
# ggplot(p_vals, aes(x = p_more, y = group)) + 
#   geom_density_ridges(scale = 3, alpha = 0.5) + 
#   scale_x_continuous(limits = c(0, 1)) +
#   labs(y = "", x = "mmm") +
#   theme_minimal()

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


# #plot raw winds
# par(mfrow = c(3,3))
# for(i in unique(sig_data$stratum)){
# 
#   data <- sig_data[sig_data$stratum == i,]
# 
#   plot(density(data$wind_speed), main = data$common_name[1], ylab = "", xlab = "wind speed (m/s)")
#   abline(v = data[data$used == 1, "wind_speed"], col = "red")
# 
#   if(i == "2009 2009-SULA12 2009-SULA12_3_1_22" ){
#     legend("topright",legend = "used", col = "red", lty = 1,bty = "n")
#   }
# }



# sig_strata <- ggplot(sig_data, aes(x = wind_speed, y = stratum, height = stat(density))) + 
#   geom_density_ridges(
#     stat = "binline", bins = 20, scale = 0.98, alpha = 0.3,
#     draw_baseline = FALSE
#   ) +
#   scale_x_continuous(limits = c(-4, 25)) +
#   labs(y = "Density", x = "Coefficient of variation (%)") +
#   theme_minimal()


# 
# 
# #calculate proportion of sig p values for each species
# sig <- p_vals %>% 
#   group_by(group) %>% 
#   summarize(sum_non_sig = sum(p_more > 0.05),
#             n_strata = n()) %>% 
#   mutate(prop_nonsig = sum_non_sig / n_strata)
# 
# #one plot per species, one curve per year
# ggplot(p_vals) +
#   stat_density(aes(x = p_more, group = year_f, color = year_f), adjust = 1.3,
#                geom = "line", position = "identity", size = 0.55, alpha = 0.7) +
#   xlab("P-value") + 
#   ylab("Density") +
#   theme_ipsum() +
#   facet_wrap(~ group, ncol = 3) +
#   guides(color = guide_legend(title="Year"))
# 
# #all species in one plot
# ggplot(p_vals) +
#   stat_density(aes(x = p_more, group = group, color = group), adjust = 1.3,
#                geom = "line", position = "identity", size = 0.55, alpha = 0.8) +
#   xlab("P-value") + 
#   ylab("Density") +
#   theme_ipsum() +
#   guides(color = guide_legend(title="Year"))
# 
# 
# #extract strata with p-value less than 0.05
# sig <- p_vals %>% 
#   filter(p_more <= 0.05)
# 
# sig_data <- ann_30 %>% 
#   filter(stratum %in% sig$stratum)
# 
# #plot raw winds
# par(mfrow = c(3,3))
# for(i in unique(sig_data$stratum)){
#   
#   data <- sig_data[sig_data$stratum == i,]
#   
#   plot(density(data$wind_speed), main = data$common_name[1], ylab = "", xlab = "wind speed (m/s)")
#   abline(v = data[data$used == 1, "wind_speed"], col = "red")
#   
#   if(i == "2009 2009-SULA12 2009-SULA12_3_1_22" ){
#     legend("topright",legend = "used", col = "red", lty = 1,bty = "n")
#   }
# }
# 
# #plot permutation results
# p_sig <- p_vals %>% 
#   filter(stratum %in% sig$stratum)
# 
# #open random stats
# load("R_files/rnd_stats_1000_perm_df.RData") #rnd_stat
# 
# par(mfrow = c(3,3))
# for(i in unique(p_sig$stratum)){
#   
#   data <- p_sig[p_sig$stratum == i,]
#   rnds <- rnd_stat[rnd_stat$stratum == i,] 
#     
#   plot(density(rnds$max_minus_obs), main = data$common_name[1], ylab = "", xlab = "max_minus_obs_wind")
#   abline(v = data$max_minus_obs, col = "red")
#   
#   if(i == "73500_1_7_15" ){
#     legend("topleft",legend = "observed", col = "red", lty = 1,bty = "n")
#   }
# }
# 
# 
# 
# ## ggplot?
# ggplot(sig_data) +
#   stat_density(aes(x = wind_speed), adjust = 1.3, geom = "line") +
#   xlab("Wind speed (m/s)") + 
#   ylab("Density") +
#   theme_ipsum() +
#   facet_wrap(~ stratum, ncol = 3) +
#   geom_vline(xintercept = sig_data$wind_speed[sig_data$used == 1, "wind_speed"], linetype = "dotted", 
#              color = "blue", size = 1.5) 
# 
