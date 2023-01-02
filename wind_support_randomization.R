#trying the seabird permutation on wind support
#OCT 20.2022
#Elham Nourani, PhD. enourani@ab.mpg.de

#code copied from seabirds_public.R and wind speed replaced with wind support


library(tidyverse)
library(lubridate)
library(sf)
library(parallel)
library(move)
library(ggridges)

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/used_alt_annotated.RData") #ann_30


observed_stat <- ann_30 %>% 
  group_by(group, year, stratum) %>% #group by species-colony (group) instead of just species
  arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached detach("package:plyr", unload=TRUE)
  summarize(max_minus_obs = max(wind_support) - head(wind_support,1))

#randomize and calculate the same statistic
#shuffle all wind speed values within each year, then recalculate the statistic within each stratum

permutations <- 50

#prep cluster
mycl <- makeCluster(detectCores() - 7, setup_strategy = "sequential")
clusterExport(mycl, c("permutations", "ann_30")) 

clusterEvalQ(mycl, {
  library(dplyr)
})

rnd_stat <- parLapply(cl = mycl, X = c(1:permutations), fun = function(x){ #if not using multiple cores, call lapply instead of ParLapply
  
  ann_30 %>% 
    group_by(group,year) %>% 
    mutate(wind_support = sample(wind_support, replace = F)) %>% 
    group_by(group,year,stratum) %>% 
    arrange(desc(used), .by_group = TRUE) %>%
    summarize(max_minus_obs = max(wind_support) - head(wind_support,1)) %>% 
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

X11(width = 5, height = 5)

perm_sig <- ggplot(p_vals, aes(x = p_more, y = reorder(as.factor(group), desc(as.factor(group))), height = stat(density))) + 
  geom_density_ridges(
    stat = "binline", bins = 20, scale = 0.98, alpha = 0.3,
    draw_baseline = FALSE
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "", x = "Significance") +
  theme_minimal()

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
