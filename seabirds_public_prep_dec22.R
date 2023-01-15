#scripts for reproducing the results of Nourani et al. 2023, Current Biology
#session info is provided at the end of the script
#Elham Nourani, PhD. Dec. 2022. enourani@ab.mpg.de

#--------------------------------------------------------------------------

library(tidyverse)
library(ggridges)
library(lubridate)
library(parallel)
library(oce)
library(sp)
library(ggnewscale)
library(ggimage)

#----- STEP 1: input data --------------------------------------------
#data is available via the Edmond repository: "seabirds' wind niche", https://doi.org/10.17617/3.8U7EHD
ann_30 <- readRDS("used_alt_annotated.rds") #this file contains used and alternative steps (hourly). This will be used for permutation tests. Alternative steps were generated following the workflow in: https://github.com/mahle68/global_seascape_public/blob/main/step_generation.R
p_vals <- readRDS("p_vals_1000_perm.rds") #this file contains p_values resulting from permutation test done in step 3 below.
lm_input <- readRDS("species_summary_data_with_ranges.RDS") #this file contains summary information for each species including morphology, wind conditions experienced, and wind conditions at breeding range
ann <- readRDS("annotated_data_hrly.rds") #this file contains sub-sampled hourly data annotated with atmospheric information (this is the full dataset that was used for generating alternative steps)
raw_wind <- readRDS("wind_at_breeding_range.rds") #this file contains wind speed data at the breeding range (during breeding season) for each species. Caution: large file!

#----- STEP 2: calculate within-stratum variances (+ plot Fig S1)-------------------------

#calculate and plot within stratum variances
data_var <- ann_30 %>%
  group_by(stratum) %>% #group by stratum to estimate variation in wind speeds within each stratum separately
  summarise(wspd_cov = (sd(wind_speed)/mean(wind_speed))*100, #coefficient of variation of wind speed
            species = head(species, 1)) %>%
  ungroup()

# ------------------------------------------------------ plot Figure S1:
X11(width = 5, height = 5)
ggplot(data_var, aes(x = wspd_cov, y = reorder(as.factor(species), desc(as.factor(species))), height = stat(density))) + 
  geom_density_ridges(
    stat = "binline", bins = 20, scale = 0.98, alpha = 0.3,
    draw_baseline = FALSE
  ) +
  scale_x_continuous(limits = c(-4, 140)) +
  labs(y = "", x = "Coefficient of variation (%)") +
  theme_minimal()

#----- STEP 3: permutation test (+ plot Figs S2, S3)--------------------------------------

#In each stratum, is the difference between selected and max available wind speed higher than expected by chance?
#for each stratum, calculate the difference between observed wind speed and max wind speed (incl. observed)

observed_stat <- ann_30 %>%
  group_by(group, year, stratum) %>% #group by species-colony (group) instead of just species
  arrange(desc(used), .by_group = TRUE) %>% #make sure plyr is detached detach("package:plyr", unload=TRUE)
  summarize(max_minus_obs = max(wind_speed) - head(wind_speed,1))

#randomize and calculate the same statistic
#shuffle all wind speed values within each year, then recalculate the statistic within each stratum

permutations <- 1000

#prepare cluster for multi-core analysis.
mycl <- makeCluster(detectCores() - 2, setup_strategy = "sequential") #choose the number of cores based on your machine
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


# P_vals is available as "p_vals_1000_perm.rds"

# ------------------------------------------------------ plot Figure S2:

X11(width = 5, height = 5)
ggplot(p_vals, aes(x = p_more, y = reorder(as.factor(species), desc(as.factor(species))), height = stat(density))) + 
  geom_density_ridges(
    stat = "binline", bins = 20, scale = 0.98, alpha = 0.3,
    draw_baseline = FALSE
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(y = "", x = "Significance") +
  theme_minimal()

# ------------------------------------------------------ plot Figure S3:

#extract strata with significant avoidance of strong winds (p-value less than 0.05)
sig <- p_vals %>%
  filter(p_more <= 0.05)

#extract data for the significant strata
sig_data <- ann_30 %>%
  filter(stratum %in% sig$stratum)

sig_data$stratum <- reorder(as.factor(sig_data$stratum),desc(sig_data$species))

#extract the used wind speeds
used_wind <- sig_data %>% 
  group_by(stratum) %>% 
  filter(used == 1) %>% 
  ungroup() %>% 
  dplyr::select(c("stratum", "species", "wind_speed")) %>% 
  as.data.frame()

X11(width = 7, height = 4)
ggplot(sig_data, aes(x = wind_speed, y = stratum)) + 
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


#----- STEP 4: linear models (+ plot Fig 4)--------------------------------------

# maximum encountered wind as a function of wing loading
m1 <- lm(max_wind ~ wing.loading..Nm.2., data = lm_input) # adjRsq = 0.31  

# wing loading as a function of median wind conditions at breeding range
m2 <- lm(wing.loading..Nm.2. ~ range_median, data = lm_input) # adjRsq = 0.3463   

# wing loading as a function of maximum wind conditions at breeding range
m3 <- lm(wing.loading..Nm.2. ~ range_max, data = lm_input) # adjRsq = -0.04219 

# max encountered wind as a function of median wind speed at the breeding range
m4 <- lm(max_wind ~ range_median, data = lm_input) # adjRsq = 0.31

# max encountered wind as a funciton of max wind speed at the breeding range
m5 <- lm(max_wind ~ range_max, data = lm_input) # adjRsq = 0.31

# ------------------------------------------------------ plot Figure 4:
#convert to long form for plotting
long_df <- lm_input %>% 
  group_by(species) %>% 
  arrange(desc(max_wind_ms)) %>% 
  slice(1) %>% 
  pivot_longer(cols = c("max_wind_ms", "range_median", "range_max"),
               names_to = "wind_source",
               values_to = "wind_speed")

#pick colors from the oce package
clr <- oce::oceColorsPalette(120)[14]
clr2 <- oce::oceColorsPalette(120)[105]

ggplot(long_df, aes(x = wing.loading..Nm.2., y = wind_speed, col = wind_source, fill = wind_source)) +
  geom_smooth(aes(y = wind_speed, group = wind_source), method = "lm", alpha = .1, level = .95) + #95% standard error
  geom_point(aes(y = wind_speed, shape = flight_style), size = 1.5, stroke = 0.8) +
  labs(x = expression("Wing loading (Nm"^-2*")"),
       y = expression("Wind speed (m s"^-1*")")) +
  scale_shape_manual(values = c(4,0,2,1)) +
  scale_color_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range median")) +
  scale_fill_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range median")) +
  theme_minimal() + 
  guides(shape = guide_legend("Flight style:"))


  
#greyscale

png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_lms_bw.png", width = 7, height = 5.5, units = "in", res = 300)

ggplot(long_df, aes(x = wing.loading..Nm.2., y = wind_speed, linetype = wind_source)) +
  geom_smooth(aes(y = wind_speed, group = wind_source), color = "black", method = "lm", alpha = .1, level = .95) + #95% standard error
  geom_point(aes(y = wind_speed, shape = flight_style), size = 1.5, stroke = 0.8) +
  labs(x = expression("Wing loading (Nm"^-2*")"),
       y = expression("Wind speed (m s"^-1*")")) +
  scale_shape_manual(values = c(4,0,2,1)) +
  scale_linetype_manual(values = c("dashed", "twodash", "dotted"), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range median")) +
  scale_fill_manual(values = c("dashed", "twodash", "dotted"), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range median")) +
  theme_minimal() + 
  guides(shape = guide_legend("Flight style:"))

dev.off()
  
#three panels instead of one. the flight style points get mixed up in one panel

X11(width = 10, height = 4)

png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/Curr_Biol/revision_2/files_to_submit/linear_models.png", width = 10, height = 4)

ggplot(long_df, aes(x = wing.loading..Nm.2., y = wind_speed)) +
  geom_smooth(aes(y = wind_speed, group = wind_source), color = "black", method = "lm", alpha = .1, level = .95, lwd = 0.5) + #95% standard error
  geom_point(aes(y = wind_speed, shape = flight_style), size = 0.9, stroke = 0.8) +
  labs(x = expression("Wing loading (Nm"^-2*")"),
       y = expression("Wind speed (m s"^-1*")")) +
  scale_shape_manual(values = c(4,0,2,1)) +
  facet_wrap(~ wind_source, labeller = labeller(wind_source = c( "max_wind_ms" = "Maximum wind speed \n encountered", 
                                                                 "range_max" = "Maximum wind speed \n at breeding range", 
                                                                 "range_median" = "Median wind speed \n at breeding range"))) +
  theme_minimal() + 
  guides(shape = guide_legend("Flight style:")) +
  theme(legend.title = element_text(size = 11),
        strip.text = element_text(size = 12,face = "bold"))

dev.off()

#----- STEP 5: test for spatio-temporal autocorrelation --------------------------------------

#temporal correlation. result: no temporal autocorrelation
acf(resid(m1))
acf(resid(m2))
acf(resid(m3))
acf(resid(m4))
acf(resid(m5))

#spatial autocorrelation. bubble plot result: no autocorrelation
spdata <- data.frame(resid = resid(m1), x = lm_input$colony_long, y = lm_input$colony_lat)
spdata <- data.frame(resid = resid(m2), x = lm_input$colony_long, y = lm_input$colony_lat)
spdata <- data.frame(resid = resid(m3), x = lm_input$colony_long, y = lm_input$colony_lat)
spdata <- data.frame(resid = resid(m4), x = lm_input$colony_long, y = lm_input$colony_lat)
spdata <- data.frame(resid = resid(m5), x = lm_input$colony_long, y = lm_input$colony_lat)


#convert to a spatial object
coordinates(spdata) <-~ x + y
bubble(spdata, "resid", col = c("blue","orange"))

#----- STEP 6: test for correlation between data quantity and max wind speeds --------------------------------------

summary_info <- ann %>% 
  group_by(sci_name, colony_name) %>% 
  summarize(n_trips = n_distinct(trip_id)) %>% 
  full_join(lm_input, by = c("sci_name", "colony_name")) %>% 
  mutate(max_wind_ms = max_wind/3.6) %>% 
  as.data.frame()

cor.test(summary_info$n_trips, summary_info$max_wind_ms) #R = 0.149

#----- STEP 7: Plot Figure 1 ----------------------------------------------------------------

#plotting this figure is computationally demanding. I suggest writing it directly to file.

#order the species names based on wing loading
new_order <- lm_input %>%
  group_by(species) %>% 
  slice(1) %>% 
  arrange(desc(wing.loading..Nm.2.)) %>% 
  pull(species)

#assign sun icon to tropical species. sun png is included in the edmond repo
lm_input <- lm_input %>% 
  mutate(image = ifelse(between(colony_lat, -23.43, 23.43), 
                        "sun.png", NA), 
         species_f = factor(species, levels = new_order))


raw_wind <- raw_wind %>% 
  mutate(species_f = factor(species, levels = new_order),
         wind_data_f = factor(wind_data, levels = c("range", "gps_pts")),
         wind_speed_ms = round(wind_speed_ms, 2))

#extract a color from the oce palette
clr <- oce::oceColorsPalette(120)[14]

#write the plot to disk
pdf("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/Curr_Biol/revision_2/files_to_submit/wind_distribution.pdf", width = 6.5, height = 10)

ggplot(raw_wind, aes(x = wind_speed_ms, y = species_f)) + 
  stat_density_ridges(data = raw_wind[raw_wind$wind_data == "range",], color = "#A9A9A9", fill = "#A9A9A9",
                      jittered_points = TRUE, rel_min_height = .01,
                      point_shape = "|", point_size = 1, point_alpha = 0.8, size = 0.2,
                      calc_ecdf = F, panel_scaling = F, alpha = 0.5,
                      scale = 1.5) +
  new_scale_fill() +
  stat_density_ridges(data = raw_wind[raw_wind$wind_data == "gps_pts",], aes( fill = stat(x), point_color = stat(x)),
                      jittered_points = TRUE, rel_min_height = .01,
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2,
                      geom = "density_ridges_gradient", calc_ecdf = F, panel_scaling = F, 
                      scale = 1.5) +
  scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.6), limits = c(0,23), 
                       na.value = "white", guide = 'none') +
  scale_color_gradientn(aesthetics = "point_color",  colours = alpha(oce::oceColorsPalette(120)), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
  new_scale_color() +
  scale_x_continuous(limits = c(-1, 28)) +
  geom_point(data = raw_wind %>% group_by(species_f) %>% slice(1),
             aes(x = -0.9, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = "black") +
  scale_shape_manual(values = c("Dynamic soaring" = 4,"Flapping" = 0, "Thermal soaring" = 2, "Wind soaring" = 1)) +
  geom_image(data = lm_input, aes( x = 22, y = as.numeric(species_f) + 0.5, image = image),asp = 0.5, size = 0.05) +
  labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
  theme_minimal() +
  guides(shape = guide_legend("Flight style:")) +
  theme(legend.position = "bottom",legend.title = element_text(size = 10), 
        legend.text=element_text(size = 7))

dev.off()

#----- STEP 8: Plot Figure 3 ----------------------------------------------------------------
#based on https://semba-blog.netlify.app/10/29/2018/animating-oceanographic-data-in-r-with-ggplot2-and-gganimate/

#SESSION INFO ----------------------------------------------------------------
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Pop!_OS 22.04 LTS
# 
# Matrix products: default
# BLAS:   /usr/local/lib/R/lib/libRblas.so
# LAPACK: /usr/local/lib/R/lib/libRlapack.so
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggimage_0.3.1    ggnewscale_0.4.7 sp_1.4-7         oce_1.7-2        gsw_1.0-6        lubridate_1.8.0  ggridges_0.5.3   forcats_0.5.1    stringr_1.4.1    dplyr_1.0.9      purrr_0.3.4     
# [12] readr_2.1.2      tidyr_1.2.0      tibble_3.1.8     ggplot2_3.3.6    tidyverse_1.3.1 
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.8.3       lattice_0.20-45    png_0.1-7          assertthat_0.2.1   digest_0.6.29      utf8_1.2.2         R6_2.5.1           cellranger_1.1.0   plyr_1.8.7         backports_1.4.1   
# [11] reprex_2.0.1       ggfun_0.0.7        httr_1.4.2         pillar_1.8.1       yulab.utils_0.0.5  rlang_1.0.4        readxl_1.4.0       rstudioapi_0.13    magick_2.7.3       Matrix_1.4-1      
# [21] labeling_0.4.2     splines_4.2.0      rgdal_1.5-31       munsell_0.5.0      broom_1.0.0        compiler_4.2.0     modelr_0.1.8       gridGraphics_0.5-1 pkgconfig_2.0.3    mgcv_1.8-40       
# [31] tidyselect_1.1.2   fansi_1.0.3        crayon_1.5.1       tzdb_0.3.0         dbplyr_2.1.1       withr_2.5.0        grid_4.2.0         nlme_3.1-157       jsonlite_1.8.0     gtable_0.3.0      
# [41] lifecycle_1.0.1    DBI_1.1.2          magrittr_2.0.3     scales_1.2.1       cli_3.3.0          stringi_1.7.8      farver_2.1.1       fs_1.5.2           xml2_1.3.3         ellipsis_0.3.2    
# [51] generics_0.1.3     vctrs_0.4.1        tools_4.2.0        ggplotify_0.1.0    glue_1.6.2         hms_1.1.1          abind_1.4-5        colorspace_2.0-3   rvest_1.0.2        haven_2.5.0 