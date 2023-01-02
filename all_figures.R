#figures for seabirds and storms manuscript.
#follows from seabirds_public.R

#Aug 22, 2022
#Elham Nourani, PHD

library(tidyverse)
library(ggridges)
library(ggnewscale)
library(oce)
library(ggimage)

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")
source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")



# Fig 2: encountered wind vs wind at breeding range -------------------------------------

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/wind_raw.rdata")
lm_input <- readRDS("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/species_summary_data_with_ranges.RDS")

new_order <- lm_input %>%
  group_by(species) %>% 
  slice(1) %>% 
  arrange(desc(wing.loading..Nm.2.)) %>% 
  pull(species)

#assign sun icon to tropical species
lm_input <- lm_input %>% 
   mutate(image = ifelse(between(colony.lat, -23.43, 23.43), 
                         "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/sun.png", NA),
          species_f = factor(species, levels = new_order))
   

wind_raw <- wind_raw %>% 
  mutate(species_f = factor(species, levels = new_order),
         wind_data_f = factor(wind_data, levels = c("range", "gps_pts")),
         wind_speed_ms = round(wind_speed_ms, 2))



#use a sample of the data to nail the code
sample <- wind_raw %>% 
  group_by(species, wind_data) %>%
  slice(1:200) %>% 
  mutate(wind_data_f = factor(wind_data, levels = c("range", "gps_pts")))

clr <- oce::oceColorsPalette(120)[14]

#https://stackoverflow.com/questions/59152556/ggridges-color-gradient-per-group?rq=1


X11(width = 6.5, height = 10)
ggplot(sample, aes(x = wind_speed_ms, y = species_f)) + 
   stat_density_ridges(data = sample[sample$wind_data == "range",], color = "#A9A9A9", fill = "#A9A9A9",
                       jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1, point_alpha = 0.8, size = 0.2,
                       calc_ecdf = F, panel_scaling = F, alpha = 0.5,
                       scale = 1.5) +
   new_scale_fill() +
   stat_density_ridges(data = sample[sample$wind_data == "gps_pts",], aes( fill = stat(x), point_color = stat(x)),
                       jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2,
                       geom = "density_ridges_gradient", calc_ecdf = F, panel_scaling = F, 
                      scale = 1.5) +
   scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.5), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   scale_color_gradientn(aesthetics = "point_color",  colours = alpha(oce::oceColorsPalette(120)), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   new_scale_color() +
   scale_x_continuous(limits = c(-1, 28)) +
   geom_point(data = sample %>% group_by(species_f) %>% slice(1), 
              aes(x = -0.9, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
   scale_shape_manual(values = c("Dynamic soaring" = 4,"Flapping" = 0, "Thermal soaring" = 2, "Wind soaring" = 1)) +
   geom_image(data = lm_input, aes( x = 22, y = as.numeric(species_f) + 0.5, image = image),asp = 0.5, size = 0.05) +
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7))


#plot the complete version to disk
png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/distr_plot_sun.png", width = 6.5, height = 10, units = "in", res = 300)

ggplot(wind_raw, aes(x = wind_speed_ms, y = species_f)) + 
   stat_density_ridges(data = wind_raw[wind_raw$wind_data == "range",], color = "#A9A9A9", fill = "#A9A9A9",
                       jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1, point_alpha = 0.8, size = 0.2,
                       calc_ecdf = F, panel_scaling = F, alpha = 0.5,
                       scale = 1.5) +
   new_scale_fill() +
   stat_density_ridges(data = wind_raw[wind_raw$wind_data == "gps_pts",], aes( fill = stat(x), point_color = stat(x)),
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
   geom_point(data = wind_raw %>% group_by(species_f) %>% slice(1),
              aes(x = -0.9, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
   scale_shape_manual(values = c("Dynamic soaring" = 4,"Flapping" = 0, "Thermal soaring" = 2, "Wind soaring" = 1)) +
   geom_image(data = lm_input, aes( x = 22, y = as.numeric(species_f) + 0.5, image = image),asp = 0.5, size = 0.05) +
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7))


dev.off()

# Fig 3: Wind stills -------------------------------------

#load files prepared in wind_fields.R

load( "R_files/trips_df_for_windfields.RData") #trips_df
load("R_files/avoidance_hour_per_trip.RData") #avoidance_hour
wind_ls <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/processed/", pattern = "wind", full.names = T)
gps_ls <- setdiff(list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/wind_fields/processed/", full.names = T),wind_ls)


 
# Fig 4: linear models -------------------------------------

lm_input <- readRDS("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/species_summary_data_with_ranges.RDS")
 

long_df <- lm_input %>% 
   group_by(species) %>% 
   arrange(desc(max_wind_ms)) %>% 
   slice(1) %>% 
   pivot_longer(cols = c("max_wind_ms", "range_median", "range_max"),
                names_to = "wind_source",
                values_to = "wind_speed")
 
clr <- oce::oceColorsPalette(120)[14]
clr2 <- oce::oceColorsPalette(120)[105]


png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_lms.png", width = 7, height = 5.5, units = "in", res = 300)

ggplot(long_df, aes(x = wing.loading..Nm.2., y = wind_speed, col = wind_source, fill = wind_source)) +
  geom_smooth(aes(y = wind_speed, group = wind_source), method = "lm", alpha = .1, level = .95) + #95% standard error
  geom_point(aes(y = wind_speed, shape = flight_style), size = 1.5, stroke = 0.8) +
  labs(x = expression("Wing loading (Nm"^-2*")"),
       y = expression("Wind speed (m s"^-1*")")) +
  scale_shape_manual(values = c(4,0,2,1)) + #filled points: c(15, 17, 19)
  scale_color_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range median")) +
  scale_fill_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range median")) +
  theme_minimal() + #theme_ipsum() looks better
  #theme(axis.title.y = element_text(size = 13)) +
  guides(shape = guide_legend("Flight style:"))
 
dev.off()

#add stats 
ggplot(long_df, aes(x = wing.loading..Nm.2., y = wind_speed, col = wind_source, fill = wind_source)) +
   geom_smooth(aes(y = wind_speed, group = wind_source), method = "lm", alpha = .1, level = .95) + #95% standard error
   geom_point(aes(y = wind_speed, shape = flight_style), size = 2, stroke = 0.8) +
   labs(x = expression("Wing loading (Nm"^-2*")"),
        y = "Wind speed") +
   scale_shape_manual(values = c(4,0,2,1)) + #filled points: c(15, 17, 19)
   scale_color_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range mediam")) +
   scale_fill_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range mediam")) +
   theme_minimal() + #theme_ipsum() looks better
   theme(axis.title.y = element_text(size = 13)) +
   guides(shape = guide_legend("Flight style:"))

###alternative. have one panel per plot

#library(ggpp)
#library(gginnards)
library(ggpmisc)

ggplot(long_df, aes(x = wing.loading..Nm.2., y = wind_speed, col = wind_source, fill = wind_source)) +
   stat_poly_line() +
   stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                  after_stat(rr.label), sep = "*\", \"*"))) +
   #geom_smooth(aes(y = wind_speed, group = wind_source), method = "lm", alpha = .1, level = .95) + #95% standard error
   geom_point(aes(y = wind_speed, shape = flight_style), size = 2, stroke = 0.8) +
   labs(x = expression("Wing loading (Nm"^-2*")"),
        y = "Wind speed") +
   scale_shape_manual(values = c(4,0,2,1)) + #filled points: c(15, 17, 19)
   scale_color_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range mediam")) +
   scale_fill_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range mediam")) +
   theme_minimal() + #theme_ipsum() looks better
   theme(axis.title.y = element_text(size = 13)) +
   guides(shape = guide_legend("Flight style:"))  +
   facet_wrap(~wind_source)

 
 