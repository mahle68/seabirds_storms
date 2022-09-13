#figures for seabirds and storms manuscript.
#follows from seabirds_public.R

#Aug 22, 2022
#Elham Nourani, PHD

library(tidyverse)
library(ggridges)
library(ggnewscale)
library(oce)




# Fig 2: encountered wind vs wind at breeding range -------------------------------------

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/wind_raw.rdata")
lm_input <- readRDS("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/data_public/species_summary_data_with_ranges.RDS")

new_order <- lm_input %>%
  group_by(species) %>% 
  slice(1) %>% 
  arrange(desc(wing.loading..Nm.2.)) %>% 
  pull(species) 

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


#ridgeline plots with the original color scheme

ggplot(sample[sample$wind_data == "gps_pts",], aes(x = wind_speed_ms, y = species_f, fill = stat(x), point_color = stat(x))) + 
   geom_density_ridges_gradient(rel_min_height = .01, jittered_points = TRUE,point_shape = "|", point_size = 1.2) +
   scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   scale_colour_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                          na.value = "white", guide = 'none') +
   theme_minimal()


ggplot(sample, aes(x = wind_speed_ms,y = species_f, fill = stat(x), point_color = stat(x), group = wind_data)) + 
   geom_density_ridges_gradient(rel_min_height = .01, jittered_points = TRUE,point_shape = "|", point_size = 1.2) +
   scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   scale_colour_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                          na.value = "white", guide = 'none') +
   theme_minimal()


ggplot(sample[sample$wind_data == "gps_pts",], aes(x = wind_speed_ms, y = species_f, fill = stat(x), point_color = stat(x))) + 
   geom_point(data = sample[sample$wind_data == "range",], aes(x = wind_speed_ms, y = species_f), shape = "|", show.legend = F,
              color = "#A9A9A9", size = 1.2, alpha = 0.8) +
   geom_density_ridges_gradient() +
   #stat_density_ridges(jittered_points = TRUE, rel_min_height = .01,
   #                    point_shape = "|", point_size = 1.2, point_alpha = 0.8, size = 0.25,
   #                    geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, scale = 3, #only the median line
   #                    quantiles = 0.5, quantile_lines = T) +
   scale_x_continuous(limits = c(-1, 30)) +
   scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   scale_colour_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                         na.value = "white", guide = 'none') + 
   geom_point(data = sample, aes(x = -0.8, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
   scale_shape_manual(values = c("Dynamic soaring" = 4,"Flapping" = 0, "Thermal soaring" = 2, "Wind soaring" = 1)) +
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7))


ggplot(sample, aes(x = wind_speed_ms, y = species_f)) + 

   stat_density_ridges(data = sample[sample$wind_data == "range",], color = "#A9A9A9", fill = "#e6e3e3",
                       jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1.3, point_alpha = 0.8, size = 0.2,
                       calc_ecdf = TRUE, panel_scaling = F, alpha = 0.5,
                       scale = 1.5) +
   new_scale_fill() +
   stat_density_ridges(data = sample[sample$wind_data == "gps_pts",], aes( fill = stat(x), point_color = stat(x)),
                       jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1.3, size = 0.2,
                       geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, #only the median line
                      scale = 1.5) +
   scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   scale_color_gradientn( aesthetics = "point_color",  colours = alpha(oce::oceColorsPalette(120), alpha = 1), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   new_scale_color() +
   scale_x_continuous(limits = c(-1, 28)) +
   geom_point(data = sample, aes(x = -0.9, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
   scale_shape_manual(values = c("Dynamic soaring" = 4,"Flapping" = 0, "Thermal soaring" = 2, "Wind soaring" = 1)) +
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7))


#plot the complete version to disk
png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/distr_plot_jitter_only.png", width = 6.5, height = 10, units = "in", res = 300)

ggplot(wind_raw[wind_raw$wind_data == "gps_pts",], aes(x = wind_speed_ms, y = species_f, fill = stat(x))) + 
   geom_point(data = wind_raw[wind_raw$wind_data == "range",], aes(x = wind_speed_ms, y = species_f), shape = "|", show.legend = F,
              color = "#A9A9A9", size = 1.3, alpha = 0.8) +
   stat_density_ridges(jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1.3, point_alpha = 0.8, size = 0.25,
                       geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, #only the median line
                       quantiles = 0.5, quantile_lines = T, scale = 3) +
   scale_x_continuous(limits = c(-1, 28)) +
   scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   geom_point(data = wind_raw, aes(x = -0.9, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
   scale_shape_manual(values = c("Dynamic soaring" = 4,"Flapping" = 0, "Thermal soaring" = 2, "Wind soaring" = 1)) +
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7))

dev.off()

png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/distr_plot_gps&range15.png", width = 6.5, height = 10, units = "in", res = 300)

ggplot(wind_raw, aes(x = wind_speed_ms, y = species_f)) + 
  stat_density_ridges(data = wind_raw[wind_raw$wind_data == "range",], color = "#A9A9A9", fill = "#e6e3e3",
                       jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1.3, point_alpha = 0.8, size = 0.25,
                       calc_ecdf = TRUE, panel_scaling = F, alpha = 0.5,
                       scale = 1.5) +
   new_scale_fill() +
   stat_density_ridges(data = wind_raw[wind_raw$wind_data == "gps_pts",], aes( fill = stat(x)),
                       jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1.3, point_alpha = 0.6, size = 0.25,
                       geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, 
                       scale = 1.5) +
   scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   scale_x_continuous(limits = c(-1, 28)) +
   geom_point(data = wind_raw, aes(x = -0.9, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
   scale_shape_manual(values = c("Dynamic soaring" = 4,"Flapping" = 0, "Thermal soaring" = 2, "Wind soaring" = 1)) +
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7))

dev.off()

#-------------------------------------------------------------------------------------

#the above with the breeding rnage wind densities

ggplot(sample[sample$wind_data == "gps_pts",])+ 
   stat_density_ridges(aes(wind_speed_ms, y = species_f, fill = stat(x)), jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1.3, point_alpha = 0.8, size = 0.25,
                       geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, #only the median line
                       quantiles = 0.5, quantile_lines = T, scale = 3) +
   stat_density_ridges(aes(wind_speed_ms, y = species_f, fill = "gray"), jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 1.3, point_alpha = 0.8, size = 0.25,
                       geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, #only the median line
                       quantiles = 0.5, quantile_lines = T, scale = 3) +
   scale_x_continuous(limits = c(-1, 28)) +
   scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   geom_point(data = sample, aes(x = -0.9, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
   scale_shape_manual(values = c("Dynamic soaring" = 4,"Flapping" = 0, "Thermal soaring" = 2, "Wind soaring" = 1)) +
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7))

 #############################
#ridgeline plots
 ggplot(sample) +
   geom_density_ridges(aes(x = wind_speed_ms, y = species_f,
                           group = interaction(wind_data_f, species_f),
                           fill = wind_data_f, color = wind_data_f),
                       jittered_points = TRUE,
                       position = position_points_jitter(width = 0.05, height = 0),
                       point_shape = '|', point_size = 1.5, point_alpha = 1, alpha = 0.5,
                       rel_min_height = .01, show.legend = F) +
   scale_fill_manual(values = c("#A9A9A9",clr)) +
   scale_color_manual(values = c( "#A9A9A9",clr)) +
   geom_point(data = sample, aes(x = -0.9, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
   scale_x_continuous(limits = c(-1, 30)) +
   scale_shape_manual(values = c(4,0,2,1)) +
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7)) +
   theme(plot.margin = margin(1,0.3,0,0, "cm"))
 
 
 #save
 png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/distr_plot4.png", width = 6.5, height = 10, units = "in", res = 300)
 ggplot(wind_raw) +
   geom_density_ridges(aes(x = wind_speed_ms, y = species_f,
                           group = interaction(wind_data_f, species_f),
                           fill = wind_data_f, color = wind_data_f),
                        jittered_points = TRUE,
                       position = position_points_jitter(width = 0.05, height = 0),
                       point_shape = '|', point_size = 1.1, point_alpha = 1, alpha = 0.4,
                       rel_min_height = .01, show.legend = F) +
   scale_fill_manual(values = c("#A9A9A9", clr)) +
   scale_color_manual(values = c("#A9A9A9", clr)) +
   geom_point(data = wind_raw, aes(x = -0.9, y = species_f, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
   scale_x_continuous(limits = c(-1, 28)) +
   scale_shape_manual(values = c(4,0,2,1)) +
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7)) +
   theme(plot.margin = margin(1.3,0.3,0,0, "cm"))
 
 dev.off()
 
 
#original plot with jittered points
png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/distr_plot5.png", width = 7.5, height = 7.5, units = "in", res = 300)
 
 ggplot(wind_raw[wind_raw$wind_data == "gps_pts",], aes(x = wind_speed_ms, y = species, fill = stat(x))) + 
    geom_point(data =  wind_raw[wind_raw$wind_data == "range",], aes(x = wind_speed_ms, y = species, shape = "|"), size = 0.8, color = "gray") +
   stat_density_ridges(jittered_points = TRUE, rel_min_height = .01,
                       point_shape = "|", point_size = 0.8, point_alpha = 0.5, size = 0.25,
                       geom = "density_ridges_gradient", calc_ecdf = TRUE, panel_scaling = F, #only the median line
                       quantiles = 0.5, quantile_lines = T, scale = 3) +
   geom_point(data = wind_raw[wind_raw$wind_data == "gps_pts",], aes(x = -0.8, y = species, shape = flight_style_F), size = 1.8, stroke = 0.4, color = clr) +
    scale_shape_manual(values = c("Dynamic soaring" = 4, "Flapping" = 0, "Thermal soaring" = 2, "Wind soaring" = 1)) +
    #scale_shape_manual(values = c(4,0,2,1)) +
   scale_x_continuous(limits = c(-0.8, 28)) +
   scale_fill_gradientn(colours = alpha(oce::oceColorsPalette(120), alpha = 0.8), limits = c(0,23), 
                        na.value = "white", guide = 'none') +
   
   labs(y = "", x = expression("Wind speed (m s"^-1*")")) +
   theme_minimal() +
   guides(shape = guide_legend("Flight style:")) +
   theme(legend.position = "bottom",legend.title = element_text(size = 10), 
         legend.text=element_text(size = 7))
 
 dev.off()
 
# Fig 3: Wind stills -------------------------------------
 
 
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


png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper prep/figs/wind_lms.png", width = 7, height = 5.5, units = "in", res = 300)

ggplot(long_df, aes(x = wing.loading..Nm.2., y = wind_speed, col = wind_source, fill = wind_source)) +
  geom_smooth(aes(y = wind_speed, group = wind_source), method = "lm", alpha = .1, level = .95) + #95% standard error
  geom_point(aes(y = wind_speed, shape = flight_style), size = 1.5, stroke = 0.8) +
  labs(x = expression("Wing loading (Nm"^-2*")"),
       y = expression("Wind speed (m s"^-1*")")) +
  scale_shape_manual(values = c(4,0,2,1)) + #filled points: c(15, 17, 19)
  scale_color_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range mediam")) +
  scale_fill_manual(values = c(clr2, "#A9A9A9", clr), name = "Wind data:", labels = c("Max encountered", "Breeding range max", "Breeding range mediam")) +
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

 
 