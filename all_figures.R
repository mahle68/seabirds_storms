#figures for seabirds and storms manuscript.
#follows from seabirds_public.R

#Aug 22, 2022
#Elham Nourani, PHD

library(tidyverse)
library(ggridges)

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

 

 
long_df <- lm_input %>% 
   group_by(species) %>% 
   arrange(desc(max_wind_ms)) %>% 
   slice(1) %>% 
   pivot_longer(cols = c("max_wind_ms", "range_median", "range_max"),
                names_to = "wind_source",
                values_to = "wind_speed")
 
 clr2 <- oce::oceColorsPalette(120)[96]

ggplot(long_df, aes(x = wing.loading..Nm.2., y = wind_speed, col = wind_source, fill = wind_source)) +
  geom_smooth(aes(y = wind_speed, group = wind_source), method = "lm", alpha = .1) +
  geom_point(aes(y = wind_speed, shape = flight_style), size = 2, stroke = 0.8) +
  labs(x = expression("Wing loading (Nm"^-2*")"),
       y = "Wind speed") +
  scale_shape_manual(values = c(4,0,2,1)) + #filled points: c(15, 17, 19)
  scale_color_manual(values = c(clr2, "#A9A9A9", clr)) +
  scale_fill_manual(values = c(clr2, "#A9A9A9", clr)) +
  theme_minimal() + #theme_ipsum() looks better
  theme(axis.title.y = element_text(size = 13)) +
  guides(shape = guide_legend("Flight style:"))
 



 