# script for plotting the globe to add as inset to the seabird storm figures
# Nov12. 2021. Konstanz, DE. Elham Nourani
# https://egallic.fr/en/maps-with-r/

library(rworldmap)
library(dplyr)
library(ggplot2)
library(geosphere)
library(gpclib)

# World map
worldMap <- getMap()
world.points <- fortify(worldMap)
world.points$region <- world.points$id

region_df <- data.frame(x = c(st_bbox(region)[1],st_bbox(region)[3],st_bbox(region)[3],st_bbox(region)[1]),
                        y = c(st_bbox(region)[2],st_bbox(region)[2],st_bbox(region)[4],st_bbox(region)[4]))

world.df <- world.points[,c("long","lat","group", "region")]

worldmap <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45)

worldmap
region_ortho <- st_as_sfc(st_bbox(region)) %>% 
  st_transform(ortho) 
  
#this works: https://egallic.fr/en/maps-with-r/
worldmap <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) +
  #coord_sf(datum = st_crs(ortho)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  coord_map("ortho", orientation = c(wanderer2$location.lat+30, wanderer2$location.long, 0)) +
  geom_sf(data = region_ortho, fill = NA, color = "red", size = 0.5) + 
  theme_void()
worldmap

worldmap <- ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group)) +
  geom_polygon(data = region_df, aes(x = x, y = y), fill = NA, color = "red", size = 0.3) + 
  coord_map("ortho", orientation = c(wanderer2$location.lat+37, wanderer2$location.long, 0)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  theme_minimal() +
  theme(axis.text = element_blank()) +
  labs(x = NULL, y = NULL)
worldmap



# write a function to plot the wind field along with the inset globe
#create main map


plto_with_inset <- function(wind, location, region_sf,  world.df){
  
  region_df <- data.frame(x = c(st_bbox(region_sf)[1],st_bbox(region_sf)[3],st_bbox(region_sf)[3],st_bbox(region_sf)[1]),
                          y = c(st_bbox(region_sf)[2],st_bbox(region_sf)[2],st_bbox(region_sf)[4],st_bbox(region_sf)[4]))
    
plot_w2 <- ggplot() +
  geom_raster(data = wind, aes(x = lon, y = lat, fill = wind_speed))+
  geom_segment(data = wind, 
               aes(x = lon, xend = lon+u10/10, y = lat, 
                   yend = lat+v10/10), arrow = arrow(length = unit(0.12, "cm")), size = 0.3)+
  geom_sf(data = world, fill = "grey85")+
  geom_point(data = location, aes(x = location.long, y = location.lat), 
             size = 2, colour = "red") +
  #coord_sf(xlim = c(-62, 11), ylim =  c(-51, -25))+
  coord_sf(xlim = st_bbox(region_sf)[c(1,3)], ylim = st_bbox(region_sf)[c(2,4)])+
  scale_fill_gradientn(colours = oce::oceColorsPalette(120), limits = c(0,23), 
                       na.value = "white", name = "Speed\n (m/s)")+
  theme_bw()+
  theme(axis.text = element_text(size = 12, colour = 1),
        legend.text = element_text(size = 10, colour = 1), 
        legend.title = element_text(size = 12, colour = 1),
        legend.position = c(0.08,0.23),
        legend.background = element_rect(colour = 1, fill = "white"),
        plot.title = element_text(face = "italic"))+
  labs(x = NULL, y = NULL, title = paste(location$sci_name, paste(location$timestamp, "UTC", sep = " "), sep = "   "))

#create inset map
worldmap_w2 <- ggplot() + 
  geom_polygon(data = region_df, aes(x = long, y = lat, group = group), fill = "grey80") +
  geom_polygon(data = region_w2_df, aes(x = x, y = y), fill = NA, color = "black", size = 0.3) + 
  coord_map("ortho", orientation = c(location$location.lat+37, location$location.long, 0)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.background = element_rect(fill = "white")) +
  labs(x = NULL, y = NULL)
#worldmap_w2


final_w2 <- ggdraw() +
  draw_plot(plot_w2) +
  draw_plot(worldmap_w2, x = 0.69, y = 0.70, width = 0.25, height = 0.25)

}