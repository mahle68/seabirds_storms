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
