#script for plotting the globe
#https://gist.github.com/fzenoni/ef23faf6d1ada5e4a91c9ef23b0ba2c1

# Download earth data first
# https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/physical/ne_110m_land.zip

library(sf)
library(lwgeom)
library(dplyr)
library(ggplot2)
library(mapview)

# Read the data
mini_world <- read_sf("/home/enourani/ownCloud/Work/GIS_files/ne_110m_land/ne_110m_land.shp")

# Define the orthographic projection
# Choose lat_0 with -90 <= lat_0 <= 90 and lon_0 with -180 <= lon_0 <= 180
lat <- 45
lon <- 2
ortho <- paste0('+proj=ortho +lat_0=', lat, ' +lon_0=', lon, ' +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs')

# Define the polygon that will help you finding the "blade"
# to split what lies within and without your projection
circle <- st_point(x = c(0,0)) %>% st_buffer(dist = 6371000) %>% st_sfc(crs = ortho)

# Project this polygon in lat-lon
circle_longlat <- circle %>% st_transform(crs = 4326)

# circle_longlat cannot be used as it is
# You must decompose it into a string with ordered longitudes
# Then complete the polygon definition to cover the hemisphere
if(lat != 0) {
  circle_longlat <- st_boundary(circle_longlat)
  
  circle_coords <- st_coordinates(circle_longlat)[, c(1,2)]
  circle_coords <- circle_coords[order(circle_coords[, 1]),]
  circle_coords <- circle_coords[!duplicated(circle_coords),]
  
  # Rebuild line
  circle_longlat <- st_linestring(circle_coords) %>% st_sfc(crs = 4326)
  
  if(lat > 0) {
    rectangle <- list(rbind(circle_coords,
                            c(X = 180, circle_coords[nrow(circle_coords), 'Y']),
                            c(X = 180, Y = 90),
                            c(X = -180, Y = 90),
                            c(X = -180, circle_coords[1, 'Y']),
                            circle_coords[1, c('X','Y')])) %>% 
      st_polygon() %>% st_sfc(crs = 4326)
  } else {
    rectangle <- list(rbind(circle_coords,
                            c(X = 180, circle_coords[nrow(circle_coords), 'Y']),
                            c(X = 180, Y = -90),
                            c(X = -180, Y = -90),
                            c(X = -180, circle_coords[1, 'Y']),
                            circle_coords[1, c('X','Y')])) %>% 
      st_polygon() %>% st_sfc(crs = 4326)
  }
  
  circle_longlat <- st_union(st_make_valid(circle_longlat), st_make_valid(rectangle))
}

# This visualization shows the visible emisphere in red
ggplot() +
  geom_sf(data = mini_world) +
  geom_sf(data = circle_longlat, color = 'red', fill = 'red', alpha = 0.3)

# A small negative buffer is necessary to avoid polygons still disappearing in a few pathological cases
# I should not change the shapes too much
#visible <- st_intersection(st_make_valid(mini_world), st_buffer(circle_longlat, -0.09)) %>%
#  st_transform(crs = ortho) #this returns an empty polygon, so skip it

visible <- circle_longlat %>%
  st_transform(crs = ortho) 

# DISCLAIMER: This section is the outcome of trial-and-error and I don't claim it is the best approach 
# Resulting polygons are often broken and they need to be fixed
# Get reason why they're broken
broken_reason <- st_is_valid(visible, reason = TRUE)

# First fix NA's by decomposing them
# Remove them from visible for now
na_visible <- visible[is.na(broken_reason),]
visible <- visible[!is.na(broken_reason),]

# Open and close polygons
na_visible <- st_cast(na_visible, 'MULTILINESTRING') %>% 
  st_cast('LINESTRING', do_split=TRUE)
na_visible <- na_visible %>% mutate(npts = npts(geometry, by_feature = TRUE))
# Exclude polygons with less than 4 points
na_visible <- na_visible %>%
  filter(npts >=4) %>%
  dplyr::select(-npts) %>%
  st_cast('POLYGON')

# Fix other broken polygons
broken <- which(!st_is_valid(visible))
for(land in broken) {
  result = tryCatch({
    # visible[land,] <- st_buffer(visible[land,], 0) # Sometimes useful sometimes not
    visible[land,] <- st_make_valid(visible[land,]) %>% 
      st_collection_extract()  
  }, error = function(e) {
    visible[land,] <<- st_buffer(visible[land,], 0)
  })
}

# Bind together the two tables
visible <- rbind(visible, na_visible)

# Final plot
ggplot() +
  geom_sf(data = circle,
          #fill = 'aliceblue') + # if you like the color
          fill = NA) +
  geom_sf(data=st_collection_extract(visible)) +
  coord_sf(crs = ortho)


#####################
library(rworldmap)
library(dplyr)
library(ggplot2)
library(mapproj)
library(sf)
library(rnaturalearth)
## MWE for how sf objects create rendering artifacts in orthographic projection
## and how sf objects to not work within coord_map() in ggplot2

# download world data from natural earth
world <- ne_countries(scale = 'small', returnclass = 'sf')


class(world) # returns "sf" "data.frame"

# if we plot this with the standard plotting we get a world map
ggplot() + 
  geom_sf(data=world) 

# if we now add an orthographic projection the horizontal lines become apparent
ortho <- "+proj=ortho +lat_0=-35 +lon_0=170 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs"
ggplot() + 
  geom_sf(data=world) + 
  coord_sf(crs = ortho)


# now lets try to use the coord_map() function instead of coord_sf
# this results in an error because we need to use coord_sf with geom_sf
ggplot() + 
  geom_sf(data=world) + 
  coord_map("ortho", orientation=c(-35, 175, 0))


# if we want to use the world data with the coord_map() function
# we need to convert it to a dataframe. First we need to convert 
# the sf object to a SpatialPolygonsDataframe so we can fortify it
world_fort <- world %>%
  as("Spatial") %>%
  fortify()

class(world_fort) # returns "data.frame"

ggplot() + 
  geom_polygon(data = world_fort, aes(x = long, y = lat, group = group), fill = "grey60") +
  coord_map("ortho", orientation=c(-35, 175, 0)) +
  theme_bw()


################################################
### World map ------------------------------------
world <- rnaturalearth::ne_countries(scale = 'small', returnclass = 'sf')

# Fix polygons to ortho projection, following from @fzenoni: https://github.com/r-spatial/sf/issues/1050
world  <- st_cast(world, 'MULTILINESTRING') %>%
  st_cast('LINESTRING', do_split=TRUE) %>%
  mutate(npts = npts(geometry, by_feature = TRUE)) %>%
  st_cast('POLYGON')


# Define the orthographic projection
# Choose lat_0 with -90 <= lat_0 <= 90 and lon_0 with -180 <= lon_0 <= 180
lat <- -6
lon <- -57
ortho <- paste0('+proj=ortho +lat_0=', lat, ' +lon_0=', lon, ' +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs')


# globe border
globe <- st_graticule(ndiscr = 10000, margin = 10e-6) %>%
  st_transform(crs = ortho) %>%
  st_convex_hull() %>%
  summarise(geometry = st_union(geometry))  




### Plot  ------------------------------------

temp_plot <- ggplot() +
  geom_sf(data=globe, fill="gray98", color="gray98") +
  geom_sf(data=world, fill="gray90", color="gray80") +
  geom_sf(data=brazil, fill="gray75", color="gray80") +
  geom_sf(data=amazon, fill="#306844", color="#306844") +
  theme_map()


## zoom
# + coord_sf(datum =ortho, xlim = c(st_bbox(brazil2)[[1]], st_bbox(brazil2)[[3]]), ylim = c(st_bbox(brazil2)[[2]], st_bbox(brazil2)[[4]]), expand = FALSE) 
# brazil2 <- st_transform(brazil, crs = ortho)


############################################
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
