#figs for ASAB winter meeting 2022 poster
#Nov 28. 2022
#Elham Nourani, Konstanz, DE


#the globe (from wind_fields.R)

X11(width = 10, height = 10) #in inches


png("/home/enourani/ownCloud/Work/conferences/ASAB2022_winter/figs/the_globe.png", width = 10, height = 10, units = "in", res = 300)

ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group), fill = "black") +
  geom_polygon(data = ext, aes(x = x, y = y), fill = NA, color = "black", size = 0.3) + 
  coord_map("ortho", orientation = c(point$location.lat+37, point$location.long, 0)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  theme_minimal() +
  theme(axis.text = element_blank()) +
  labs(x = NULL, y = NULL) 

dev.off()



ggplot() + 
  geom_polygon(data = world.df, aes(x = long, y = lat, group = group), fill = "black") +
  geom_polygon(data = ext, aes(x = x, y = y), fill = NA, color = "black", size = 0.3) + 
  coord_map("ortho", orientation = c(point$location.lat+37, point$location.long, 0)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  theme_minimal() +
  theme(axis.text = element_blank()) +
  labs(x = NULL, y = NULL) 



