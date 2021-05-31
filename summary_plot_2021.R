#summary plots to show the max wind speed experienced by each seabird species
#update 31.05.2021. Radolfzell, Germany
#Elham Nourani, PhD


library(scales)

#open grouping done using clustering in PCoA seabirds.R
load("R_files/clustering_5k.R") #k5
load("R_files/clustering_4k.R") #k4

k4_df <- as.data.frame(k4$cluster) %>% 
  mutate(species = row.names(.))

ks_df <- as.data.frame(k5$cluster) %>% 
  mutate(species = row.names(.)) %>% 
  full_join(k4_df, by = "species") %>% 
  dplyr::rename(k5 = 1,
                k4 = 3)


load("R_files/morph_wind_18spp.RData") #data
  
data_plot <- data %>% 
  mutate(species = rownames(.)) %>% 
  full_join(ks_df, by = "species") 
  
data_plot[17, "species"] <- "Galapagos albatross"
data_plot[11, "species"] <- "A. Yellow-nosed Albatross"

data_plot$species <- as.factor(data_plot$species)
data_plot <- data_plot[-12,] #leave out great shearwater

#assign colors

#color palette
my_pal <- colorRampPalette(c( "indianred3","darkgoldenrod1","lightpink2", "darkseagreen3","lightslateblue")) #colors palette
Cols_5 <- my_pal(5) #add transparency. 50% is "80". 70% is "B3". 80% is "CC". 90% is "E6"
Cols_4 <- my_pal(4) #paste0(my_pal(4), "80") #add transparency. 50% is "80". 70% is "B3". 80% is "CC". 90% is "E6"

#convert complex numbers to numerics

data_plot$color4 <- as.factor(data_plot$k4)
levels(data_plot$color4) <- Cols_4

data_plot$color5 <- as.factor(data_plot$k5)
levels(data_plot$color5) <- Cols_5



#order data by k groups
data_plot <- data_plot %>% 
  arrange(k5) %>% 
  mutate(color4 = as.character(color4),
         color5 = as.character(color5),
         species = factor(species, unique(species))) #re-order species levels to follow k5 categories

#jpeg("/home/enourani/ownCloud/Work/conferences/seabirds_2021/max_winds.jpg", width = 6, height = 3.5, units = "in", res = 500)


X11(width = 6, height = 3.5)
par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(4.5, 6.7, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(30,90), 
     ylim = c(1,18), xlab = "Max wind speed (km/hr)", ylab = "")

points(species ~ max_wind, data = data_plot, col = data_plot$color5, 
       pch = 20, cex = 2)

axis(side = 1, at = c(30,40,50, 60,70, 80, 90), line = 0, labels = c(30,40,50, 60,70, 80, 90),
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, 
     las = 1)
axis(side= 2, at= c(1:n_distinct(data_plot$species)), #line=-4.8, 
     labels = levels(factor(data_plot$species)),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 


legend("topleft", legend = c("dynamic soaring", "flap-gliding", "gliding-soaring / shearing", "soaring"), 
       col= c("yellowgreen", "firebrick1","deepskyblue","darkorchid3"), pch = 20, cex = 0.8, pt.cex = 1.3, bg = "white", bty = "n" )




dev.off()


#### use wing loading to color the plot

#order data by k groups
data_wl <- data_plot %>% 
  arrange(wing.loading..Nm.2.) %>% 
  mutate(species = factor(species, unique(species))) #re-order species levels to follow k5 categories

X11(width = 6, height = 3.5)
par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.7,0,0),
    mar = c(4.5, 6.7, 0.5, 1),
    bty = "l"
)

plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(30,90), 
     ylim = c(1,18), xlab = "Max wind speed (km/hr)", ylab = "")

points(species ~ max_wind, data = data_wl, 
       pch = 20, cex = 2)


axis(side = 1, at = c(30,40,50, 60,70, 80, 90), line = 0, labels = c(30,40,50, 60,70, 80, 90),
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, 
     las = 1)
axis(side= 2, at= c(1:n_distinct(data_plot$species)), #line=-4.8, 
     labels = levels(factor(data_plot$species)),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 
