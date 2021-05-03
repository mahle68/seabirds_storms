
library(scales)

load("R_files/pcoa_18spp.RData") #data

  
data <- data %>% 
  mutate(species = rownames(.))

data[17, "species"] <- "Galapagos albatross"
data[11, "species"] <- "A. Yellow-nosed Albatross"

data$species <- as.factor(data$species)

data <- data[-3,] #leave white-tailed out for now

jpeg("/home/enourani/ownCloud/Work/conferences/seabirds_2021/max_winds.jpg", width = 6, height = 3.5, units = "in", res = 500)


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

points(species ~ max_wind, data = data[data$flight.type == "dynamic soaring",], col = alpha("yellowgreen",0.7), 
       pch = 20, cex = 2)
points(species ~ max_wind, data = data[data$flight.type == "flap-gliding",], col = alpha("firebrick1",0.7), 
       pch = 20, cex = 2)
points(species ~ max_wind, data = data[data$flight.type == "gliding-soaring / shearing",], col = alpha("deepskyblue",0.7), 
       pch = 20, cex = 2)
points(species ~ max_wind, data = data[data$flight.type == "soaring",], col = alpha("darkorchid3",0.7), 
       pch = 20, cex = 2)

axis(side = 1, at = c(30,40,50, 60,70, 80, 90), line = 0, labels = c(30,40,50, 60,70, 80, 90),
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, 
     las = 1)
axis(side= 2, at= c(1:n_distinct(data$species)), #line=-4.8, 
     labels = levels(factor(data$species)),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2) # text perpendicular to axis label 


legend("topleft", legend = c("dynamic soaring", "flap-gliding", "gliding-soaring / shearing", "soaring"), 
       col= c("yellowgreen", "firebrick1","deepskyblue","darkorchid3"), pch = 20, cex = 0.8, pt.cex = 1.3, bg = "white", bty = "n" )




dev.off()

