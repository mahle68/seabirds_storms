
library(scales)

load("R_files/pcoa_18spp.RData")

data <- data %>% 
  mutate(species = factor(rownames(.)))

pdf("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/summar_figs/max_winds.pdf", width = 7, height = 7)
par(mar = c(11, 5, 4, 2))
plot(0, xlim = c(0,18), ylim = c(30,95), labels = F, tck = 0, ann = F)
points(max_wind ~ species, data = data[data$flight.type == "dynamic soaring",], col = alpha("yellowgreen",0.7), 
       pch = 20, cex = 2)
points(max_wind ~ species, data = data[data$flight.type == "flap-gliding",], col = alpha("firebrick1",0.7), 
       pch = 20, cex = 2)
points(max_wind ~ species, data = data[data$flight.type == "gliding-soaring / shearing",], col = alpha("deepskyblue",0.7), 
       pch = 20, cex = 2)
points(max_wind ~ species, data = data[data$flight.type == "soaring",], col = alpha("darkorchid3",0.7), 
       pch = 20, cex = 2)

axis(side = 1, at = 1:length(levels(data$species)), line = 0, labels = levels(data$species), 
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, las = 2, font =  3)
axis(side = 2, at = c(30,40,50, 60,70, 80), line = 0, labels = c(30,40,50, 60,70, 80),
     tick = T , col.ticks = 1, col = NA, lty = NULL, tck = -.015, 
     las = 2, font =  3)

legend("topright", legend = c("dynamic soaring", "flap-gliding", "gliding-soaring / shearing", "soaring"), 
       col= c("yellowgreen", "firebrick1","deepskyblue","darkorchid3"), pch = 20, cex = 0.8, pt.cex = 1.3, bg = "white", bty = "n" )

mtext("Max wind speed (km/h)", 2, line = 2.5 ,las = 0, cex = 1, font = 3)
mtext("Maximum wind speed encountered by each species", 3, line = 1.3, cex = 1.1, font = 4)
    mtext("(ordered by body mass)", 3, line = 0.2, cex = 1, font = 4)
dev.off()