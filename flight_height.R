#plot flight height and wind speed for the species that have flight height data
#Nov5.2021. Konstanz, DE.
#Elham Nourani, PhD. 

library(tidyverse)
library(lubridate)
library(scale)

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/")

#manually looked at which species have altitude data. only 3 (4, but the "Frigatebirds breeding at Iguana Island, Panama" is 92% NA)

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/final_list_track_split/MB_FREG_cyclone2012-16_splitTrip_TrackParam.Rdata") #MB_FREG_cyclone2012_16_splitTrip_TrackParam
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/final_list_track_split/MB_WAVEAL_splitTrip.Rdata") #WaveAlba_split
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/final_list_track_split/RFBO_GalapGenovesa_HW_split.Rdata") #RFBO_GalapGenovesa_HW_split

#plot 
plot(height.above.ellipsoid ~ windSpeed_ms, data = WaveAlba_split %>%  filter(height.above.ellipsoid < 1000), main = "waved albatross")

#annotate all with wind 

frigate <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/flight_height/MPIAB PNIC hurricane frigate tracking-5622463111006907301/MPIAB PNIC hurricane frigate tracking-5622463111006907301.csv") %>% 
  mutate(wind_speed_ms = sqrt(ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.^2 + ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component.^2),
         height.raw.n = as.numeric(height.raw)) 

plot(height.raw.n ~ wind_speed_ms, data = frigate  %>%  filter(height.raw.n < 1000), main = "magnificent frigatebird")
plot(height.raw.n ~ wind_speed_ms, data = frigate, main = "magnificent frigatebird")


RFB <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/flight_height/Red footed boobies (Weimerskirch)-6779291706504551434/Red footed boobies (Weimerskirch)-6779291706504551434.csv") %>% 
  mutate(wind_speed_ms = sqrt(ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.^2 + ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component.^2))

plot(height.above.ellipsoid  ~ wind_speed_ms, data = RFB   %>%  filter(height.above.ellipsoid  < 1000), main = "red footed booby")



pdf("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/figures/flight_height.pdf", height = 7, width = 5)

X11(height = 7, width = 5)
par(mfrow = c(3,1))

plot(height.above.ellipsoid ~ windSpeed_ms, data = WaveAlba_split %>%  filter(height.above.ellipsoid < 1000), main = "Galapagos albatross", 
     cex = 1, pch = 20, col = alpha("black", 0.3), xlab = "")
plot(height.raw.n ~ wind_speed_ms, data = frigate  %>%  filter(height.raw.n < 3000), main = "magnificent frigatebird",
     cex = 1, pch = 20, col = alpha("black", 0.3), xlab = "")
plot(height.above.ellipsoid  ~ wind_speed_ms, data = RFB   %>%  filter(height.above.ellipsoid  < 1000), main = "red footed booby",
     cex = 1, pch = 20, col = alpha("black", 0.3))

dev.off()


#update: May2.2022: correlation tests for altitude instead of making figures
#
cor.test(WaveAlba_split$height.above.ellipsoid, WaveAlba_split$windSpeed_ms) #-0.05
cor.test(frigate$height.raw.n, frigate$wind_speed_ms) #-0.04
cor.test(RFB$height.above.ellipsoid, RFB$wind_speed_ms) #0.002
