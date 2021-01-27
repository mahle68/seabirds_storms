# Summary box plots for the seabird project. arranged based on body mass and flight type. 
# Jan 26. 2021. Elham Nourani

library(readxl)
library(plyr)
librar(scales)

#------ PLOT 1: boxplots used vs available ----
#open data: used vs available. 2 hrly steps. no subset for wind strength. filer out sitting if povided

# Peter Ryan data
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/PR_data_ssf_ann.RData") #from ssf_windy_tracks_PR_data.R
PR <- ann_cmpl 
PR$species <- as.factor(PR$species)
levels(PR$species) <- c("Great Shearwater","Tristan Albatross","Sooty Albatross","Grey Petrel", "Atlantic Petrel", "Soft-plumaged Petrel",
                        "Atlantic Yellow-nosed Albatross")
rm(ann_cmpl)
 
#movebank data
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/Mvbnk_2hr_20n_ann.RData") #from mvbnk_ssf.R
MV <- ann_cmpl
MV$species <- as.factor(MV$species)
levels(MV$species) <- c("Scopoli's shearwater","Frigatebird","Manx Shearwater","masked booby")

rm(ann_cmpl)

#wandering albaross data (only flying)
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_all_ssf_ann_2hr.RData") #from WAAL_all.R
AL <- ann_cmpl %>% 
  mutate(species = "Wandering albatross")
rm(ann_cmpl)

#Red_footed booby
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/RFB_data_ssf_ann_2hr.RData") #fom SSF_RFboobie.R
RF <- ann_cmpl %>% 
  mutate(species = "Red-footed booby")

#put all species together
cols <- intersect(intersect(colnames(PR), colnames(MV)), intersect(colnames(AL), colnames(RF)))

all <- PR[,cols] %>% 
  full_join(MV[,cols]) %>% 
  full_join(AL[,cols]) %>% 
  full_join(RF[,cols])

#species in oder of weigh and flight type
species <-  read_excel("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/analyzed data.xlsx")
colnames(species)[c(1,2,3,4,7,8)] <- c("scientific_name", "species","wing_span","min_body_mass", "data_owner", "flight_type")

ds_sp <- species[species$flight_type == "dynamic soaring", "species"]$species
ds <- all[all$species %in% ds_sp,]
ds$species <- factor(ds$species, levels = species[species$flight_type == "dynamic soaring", "species"]$species)

non_ds_sp <- species[species$flight_type != "dynamic soaring", "species"]$species
non_ds <- all[all$species %in% non_ds_sp,]
non_ds$species <- factor(non_ds$species, levels = species[species$flight_type != "dynamic soaring", "species"]$species)

#for dynamic soarers

pdf("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/summar_figs/dyn_soarers.pdf", width = 12, height = 10)
#X11(width = 12, height = 10)
par(mfrow= c(3,1), oma = c(3,0,3,0))
for(i in c("wind_support_kmh", "cross_wind_kmh","wind_speed_kmh")){
  
  boxplot(ds[,i] ~ ds[,"species"], data = ds, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")

  legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")

  boxplot(ds[ds$used == 1,i] ~ ds[ds$used == 1,"species"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ds$species)) - 0.15)
  boxplot(ds[ds$used == 0,i] ~ ds[ds$used == 0,"species"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ds$species)) + 0.15)

}
mtext(paste0("Used and available wind conditions (2-hrs) for dynamic soarers (ordered by body mass)"), side = 3, outer = T, cex = 1.3)
dev.off()

#for non_dynamic soarers
pdf("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/summar_figs/non_dyn_soarers.pdf", width = 12, height = 10)
#X11(width = 12, height = 10)
par(mfrow= c(3,1), oma = c(3,0,3,0))
for(i in c("wind_support_kmh", "cross_wind_kmh","wind_speed_kmh")){
  
  boxplot(non_ds[,i] ~ non_ds[,"species"], data = non_ds, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")
  
  legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  
  boxplot(non_ds[non_ds$used == 1,i] ~ non_ds[non_ds$used == 1,"species"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(non_ds$species)) - 0.15)
  boxplot(non_ds[non_ds$used == 0,i] ~ non_ds[non_ds$used == 0,"species"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(non_ds$species)) + 0.15)
  
}
mtext(paste0("Used and available wind conditions (2-hrs) for non-dynamic soarers (ordered by body mass)"), side = 3, outer = T, cex = 1.3)
dev.off()


#------ PLOT 2: maximum wind speeds encountered ----
#max wind speed data
order <- species #using species in the pipeilne is problematic because variable name is species. so, rename
max_wind <- all %>% 
  filter(used == 1) %>% 
  group_by(species) %>% 
  summarise(maximum_wind_speed = max(wind_speed_kmh,na.rm = T)) %>% 
  full_join(species, by = "species") %>% 
  mutate(species = factor(species, levels = order$species)) #order the factor levels based on body mass


#USE RAW DATA FOR THIS

#load data
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/movebank_data_spli_trip.RData") #MV; called data
data$species <- as.factor(data$species)
levels(data$species) <- c("Scopoli's shearwater","Frigatebird","Manx Shearwater","masked booby")
max_mv <- data %>% group_by(species) %>%  summarize(max_wind_speed = max(windSpeed_kmh)) 
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/waal_all.RData") #waal; called waal
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/Peter_Ryan_data_annotated_SplitTrip_filterWind50kmh.RData") #PR_tot
max_pr <- PR_tot %>% 
  mutate(WindSpeed_kmh = WindSpeed_ms *3.6) %>% 
  group_by(common_name) %>%  
  summarize(max_wind_speed = max(WindSpeed_kmh)) %>% 
  rename(species = common_name)

RFB <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/RFBO_2012_ParamWind.csv", 
         stringsAsFactors = F) #booby

MAX <- rbind(max_pr,max_mv, data.frame(species = c("Red-footed booby", "Wandering albatross"), 
                                       max_wind_speed = c(max(RFB$windSpeed_kmh),max(waal$windSpeed_kmh, na.rm =T))))

species[species$species == "Atlantic Yellow-nosed Albatross","species"] <- "A. Yellow-nosed Albatross"

species <- species %>% 
  full_join(MAX, by = "species") %>% 
  mutate(species = factor(species, levels = order$species)) %>% 
  transform(species = revalue (species,c( "Atlantic Yellow-nosed Albatross" = "A. Yellow-nosed Albatross")))



pdf("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/SSF_process_figures/summar_figs/max_wind_speeds.pdf", width = 7, height = 7)
par(mar = c(11, 5, 4, 2))
plot(0, xlim = c(0,14), ylim = c(30,95), labels = F, tck = 0, ann = F)
points(max_wind_speed ~ species, data = species[species$flight_type == "dynamic soaring",], col = alpha("yellowgreen",0.7), 
       pch = 20, cex = 2)
points(max_wind_speed ~ species, data = species[species$flight_type == "flap-gliding",], col = alpha("firebrick1",0.7), 
       pch = 20, cex = 2)
points(max_wind_speed ~ species, data = species[species$flight_type == "gliding-soaring / shearing",], col = alpha("deepskyblue",0.7), 
       pch = 20, cex = 2)
points(max_wind_speed ~ species, data = species[species$flight_type == "soaring",], col = alpha("darkorchid3",0.7), 
       pch = 20, cex = 2)

axis(side = 1, at = c(seq(1,13,by =1)), line = 0, labels = levels(species$species), 
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