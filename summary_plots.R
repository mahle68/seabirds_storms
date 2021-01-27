# Summary box plots for the seabird project. arranged based on body mass and flight type. 
# Jan 26. 2021. Elham Nourani

library(readxl)


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
rm(ann_cmpl)

#wandering albaross data (only flying)
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/WAAL_all_ssf_ann_2hr.RData") #from WAAL_all.R
AL <- ann_cmpl
rm(ann_cmpl)

#Red_footed booby
load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/RFB_data_ssf_ann_2hr.RData") #fom SSF_RFboobie.R
RF <- ann_cmpl


#species in oder of weigh and flight type
species <- read_xl
  
#for plotting
ann_cmpl$species <- factor(ann_cmpl$species)


#plot
X11(width = 15, height = 10);par(mfrow= c(3,1), oma = c(0,0,3,0))
for(i in c("wind_support_kmh", "cross_wind_kmh","wind_speed_kmh")){
  #  for(j in c("tmpz", "twz")){ 
  
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")
  
  legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  
  boxplot(ann_cmpl[ann_cmpl$used == 1,i] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0,i] ~ ann_cmpl[ann_cmpl$used == 0,"species"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  #  } 
}
mtext("40-yr averages at each point", side = 3, outer = T, cex = 1.3)

#------ PLOT 2: scatter plot max wind speeds. color = fligh type ----

#open data. all data annoated with wind. no need for available poins

load("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/movebank_data_spli_trip.RData")




