#script for preparing the data for PCoA. No need to subsample the data. just make sure that only flying points are used in the PCoA.
#separate species from different locations.
#Elham Nourani. PhD.
#Apr. 14. 2021. Radolfzell, DE. 




#data from summary_plots.R


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
