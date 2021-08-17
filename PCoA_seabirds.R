#PCoA analysis for seabirds
# Apr 13. 2021. Radolfzell, DE.
#Elham Nourani, PhD.
#starting here: https://archetypalecology.wordpress.com/2018/02/19/principal-coordinates-analysis-pcoa-in-r/:~:text=Principal%20coordinates%20analysis%20(PCoA%3B%20also,dis)similarity%20matrix%20as%20input.
#and https://ourcodingclub.github.io/tutorials/ordination/

library(vegan)
library(ape)
library(tidyverse)
library(lme4)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization. to install: if(!require(devtools)) install.packages("devtools"); #devtools::install_github("kassambara/factoextra")
#detach(package:plyr)
library(corrplot)
library(corrr)


setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms")


## STEP 1: PCA ####
#species-specific
#------------------ prep input data 
species <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/final_datasets.csv") %>% 
  filter(scientific.name %in% c( "Ardenna gravis", "Diomedea dabbenena", "Diomedea exulans", "Fregata magnificens", "Fregata minor", "Morus bassanus", "Morus capensis", "Phaethon lepturus", "Phaethon rubricauda", 
                          "Phoebastria irrorata", "Phoebetria fusca", "Procellaria cinerea", "Pterodroma incerta", "Pterodroma mollis",  "Sula dactylatra", "Sula granti",  "Sula sula", "Thalassarche chlororhynchos"))

morph <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/morphology_data/morphometrics.csv")

data <- species %>% 
  dplyr::select(c(1:3)) %>%
  inner_join(morph[,c(2:7)], by = "species") %>% 
  column_to_rownames("species") %>% 
  dplyr::select(-"scientific.name")

save(data, file = "R_files/morph_18spp.RData")

#------------------ wing loading and aspect ratio correlation
cor(data[-12,c("wing.loading..Nm.2.", "aspect.ratio")]) #0.77

data_z <- scale(data[-12,c("wing.loading..Nm.2.", "aspect.ratio")])

#------------------ run pca and investigate results
##pca in datacamp
#https://www.datacamp.com/community/tutorials/pca-analysis-r

pca_out <- prcomp(data_z, center = F, scale. = F) 

summary(pca_out)

library(ggbiplot)
ggbiplot(pca_out, labels = rownames(data)[-12], groups = data$flight.type[-12], ellipse = T) +
  theme_minimal()+
  theme(legend.position = "bottom")

#how much variation in the original data each PC accounts for
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
pca.var <- pca_out$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) #https://www.youtube.com/watch?v=0Jp4gsfOLMs


var <- get_pca_var(pca_out)
corrplot(var$cos2, is.corr=FALSE)

summary(pca_out)

fviz_pca_var(pca_out, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#contribution of variables to PCs
fviz_contrib(pca_out, choice = "var", axes = 1, top = 10) # contribution to the first PC is 50-50??
corrplot(var$contrib, is.corr=FALSE)   

#extract the PC values to use in the LM
loadings <- pca_out$rotation
axes <- as.data.frame(predict(pca_out, newdata = data_z))

#append to data
data_spp <- rownames_to_column(data, var = "species")

data_pca <-  axes %>%
  rownames_to_column("species") %>% 
  full_join(data_spp, by = "species")

save(data_pca, file = "R_files/data_morph_PCs.RData") #two variables: wing loading and aspect ratio

## STEP 2: summarize wind and timing of breeding ####
#species-colony-specific

#annotated tracking data (one hourly sub-sample; flying only; breeding only (as far as I know); adults only (as far as I know :/ ))
files_ls <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/all_spp_for_pcoa/", pattern = ".csv",recursive = T, full.names = T)

ann <- lapply(files_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  drop_na(ECMWF.ERA5.SL.Sea.Surface.Temperature) %>%   #remove points over land
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.,
         v10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component.,
         wh = ECMWF.ERA5.SL.Significant.Wave.Height,
         airp = ECMWF.ERA5.SL.Mean.Sea.Level.Pressure) %>%
  mutate(delta_t = sst - t2m,
         wind_speed_ms = sqrt(u10m^2 + v10m^2),
         wind_speed_kmh = sqrt(u10m^2 + v10m^2) * 3.6,
         month = month(timestamp),
         yday = yday(timestamp))

save(ann, file = "R_files/ann_18spp.RData")


wind_br <- ann %>% 
  group_by(sci_name, colony.name) %>% #make sure plyr is detached: detach("package:plyr", unload=TRUE)
  summarise(min_wind = min(wind_speed_kmh),
            avg_wind = mean(wind_speed_kmh),
            median_wind = median(wind_speed_kmh),
            quant95_wind = quantile(wind_speed_kmh, probs = 0.95),
            max_wind = max(wind_speed_kmh),
            colony.lat = head(colony.lat,1),
            colony.long = head(colony.long,1),
            median_breeding_m = median(month),
            median_breeding_yday = median(yday))

save(wind_br, file = "R_files/wind_br_18spp.RData")

## STEP 2: summarize timing of breeding ####
#species-colony-specific
# (median of all timestamps)
files_ls <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/all_spp_for_pcoa/", pattern = ".csv",recursive = T, full.names = T)

timing <- lapply(files_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  drop_na(ECMWF.ERA5.SL.Sea.Surface.Temperature) %>%   #remove points over land
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  mutate(month = month(timestamp),
         yday = yday(timestamp)) %>% 
  group_by(sci_name, colony.name) %>% 
  summarise(median_breeding_m = median(month),
            median_breeding_yday = median(yday))

#add wind and morphology variables
lm_input <- timing %>% 
  left_join(species[,c(1,2)], by = c("sci_name" = "scientific.name")) %>% 
  full_join(data_pca, by = "species") %>% 
  dplyr::select(c(1,5,8,2,14,15,3,4,16:20,9:13,6,7)) 

save(lm_input, file = "R_files/lm_input_20spp_col.RData")

#also save as csv and send to Emily
write.csv(lm_input[,-c(19,20)], "R_files/wind_data.csv")


## STEP 4: LM and GAM ####

#on the one-row-per-species dataset
load("R_files/data_w_PCs.RData") #data_pca

data_pca[data_pca$flight.type == "gliding-soaring / shearing", "flight.type"] <- "dynamic soaring"

#LM




#z-transform and correlation
data_sc <- data_pca %>% 
  mutate_at(c("body.mass..kg.","wing.span..m.","wing.area..m2.","wing.loading..Nm.2.","aspect.ratio", "colony.lat", "colony.long"),
          list(z = ~as.numeric(scale(.)))) 

data_sc %>% 
dplyr::select(c("body.mass..kg.","wing.span..m.","wing.area..m2.","wing.loading..Nm.2.","aspect.ratio")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #PC1 and PC2 are correlated


m1 <- lm(max_wind ~ colony.long + colony.lat + wing.loading..Nm.2._z, data = data_sc) #wing loading and aspect ratio are correlated
m1b <- lm(max_wind ~ colony.lat + wing.loading..Nm.2._z, data = data_sc) #wing loading and aspect ratio are correlated


m2 <- lm(max_wind ~ colony.long + colony.lat + PC1, data = data_sc) 


m3 <- lm(max_wind ~ colony.long + colony.lat + PC1, data = data_sc) #PC1 is only wing loading and aspect ratio



#plot residual vs. the values
plot(data_sc$max_wind[-18], resid(m1))
plot(data_sc$colony.lat[-18], resid(m1))

m2 <- lm(max_wind ~ PC1, data = data_sc)

m3 <- lmer(max_wind ~ avg_wsp + avg_mass + (1|colony.lat), data = data)

m2 <- lm(max_wind ~ flight.type * avg_mass, data = data) #Rsquared = 0.50

m4 <- lm(max_wind ~ flight.type, data = data) #54% of variation is explained

m5 <-  lm(max_wind ~ colony.long + colony.lat, data = data) #Rsquared = 0.3465 
  
boxplot(data$max_wind ~ data$flight.type)


library(mgcv)
all_data$yday <- yday(all_data$timestamp)


g1 <- gamm(wind_speed_ms ~ s(location.lat, location.long, k = 100) +
             s(yday, bs = "cc") +
             avg_mass + flight.type , method = "REML", data = all_data) #, 
#weights = varPower(form = ~ lat))


g2 <- gam(max_wind ~ s(colony.lat_z, colony.long_z),
          data = data_sc)


###############################
## STEP 2: clustering #### 
#do i need this step?

#https://uc-r.github.io/kmeans_clustering

load("R_files/morph_wind_20spp_col.RData") #data

#data_z <- scale(data[-12,c(8:12)])
#data_z <- scale(data[-12,c(8,11,12)]) #only mass, wing loading and aspect ratio

#Only include variables related to wing shape: wing loading, aspect ratio, wing span
#data_z <- scale(data[-12,c(9,11,12)]) #remove great shearwater. no wing loading info
data_z <- scale(data[-12,c(11,12)])
#data_z <- scale(data[-12,c(8,11,12)])


#optimal number of clusters
fviz_nbclust(data_z, kmeans, method = "wss") #elbow method
fviz_nbclust(data_z, kmeans, method = "silhouette") #silhouette method


distance <- get_dist(data_z, method = "manhattan")

fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))


k2 <- kmeans(data_z, centers = 2, nstart = 25)
str(k2)

fviz_cluster(k2, data = data_z, ggtheme = theme_bw())


k3 <- kmeans(data_z, centers = 3, nstart = 25)
str(k3)

fviz_cluster(k3, data = data_z, ggtheme = theme_bw())

k4 <- kmeans(data_z, centers = 4, nstart = 25)
str(k4)

fviz_cluster(k4, data = data_z, ggtheme = theme_bw())

k5 <- kmeans(data_z, centers = 5, nstart = 25)
str(k5)

fviz_cluster(k5, data = data_z, ggtheme = theme_bw())



k6 <- kmeans(data_z, centers = 6, nstart = 25)
str(k6)

fviz_cluster(k6, data = data_z, ggtheme = theme_bw())



save(k4, file = "R_files/clustering_4k.R")
save(k5, file = "R_files/clustering_5k.R")




#clustering for wind and geography
data_w <- scale(data[-12,c(2:4)]) #wind
fviz_nbclust(data_w, kmeans, method = "silhouette") #silhouette method


dist_w <- get_dist(data_w, method = "manhattan")


k5_w <- kmeans(data_w, centers = 2, nstart = 25)
str(k5_w)

fviz_cluster(k5_w, data = data_z, ggtheme = theme_bw())

