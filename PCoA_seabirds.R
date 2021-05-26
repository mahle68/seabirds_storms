#PCoA analysis for seabirds
# Apr 13. 2021. Radolfzell, DE.
#Elham Nourani, PhD.
#starting here: https://archetypalecology.wordpress.com/2018/02/19/principal-coordinates-analysis-pcoa-in-r/:~:text=Principal%20coordinates%20analysis%20(PCoA%3B%20also,dis)similarity%20matrix%20as%20input.
#and https://ourcodingclub.github.io/tutorials/ordination/

library(vegan)
library(ape)
library(tidyverse)
library(ggbiplot)
library(lme4)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization. to install: if(!require(devtools)) install.packages("devtools")

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms")

rsd <- function(x){
  cv <- sd(x, na.rm = T)/abs(mean(x, na.rm = T))
  rsd <- cv*100
  return(rsd)
}

## STEP 1: data prep ####
#annotated tracking data (one hourly sub-sample; flying only; breeding only (as far as I know); adults only (as far as I know :/ ))
files_ls <- list.files("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/all_spp_for_pcoa/", pattern = ".csv",recursive = T, full.names = T)

pca_input <- lapply(files_ls, read.csv, stringsAsFactors = F) %>% 
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
         wind_speed_kmh = sqrt(u10m^2 + v10m^2) * 3.6)

save(pca_input, file = "R_files/pcoa_input_18spp.RData")


species <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/final_datasets.csv")

winds <- pca_input %>% 
  group_by(sci_name) %>% #make sure plyr is detached
  summarise(max_wind = max(wind_speed_kmh),
            min_wind = min(wind_speed_kmh),
            rsd_wind = rsd(wind_speed_kmh),
            quant_95 = quantile(wind_speed_kmh, probs = 0.95),
            colony.lat = head(colony.lat,1),
            colony.long = head(colony.long,1))

data <- species[species$scientific.name %in% unique(pca_input$sci_name),] %>% 
  dplyr::select(c(1:5)) %>%
  separate(wing.span, c("min_wspn","max_wspn")) %>% 
  separate(body.mass, c("min_mass","max_mass")) %>% 
  mutate_at(c("min_wspn","max_wspn","min_mass","max_mass"), as.numeric) %>% 
  full_join(winds, by = c("scientific.name" = "sci_name")) %>%
  rowwise () %>% 
  mutate(avg_wsp = mean(c(min_wspn, max_wspn), na.rm = T),
         avg_mass = mean(c(min_mass, max_mass), na.rm = T)) %>% 
  ungroup() %>% 
  column_to_rownames("species") %>% 
  dplyr::select(-c("scientific.name","min_wspn","max_wspn","min_mass","max_mass"))

save(data, file = "R_files/pcoa_18spp.RData")


## STEP 2: PCA ####

#try pca (only numeric variables)
PCA <- rda(data[,-1], scale = T)
# Now plot a bar plot of relative eigenvalues. This is the percentage variance explained by each axis
barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig)) 
# How much of the variance in our dataset is explained by the first principal component?

# Calculate the percent of variance explained by first two axes
sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2]) # 79%, this is ok.
# Also try to do it for the first three axes

# Now, we`ll plot our results with the plot function
plot(PCA)
plot(PCA, display = "sites", type = "points")
plot(PCA, display = "species", type = "text")

# You can extract the species and site scores on the new PC for further analyses:
speciesPCA <- PCA$CA$u # Site scores
varPCA <- PCA$CA$v # Species scores

# In a biplot of a PCA, species' scores are drawn as arrows 
# that point in the direction of increasing values for that variable
biplot(PCA, choices = c(1,2), type = c("text", "points"), xlim = c(-5,10)) # biplot of axis 1 vs 2
biplot(PCA, choices = c(1,3), type = c("text","points")) # biplot of axis 1 vs 3



#try pcoa
dist <- vegdist(data[,-1], method = "bray")

# PCoA is not included in vegan. 
# We will use the ape package instead
PCOA <- pcoa(dist)

# plot the eigenvalues and interpret
barplot(PCOA$values$Relative_eig[1:10])
# Can you also calculate the cumulative explained variance of the first 3 axes?

# Some distance measures may result in negative eigenvalues. In that case, add a correction:
PCOA <- pcoa(dist, correction = "cailliez")

# Plot your results
biplot.pcoa(PCOA)

# You see what`s missing? 
# Indeed, there are no species plotted on this biplot. 
# That's because we used a dissimilarity matrix (sites x sites) 
# as input for the PCOA function. 
# Hence, no species scores could be calculated. 
#However, we could work around this problem like this:
biplot.pcoa(PCOA, data[,-1])


##pca in datacamp
#https://www.datacamp.com/community/tutorials/pca-analysis-r


#all variables
pca_out <- prcomp(data[,-c(1,5,7)], center = T, scale. = T)

summary(pca_out)

ggbiplot(pca_out, labels = rownames(data), groups = data$flight.type, ellipse = T) +
  theme_minimal()+
  theme(legend.position = "bottom")

#try with three variables
data_z <- scale(data[,c(2,6,9)])
pca_out <- prcomp(data_z, center = F) #variables are already scaled


#how much variation in the original data each PC accounts for
#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
pca.var <- pca_out$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1) #https://www.youtube.com/watch?v=0Jp4gsfOLMs

#
var <- get_pca_var(pca_out)
corrplot(var$cos2, is.corr=FALSE)

summary(pca_out)

ggbiplot(pca_out, labels = rownames(data), groups = data$flight.type, ellipse = T) +
  theme_minimal()+
  theme(legend.position = "bottom")

fviz_pca_var(pca_out, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

#contribution of variables to PCs
fviz_contrib(pca_out, choice = "var", axes = 1, top = 10)
corrplot(var$contrib, is.corr=FALSE)   


####

#only morphological measurements
pca_morph <- prcomp(data[,c(8,9)], center = T, scale. = T)


ggbiplot(pca_morph, labels = rownames(data), groups = data$flight.type, ellipse = T) +
  theme_minimal()+
  theme(legend.position = "bottom")

fviz_pca_var(pca_morph, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)

## STEP 3: clustering ####

#https://uc-r.github.io/kmeans_clustering


#devtools::install_github("kassambara/factoextra")

load("R_files/pcoa_18spp.RData") #data

data_z <- scale(data[,c(8,9)])


distance <- get_dist(data_z, method = "pearson")

fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

k2 <- kmeans(data_z, centers = 4, nstart = 25)
str(k2)

fviz_cluster(k2, data = data_z, ggtheme = theme_bw())





  ## STEP 4: LM and GAM ####

#on the one-row-per-species dataset
load("R_files/pcoa_18spp.RData") #data

data <- rownames_to_column(data, var = "common_name")
data[data$flight.type == "gliding-soaring / shearing", "flight.type"] <- "dynamic soaring"

#LM

m1 <- lm(max_wind ~ colony.long + colony.lat + avg_wsp + avg_mass, data = data)

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


g2 <- gam(max_wind ~ s(colony.lat, colony.long) + avg_wsp + avg_mass,
          data = data)
