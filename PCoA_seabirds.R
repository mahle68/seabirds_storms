#PCoA analysis for seabirds
# Apr 13. 2021. Radolfzell, DE.
#Elham Nourani, PhD.
#starting here: https://archetypalecology.wordpress.com/2018/02/19/principal-coordinates-analysis-pcoa-in-r/
#:~:text=Principal%20coordinates%20analysis%20(PCoA%3B%20also,dis)similarity%20matrix%20as%20input.
#and https://ourcodingclub.github.io/tutorials/ordination/

library(vegan)
library(ape)
library(tidyverse)

setwd("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms")

#annotated tracking data
load("R_files/ssf_input_annotated_60_30_all.RData") # ann_50_all  from all_species_prep_2021.R
species <- read.csv("/home/mahle68/ownCloud/Work/Projects/seabirds_and_storms/data/final_datasets.csv")

winds <- ann_50_all %>% 
  filter(used == 1) %>% 
  group_by(sci_name) %>% 
  summarise(max_wind = max(wind_speed),
            min_wind = min(wind_speed),
            var_wind = var(wind_speed))

data <- species[species$scientific.name %in% unique(ann_50_all$sci_name),] %>% 
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

