#PCoA analysis for seabirds
# Apr 13. 2021. Radolfzell, DE.
#Elham Nourani, PhD.
#starting here: https://archetypalecology.wordpress.com/2018/02/19/principal-coordinates-analysis-pcoa-in-r/
#:~:text=Principal%20coordinates%20analysis%20(PCoA%3B%20also,dis)similarity%20matrix%20as%20input.
#and https://ourcodingclub.github.io/tutorials/ordination/

library(vegan)
library(ape)
library(dplyr)


#annotated tracking data
load("R_files/ssf_input_annotated_60_30_all.RData") # ann_50_all  from all_species_prep_2021.R
species <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/final_datasets.csv")

#create input data. rows are different species, columns are the variables of interest
data <- data.frame(body_mass = rep(NA,15),
                   min_wind = NA,
                   max_wind = NA,
                   avg_wind = NA,
                   var_wind = NA,
                   flight_type = NA,
                   colony_lat = NA,
                   wing_span = NA,
                   row.names = unique(ann_50_all$sci_name))


