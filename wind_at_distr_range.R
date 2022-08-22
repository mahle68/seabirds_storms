#script to plot the wind conditions at the species' distribution ranges
#based on code from Manos in "seabirds_and_storms/from manos/CodeForPlotsAnalysisForElham.R"
#Aug 15, 2022. Konstanz, Germany.
#Elham Nourani

library(ggplot2)

library(ggridges)
library(viridis)


library(forcats)

setwd("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/from manos/")


Group1_WindSpeeds_BySpeciesDitrbitution<-read.csv("SeabirdBreedingDistr_WindSpeed_Group1.csv")

Group2_WindSpeeds_BySpeciesDitrbitution<-read.csv("SeabirdBreedingDistr_WindSpeed_Group2.csv")

Group3_WindSpeeds_BySpeciesDitrbitution<-read.csv("SeabirdBreedingDistr_WindSpeed_Group3.csv")


all_data <-rbind.data.frame(Group1_WindSpeeds_BySpeciesDitrbitution,Group2_WindSpeeds_BySpeciesDitrbitution)

all_wind <-rbind.data.frame(all_data,Group3_WindSpeeds_BySpeciesDitrbitution)


nrow(Group1_WindSpeeds_BySpeciesDitrbitution)

saveRDS(all_wind, file = "/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/R_files/wind_at_distr_17_22.rds")

SpeciesSummaryTableFromElham<-read.csv("wind_wing_data_by_Species_FromEmily.csv")

SpeciesSummaryTableFromElham<-SpeciesSummaryTableFromElham[!duplicated(SpeciesSummaryTableFromElham$sci_name),]


SpeciesSummaryTableFromElham$max_wind<-SpeciesSummaryTableFromElham$max_wind/3.6 

SpeciesSummaryTableFromElham<-SpeciesSummaryTableFromElham[order(SpeciesSummaryTableFromElham$species),]


check4<-SpeciesSummaryTableFromElham$sci_name %in% unique(Group1_WindSpeeds_BySpeciesDitrbitution$Sci_name)

length(check4[check4==FALSE])


###############################################################################################################
#THE FOLLOWING FOR LOOP CODE MUST RUN IRRESPECTIVE OF WHAT YOU ARE PLANNING TO PLOT
#it finds stats of wind speed per species and assigns species name in the data using the table of species with wing loading etc#

WindSpeedThreshold<-10 #to find the percentage of wind speeds above this threshold

for (ID in unique(Group1_WindSpeeds_BySpeciesDitrbitution$Sci_name)) {
  
  
  Group1_WindSpeeds_BySpeciesDitrbitution$species[Group1_WindSpeeds_BySpeciesDitrbitution$Sci_name==ID]<-
    SpeciesSummaryTableFromElham$species[SpeciesSummaryTableFromElham$sci_name==ID][1]
  
  SpeciesSummaryTableFromElham$MaxERA5WndSpeed[SpeciesSummaryTableFromElham$sci_name==ID]<-
    max(Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed[Group1_WindSpeeds_BySpeciesDitrbitution$Sci_name==ID],na.rm = TRUE)
  
  
  SpeciesSummaryTableFromElham$MedianERA5WndSpeed[SpeciesSummaryTableFromElham$sci_name==ID]<-
    median(Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed[Group1_WindSpeeds_BySpeciesDitrbitution$Sci_name==ID],na.rm = TRUE)
  
  
  SpeciesSummaryTableFromElham$MeanERA5WndSpeed[SpeciesSummaryTableFromElham$sci_name==ID]<-
    mean(Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed[Group1_WindSpeeds_BySpeciesDitrbitution$Sci_name==ID],na.rm = TRUE)
  
  temp<-data.frame("WndSpeed"=Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed[Group1_WindSpeeds_BySpeciesDitrbitution$Sci_name==ID])
  
  SpeciesSummaryTableFromElham$PercentOfWindSpeedsAboveThreshold[SpeciesSummaryTableFromElham$sci_name==ID]<-
    round( (length(temp$WndSpeed[temp>WindSpeedThreshold])/nrow(temp))*100 , digits=1)
  
  
  
}

Group1_WindSpeeds_BySpeciesDitrbitution$species<-factor(Group1_WindSpeeds_BySpeciesDitrbitution$species)
levels(Group1_WindSpeeds_BySpeciesDitrbitution$species)

Group1_WindSpeeds_BySpeciesDitrbitution<-Group1_WindSpeeds_BySpeciesDitrbitution[order(Group1_WindSpeeds_BySpeciesDitrbitution$species),]



SpeciesSummaryTableFromElham$flight_style_F <- factor(SpeciesSummaryTableFromElham$flight.type)



###############################################################################################################
#The following code plots and estimates linear models for the summary statistics of wind speed per species#

#setwd("C:/Users/lembi/Desktop")

png(filename = "A) MaxERA5WndSpeed against WingLoading Groups1And2.png", #for image: = 3508 x 2480 Pixel
    width = 22.7, height = 20.9, units = "cm", #pointsize = 12,
    bg = "white",  res = 1000)


      ggplot(SpeciesSummaryTableFromElham, aes(x=SpeciesSummaryTableFromElham$wing.loading..Nm.2., y=SpeciesSummaryTableFromElham$MaxERA5WndSpeed)) +
        geom_point(aes(shape = SpeciesSummaryTableFromElham$flight_style_F), size = 2.5, stroke = 1.2,col="black") +
        geom_smooth(method="lm" , formula= y~x, color="black", fill="lightblue", se=TRUE) +
        ylab(expression("Max Wind speed (Breeding)" ~ (m~s^{-1}) )) +
        xlab(expression("Wing loading" ~ (n~m^{-2}) ))+
        labs(shape='Flight style')+
        theme_classic()

dev.off()

mod1<-lm(
  MaxERA5WndSpeed~wing.loading..Nm.2.
  ,data = SpeciesSummaryTableFromElham)

summary(mod1)

png(filename = "B) MaxERA5WndSpeed against Max WindSpeed At Birds Locations Groups1And2.png", #for image: = 3508 x 2480 Pixel
    width = 22.7, height = 20.9, units = "cm", #pointsize = 12,
    bg = "white",  res = 1000)


      
      ggplot(SpeciesSummaryTableFromElham, aes(x=SpeciesSummaryTableFromElham$max_wind, y=SpeciesSummaryTableFromElham$MaxERA5WndSpeed)) +
        geom_point(aes(shape = SpeciesSummaryTableFromElham$flight_style_F), size = 2.5, stroke = 1.2,col="black") +
        geom_smooth(method="lm" , formula= y~x, color="black", fill="lightblue", se=TRUE) +
        ylab(expression("Maximum Wind speed (Breeding)" ~ (m~s^{-1}) )) +
        xlab(expression("Maximum Wind speed (GPS)" ~ (m~s^{-1}) ))+
        labs(shape='Flight style')+
        theme_classic()

dev.off()


mod12<-lm(
  MaxERA5WndSpeed~max_wind
  ,data = SpeciesSummaryTableFromElham)

summary(mod12)

png(filename = "C) MedianERA5WndSpeed against WingLoading Groups1And2.png", #for image: = 3508 x 2480 Pixel
    width = 22.7, height = 20.9, units = "cm", #pointsize = 12,
    bg = "white",  res = 1000)



    ggplot(SpeciesSummaryTableFromElham, aes(x=SpeciesSummaryTableFromElham$wing.loading..Nm.2., y=SpeciesSummaryTableFromElham$MedianERA5WndSpeed)) +
      geom_point(aes(shape = SpeciesSummaryTableFromElham$flight_style_F), size = 2.5, stroke = 1.2,col="black") +
      geom_smooth(method="lm" , formula= y~x, color="black", fill="lightblue", se=TRUE) +
      ylab(expression("Median Wind speed (Breeding)" ~ (m~s^{-1}) )) +
      xlab(expression("Wing loading" ~ (n~m^{-2}) ))+
      labs(shape='Flight style')+
      theme_classic()

dev.off()



mod2<-lm(
  MedianERA5WndSpeed~wing.loading..Nm.2.
  ,data = SpeciesSummaryTableFromElham)

summary(mod2)

png(filename = "D) MedianERA5WndSpeed against Max Wind Speed At Birds Locations Groups1And2.png", #for image: = 3508 x 2480 Pixel
    width = 22.7, height = 20.9, units = "cm", #pointsize = 12,
    bg = "white",  res = 1000)



      ggplot(SpeciesSummaryTableFromElham, aes(x=SpeciesSummaryTableFromElham$max_wind, y=SpeciesSummaryTableFromElham$MedianERA5WndSpeed)) +
        geom_point(aes(shape = SpeciesSummaryTableFromElham$flight_style_F), size = 2.5, stroke = 1.2,col="black") +
        geom_smooth(method="lm" , formula= y~x, color="black", fill="lightblue", se=TRUE) +
        ylab(expression("Median Wind speed (Breeding)" ~ (m~s^{-1}) )) +
        xlab(expression("Maximum Wind speed (GPS)" ~ (m~s^{-1}) ))+
        labs(shape='Flight style')+
        theme_classic( )

dev.off()

mod3<-lm(
  MedianERA5WndSpeed~max_wind
  ,data = SpeciesSummaryTableFromElham)

summary(mod3)



###############################################################################################################
#The following code plots one log histogram of wind speed per species in the folder that is set as working directory#
#The x axsis of  the histograms is always the same but the y axis varies#

setwd("log_histograms/sameX_axis")

maxValueAllDataset<-ceiling( max(Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed))
for (ID in unique(Group1_WindSpeeds_BySpeciesDitrbitution$species)) {
  
        tempstr<-paste(ID,"log-Histogram.png")
        
        
        png(filename = tempstr, #for image: = 3508 x 2480 Pixel
            width = 22.7, height = 20.9, units = "cm", #pointsize = 12,
            bg = "white",  res = 1000)
        
        
        par(mar = c(6.0, 6.5, 1.5, 1.0)) # Set the margin on each side. sets the bottom, left, top and right 
        
        maxValue<-ceiling( max(Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed[Group1_WindSpeeds_BySpeciesDitrbitution$species==ID]))
        
        if(maxValue<max(Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed)){
          
          temp<-c(Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed[Group1_WindSpeeds_BySpeciesDitrbitution$species==ID],32)
          
          
        }else{
          temp<-Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed[Group1_WindSpeeds_BySpeciesDitrbitution$species==ID]
        }
        
        #USE THE LINE BELOW TO PLOT LOG HISTOGRAMS WITH X AXIS SPECIFIC TO EACH SPECIES
        # r <- hist(Group1_WindSpeeds_BySpeciesDitrbitution$WindSpeed[Group1_WindSpeeds_BySpeciesDitrbitution$species==ID],breaks = maxValue)
        
        #USE THE LINE BELOW TO PLOT ALL HISTOGRAMS WITH THE SAME X AXIS AND NUMBER OF BINS
        r <- hist(temp,breaks = maxValueAllDataset)
        
        logCounts<-log(r$counts)
        logCounts<- ifelse(is.infinite(logCounts),0,logCounts )
        
        barplot(logCounts, border="white", col="lightblue", names.arg=  r$breaks[-1]
                ,axis.lty=1
                ,xlab="",cex.lab=2, cex.axis =1.75, cex.names=1.75,las=1,axes=TRUE,ylab=""
                , axisnames=TRUE )
        
        par(mar = c(4.8, 5.2, 1.5, 1.0)) # Set the margin on each side. sets the bottom, left, top and right 
        
        mtext("Frequency (log)", side = 2, line = 2.5,col="black",cex=2)  
        
        mtext(expression("Wind speed " ~ (m~s^{-1}) ), side = 1, line = 2.5,col="black",cex=2)  
        #  axis(2,at=seq(min(log (r$counts)),max(log (r$counts)),by=1),labels=FALSE)
        
        dev.off()
        
}


###############################################################################################################
#The following code plots a single figure with log histograms of wind speed for each species#
#Both axes y and x are the same, but the freuqnecies are not normalised by the number data points of each species#

thm2 <- theme_classic() + theme(text = element_text(size = 11),
                                axis.title.x = element_text(size = 11, colour = "black",angle=0,margin = margin(1*c(0.4,0,0,0), unit = "cm"))
                                ,axis.text.x = element_text(size = 10, colour = "black",angle=0,margin = margin(1*c(0.2,0,0,0), unit = "cm"))
                                # ,axis.text.y =  element_text(size = 16, colour = "black",angle=0,margin = margin(1*c(0,0,0,0), unit = "cm"))
                                ,legend.position = "none",
                                # axis.title.x = element_blank(),
                                # axis.text = element_blank(),
                                # axis.ticks.x = element_blank(),
                                axis.text.y = element_blank(),
                                axis.title.y = element_blank(),
                                axis.ticks.y = element_blank(),
                                axis.line.y =   element_blank(),
                                strip.background = element_blank(),
                                # strip.text.y =  element_text(angle = 0),
                                strip.text.y.left = element_text(angle = 0)
                                # strip.placement = "outside"
)


#The next few lines sort the species names alphabetically so that they are plotted in alphabetical order

tt<-levels(Group1_WindSpeeds_BySpeciesDitrbitution$species)
tt<-tt[order(tt,decreasing = FALSE)]

Group1_WindSpeeds_BySpeciesDitrbitution$species <- 
  factor(Group1_WindSpeeds_BySpeciesDitrbitution$species, 
         levels = tt, ordered = TRUE)


setwd("C:/Users/lembi/Desktop")

png(filename = "testLogHists.png", #for image: = 3508 x 2480 Pixel
    width = 13.5, height = 18, units = "cm", #pointsize = 12,
    bg = "white",  res = 1000)


            ggplot(Group1_WindSpeeds_BySpeciesDitrbitution,aes(WindSpeed))+
              geom_histogram(bins=32,alpha = 0.9, fill = "lightblue")+
              
              #facet_wrap(.~species,nrow= 18,strip.position = "left" )+scale_y_continuous(trans = "log10")+
              #  facet_wrap(species ~ ., strip.position = "left", ncol = 1)+#, scales = "free_y") +
              facet_grid(species ~ ., switch = "y")+  #,shrink = FALSE
              scale_y_continuous(trans = "log10")+
              xlab(expression("Wind speed " ~ (m~s^{-1}) ))+
              thm2
            
            #theme(strip.text.y.left = element_text(angle = 0))
            # theme_classic() +
            # theme(strip.text.y = element_text(angle = 0),
            #       legend.position = "none")
            # 
            #thm2

dev.off()


###############################################################################################################
#In the following ggplot code I have tried to plot log histograms or densities that are normalised by the number#
#of data points of each species. But the results are weird, so it needs work#

ggplot(Group1_WindSpeeds_BySpeciesDitrbitution,aes(x=WindSpeed,y=species))+
  #the line below normalises by dividing each species count with the sum count of all species
  #geom_histogram(aes(y=(..count../sum(..count..))), bins=32,alpha = 0.9, fill = "lightblue")+
  #the line below normalises by appying in each group independetely size rather than the whole data size
  #geom_histogram(aes(y = stat(density*width)), bins=32,alpha = 0.9, fill = "lightblue")+
  geom_density(aes(y = stat(density), height = ..density..), bins=32,alpha = 0.9, fill = "lightblue")+
  
  #facet_wrap(.~species,nrow= 18,strip.position = "left" )+scale_y_continuous(trans = "log10")+
  #  facet_wrap(species ~ ., strip.position = "left", ncol = 1)+#, scales = "free_y") +
  facet_grid(species ~ ., switch = "y")+  #,shrink = FALSE
  scale_y_continuous(trans = "log10")+
  xlab(expression("Wind speed " ~ (m~s^{-1}) ))+
  thm2



############################################################Plot ridgeline plots#################################################
###################################################NOTE THAT THESE PLOTS ARE NOT LOG DENSITY PLOTS###############################

#In the following code I plot two ridgleine plots


library(ggridges)
library(viridis)


#The next few lines sort species alphabetically so that they appear in the ridgeline plots in alphabetical order
tt<-levels(Group1_WindSpeeds_BySpeciesDitrbitution$species)
tt<-tt[order(tt,decreasing = TRUE)]

Group1_WindSpeeds_BySpeciesDitrbitution$species <- 
  factor(Group1_WindSpeeds_BySpeciesDitrbitution$species, 
         levels = tt, ordered = TRUE)

thm <- theme_minimal() + theme(text = element_text(size = 16),
                               axis.title.x = element_text(size = 16, colour = "black",angle=0,margin = margin(1*c(0.4,0,0,0), unit = "cm"))
                               ,axis.text.x = element_text(size = 16, colour = "black",angle=0,margin = margin(1*c(0.2,0,0,0), unit = "cm"))
                               ,axis.text.y =  element_text(size = 16, colour = "black",angle=0,margin = margin(1*c(0,0,0,0), unit = "cm"))
                               ,legend.position = "none"
)

################A ridgeline plot of wind speed densities per species with median as a vertical black line

setwd("C:/Users/lembi/Desktop")
png(filename = "ERA5 Wind Speeds By Species.png", 
    width = 25.7, height = 20.9, units = "cm", 
    bg = "white",  res = 1000)


                ggplot(Group1_WindSpeeds_BySpeciesDitrbitution) +
                  geom_density_ridges(
                    # jittered_points = TRUE,
                    # point_shape = '|', point_size = 2, point_alpha = 1, alpha = 0.7,
                    quantile_lines=TRUE,
                    quantile_fun=function(x,...)median(x),
                    aes(x = WindSpeed, y = species,
                        #group = interaction(Season, Species),
                        #fill = Season), alpha = 0.6) +
                        fill = species), alpha = 0.6) +
                  scale_fill_viridis_d(name = expression("Wind speed " ~ (m~s^{-1}) ), option = "D") + #, option = "C") +
                  ylab(NULL) +
                  xlab(expression("Wind speed " ~ (m~s^{-1}) ))+
                  thm


dev.off()


################A ridgeline plot of wind speed densities per species with the first 10% of wind speeds indicated with blue
################and the last 10% with red

setwd("C:/Users/lembi/Desktop")
png(filename = "ERA5 Wind Speeds By Species - Group1And2_3.png", 
    width = 25.7, height = 20.9, units = "cm", 
    bg = "white",  res = 1000)


      ggplot(Group1_WindSpeeds_BySpeciesDitrbitution, aes(x = WindSpeed, y = species, fill = factor(stat(quantile)))) +
              stat_density_ridges(
                geom = "density_ridges_gradient",
                calc_ecdf = TRUE,
                quantiles = c(0.1,0.5, 0.9),#quantiles = c(0.025,0.5, 0.975),
                quantile_lines=TRUE,
                # quantile_fun=function(x,...)median(x),
              ) +
              scale_fill_manual(
                name = "Probability", values = c("lightblue", "#E0E0E0", "#E0E0E0","indianred1"), #name = "Probability", values = c("#0000FFA0", "#E0E0E0", "#E0E0E0","#FF0000A0"),
                labels = c("(0, 0.025]", "(0.025, 0.975]", "(0.975, 1]")
              )+
              thm

dev.off()


library(tcltk)

fileName <- tclvalue(tkgetSaveFile())

write.csv(SpeciesSummaryTableFromElham, file = fileName,row.names=FALSE)








library(forcats)



setwd("C:/Users/lembi/Desktop")
png(filename = "Log Density plot - ERA5 Wind Speeds By Species - Group1And2_2.png", 
    width = 25.7, height = 20.9, units = "cm", 
    bg = "white",  res = 1000)


ggplot(Group1_WindSpeeds_BySpeciesDitrbitution) +
  geom_density_ridges(
    
    # jittered_points = TRUE,
    # point_shape = '|', point_size = 2, point_alpha = 1, alpha = 0.7,
    quantile_lines=TRUE,
    quantile_fun=function(x,...)median(x),
    aes(x = WindSpeed, y =fct_reorder( species,WindSpeed,.fun = median), #USE fct_reorder to reorder species by median 
        #group = interaction(Season, Species),
        #fill = Season), alpha = 0.6) +
        fill = species), alpha = 0.6) +
  #scale_x_log10() + ##For log plot
  scale_fill_viridis_d(name = expression("Wind speed " ~ (m~s^{-1}) ), option = "D") + #, option = "C") +
  ylab(NULL) +
  xlab(expression("Wind speed " ~ (m~s^{-1}) ))+
  thm


dev.off()





























