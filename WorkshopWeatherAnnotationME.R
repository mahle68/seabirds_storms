install.packages('Rcpp')
library(Rcpp)

setwd("E:/SCGE WORKSHOP")

#### STEP 1: read your data
###########################
dats <- read.csv('hb-workshop.csv')	## enter correct filename
dats <- dats[,c("dev","datetime","lat","long")]

### STEP 2: define datetime format & additional data processing
################################################################
OHB$datetime <- as.POSIXct(strptime(OHB$datetime, format="%Y-%m-%d %H:%M:%S"), tz='UTC')  ### CHECK TIME ZONE!!!
#dats$datetime <- as.POSIXct(strptime(OHB$date.gmt, format="%d/%m/%y %H:%M:%S"), tz='UTC')  ### CHECK TIME ZONE!!!
## extract year and month columns for tracking data and nests
library(lubridate)
OHB$yr <- year(OHB$datetime)
OHB$mth <- month(OHB$datetime)
OHB$day <- day(OHB$datetime)
OHB$yday <- yday(OHB$datetime)
dats$season <- ifelse(dats$mth > 7,'autumn','spring')


OHB$yday2 <- paste('D',dats$yday,sep='')

# identify single migrations (dev,yr,season)
OHB$migr <- as.factor(paste(OHB$ptt, dats$yr, OHB$season, sep ='_', collapse = NULL))

OHB$dy<-paste(0,OHB$day, sep="")
OHB$mnth<-paste(0,OHB$month, sep="")
OHB$date<-paste(OHB$year,OHB$mnth,OHB$dy,sep="-")
OHB$datetime<-paste(OHB$date,OHB$time,sep=" ")
### STEP 3: define unique travel days per bird
################################################################
OHB$indday <- as.factor(paste(OHB$ptt, OHB$yday2, sep ='_', collapse = NULL))

# subset data to relevant columns after processing
data <- dats[,c("ptt","indday","yr","mth","day","season","migr","datetime","latitude","longitude")]

# order your data
data <- data[order(data$ptt,OHB$datetime),]

### STEP 4: determine first and last obs for each travel day
################################################################
starts <- aggregate(data$datetime, by=list(data$indday), FUN= function(x)min(as.character(x)))
colnames(starts)[1:2]  <- c("indday","daily.start")
starts$daily.start <- as.POSIXct(strptime(starts$daily.start,format="%Y-%m-%d %H:%M:%S"),tz="UTC")

ends <- aggregate(data$datetime, by=list(data$indday), FUN= function(x)max(as.character(x)))
colnames(ends)[1:2]  <- c("indday","daily.end")
ends$daily.end <- as.POSIXct(strptime(ends$daily.end,format="%Y-%m-%d %H:%M:%S"),tz="UTC")

# merge daily starts and ends with original dataframe
new <- merge(starts,ends)
data <- merge(data,new,all.x=T)

# extract locations at start and end of each day
startofdays <- which(data$datetime == data$daily.start)
endofdays <- which(data$datetime == data$daily.end)

starts <- data[startofdays,c("indday","longitude","latitude","datetime")]
ends <- data[endofdays,c("indday","longitude","latitude","datetime")]
colnames(starts)[1:4] <-c("indday","stlon","stlat","stdt")
colnames(ends)[1:4] <-c("indday","endlon","endlat","enddt")
daylocs <- merge(starts,ends,by="indday")

### STEP 5: CALCULATE DAILY DISTANCE AND DIRECTON 
##################################################
library(fossil)
# calculate daily distance
daylocs$daily.dist <- deg.dist(daylocs$stlon,daylocs$stlat,daylocs$endlon,daylocs$endlat)*1000

## IMPORTANT: we don?t calculate daily speed in this excercise because first and last fix per day 
## do not necessarily represent daily travel period

# calculate daily direction
daylocs$daily.dir <-  earth.bear(long1=daylocs$stlon, lat1=daylocs$stlat, long2=daylocs$endlon, lat2=daylocs$endlat)
# Adjust direction (on scale from 0 to 360) to GPS scale (on scale from -180 to 180) ##
daylocs$daily.dir <- ifelse(daylocs$daily.dir > 180, daylocs$daily.dir-360, daylocs$daily.dir)

# Calculate daily timerange
daylocs$daily.time <- as.numeric(difftime(daylocs$enddt,daylocs$stdt,units='hours'))

# append extra info from original dataset
daylocs <- merge(daylocs,unique(data[,c("indday","dev","daily.start","migr","season","yr","mth","day")]))

### STEP 6: CLASSIFY TRAVEL VS NON-TRAVEL DAYS 
##################################################
daylocs$travel <- ifelse(daylocs$daily.dist/1000 > 25,'travel','rest')

### INTERMEZZO: plot travel and non travel days on a map
#####################################
library(ggplot2);library(mapdata);library(grid)

longlimits <- c(110,140)
latlimits <- c(20,41)

# the following code produces maps with all first locs per day coloured for travel vs resting
# in different facets per bird and per season
ggplot(daylocs,aes(x=stlon,y=stlat,color=travel)) + 	borders("world",colour='grey70',fill='grey70',countries='FALSE') +
	geom_segment(aes(xend=endlon,yend=endlat),alpha=1,size = 1.4) +
 		theme_bw() + coord_fixed(xlim = longlimits, ylim = latlimits,ratio=1) +
	xlab("Longitude") +ylab("Latitude") +
	theme(plot.title=element_text(size=10),
		strip.text.x=element_text(size=10,color='white'),
		axis.title=element_text(size=10))+
    facet_grid(ptt~yr)

### STEP 7: SUBSET MIGRATION FROM TOTAL DATASET
##################################################
migration <- subset(daylocs,daylocs$travel == 'travel')

## IMPORTANT: this is a very rough classification. To do it well you need to take time to inspect histograms to find appropriate cut-off
## values for your species of interest. You may also want to classify migration vs staging based on multiple variables

### STEP 8: ANNOTATE WIND DATA FROM NCEP TO STARTING POINT OF EACH DAY
######################################################################
library(RNCEP)
migration$uwind <- NCEP.interp(variable='uwnd', level=925,lat=migration$stlat, lon=migration$stlon, dt=migration$daily.start,reanalysis2=TRUE, keep.unpacking.info=TRUE)
migration$vwind <- NCEP.interp(variable='vwnd', level=925,lat=migration$stlat, lon=migration$stlon, dt=migration$daily.start,reanalysis2=TRUE, keep.unpacking.info=TRUE)

### STEP 9: CALCULATE TAILWINDS ALONG TRAVEL DIRECTION
######################################################################
migration$tailwind <- (sqrt(migration$uwind^2 + migration$vwind^2)*cos(((atan2(migration$uwind,migration$vwind)*(180/pi))-migration$daily.dir)*(pi/180)))

## STEP 10: CORRELATE TAILWINDS WITH DAILY TRAVEL DISTANCE
##########################################################

## think of a suitable method (loess, (g)lm, R-squared, ...)


### INTERMEZZO: plot daily tracks and wind vectors on map
#####################################
scaler <- 0.3
ggplot() + 	borders("world",colour='grey70',fill='grey20',countries='FALSE') +
	geom_segment(data=migration,aes(x=stlon,y=stlat,xend=endlon,yend=endlat,col=factor(day)),size=.8,
		arrow=arrow(length=unit(0.1,"cm"),angle=30,ends = "last", type = "open")) +
  	geom_segment(data=migration,aes(x=stlon, y=stlat, xend=stlon+uwind*scaler, yend=stlat+vwind*scaler,col=factor(day)),size=.6,
		arrow=arrow(length=unit(0.1,"cm"),angle=30,ends = "last", type = "closed"),linetype='dashed') +
	theme_bw() + coord_fixed(xlim = longlimits, ylim = latlimits,ratio=1) +
	xlab("Longitude") +ylab("Latitude") +
	theme(plot.title=element_text(size=10),
		strip.text.x=element_text(size=10,color='white'),
		axis.title=element_text(size=10))


## STEP 11: OBTAIN WINDDATA FOR EACH DATE IN MOVEMENT DATA
##########################################################

#migration <- subset(migration,mth == 8)

winddata <- data.frame(dt= character(0), lon= numeric(0), lat = numeric(0), uw = numeric(0), vw = numeric(0))
for(i in 1:length(unique(migration$mth))){
	t.ss <- subset(migration, mth == unique(migration$mth)[i])

	month <- unique(t.ss$mth)
	year <- unique(t.ss$yr)
	minday <- min(t.ss$day)
	maxday <- max(t.ss$day)

	days <- seq(1,31,1)
	daystokeep <- unique(t.ss$day)
	daystoremove <- days [! days %in% daystokeep]

	uw<- NCEP.gather(variable = 'uwnd', level = 925, months.minmax=c(month), 
		years.minmax=c(year),lat.southnorth=c(0,55), lon.westeast=c(-20,15),
		reanalysis2 = TRUE,return.units = TRUE, status.bar=TRUE)
	uw<- NCEP.restrict(uw,hours2remove = c(0,6,18), days2remove = daystoremove, set2na = FALSE)
	uwdata<- NCEP.array2df(uw, var.names=NULL)
	colnames(uwdata) <- c("dt","lat","long","uw")

	vw<- NCEP.gather(variable = 'vwnd', level = 925, months.minmax=c(month),
		years.minmax=c(year),lat.southnorth=c(0,55), lon.westeast=c(-20,15),
		reanalysis2 = TRUE,return.units = TRUE, status.bar=TRUE)
	vw<- NCEP.restrict(vw,hours2remove = c(0,6,18), days2remove = daystoremove, set2na = FALSE)
	vwdata<- NCEP.array2df(vw, var.names=NULL)
	colnames(vwdata) <- c("dt","lat","long","vw")

	wind.temp <- merge(uwdata,vwdata)
	winddata <- rbind(wind.temp,winddata)
}

## STEP 12: OBTAIN DATE INFO FROM DATETIMES IN WINDDATA
##########################################################
dates <- unique(winddata$dt)
mths <- 'na'
for(i in 1:length(dates)){
	mths[i] <- unlist(strsplit(dates[i], split='_', fixed=TRUE))[2]
	}
yrs <- 'na'
for(i in 1:length(dates)){
	yrs[i] <- unlist(strsplit(dates[i], split='_', fixed=TRUE))[1]
	}
days <- 'na'
for(i in 1:length(dates)){
	days[i] <- unlist(strsplit(dates[i], split='_', fixed=TRUE))[3]
	}

cal <- as.data.frame(cbind(dates,yrs,mths,days))
colnames(cal)[1:4] <- c("dt","yr","mth","day")
winddata <- merge(winddata,cal,all.x=T)

## STEP 13: CREATE A DATE OBJECT
##########################################################
winddata$date <- paste(winddata$yr,winddata$mth,winddata$day,sep='-')
winddata$date <- as.Date(winddata$date)
migration$date <- paste(migration$yr,migration$mth,migration$day,sep='-')
migration$date <- as.Date(migration$date)
dats$date <- paste(dats$yr,dats$mth,dats$day,sep='-')
dats$date <- as.Date(dats$date)

## STEP 14: PRODUCE DAILY TRAJECTORIES IN WIND FIELDS
##########################################################
for(i in 1:length(unique(winddata$date))){
	t.ss <- subset(winddata, date == unique(winddata$date)[i])
	m.ss <- subset(dats, date == unique(winddata$date)[i])
	
	scaler <- 0.3
	p <- ggplot() + borders("world",colour='grey70',fill='grey70',countries='FALSE') +
		geom_path(data=m.ss,aes(x=long,y=lat,group=dev),size=1,col='red')+
  		geom_segment(data= t.ss, aes(x=long, y=lat, xend=long+uw*scaler, yend=lat+vw*scaler),size=.6,
				arrow=arrow(length=unit(0.1,"cm"))) +
  		theme_bw() + coord_fixed(xlim = longlimits, ylim = latlimits,ratio=1) +
  		theme(legend.position = 'right',
	  			legend.direction= 'vertical',
				panel.border	= element_rect(colour='white'),
				panel.grid.major	= element_line(colour='grey50',linetype='dashed'),
				panel.grid.minor	= element_line(colour='grey70',linetype='dashed'),
				strip.background	= element_rect(fill='white',colour='NA'),
				strip.text		= element_text(size=10,face='bold'),
				axis.text		= element_text(size=10,face='italic'),
				axis.title		= element_text(size=11))

	ggsave(plot=p,filename=as.character(paste("WindField",m.ss$date[i],m.ss$dev[i],".jpeg",sep="_")),height=6,width=6,dpi=300)
}









