#script for subsetting tracking data and retaining only segments with winds higher than x kmh
#Dec. 15, 2020. Elham Nourani. Radolfzell am Bodensee



#----------- STEP 1: open albatross data ----

#one-hourly data annotated by Sophie
data <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/data/From_Sophie/DATSETS/WAAL_allGPS_2010-2020_homo_R1h_TrackParam_Wind50kmh.csv",
                 na.strings = c("NA",""), fileEncoding="latin1") %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  #filter (FlyingSitting == "flying") %>% #filter out sitting positions
  arrange(TripID,date_time) %>% 
  #rowwise() %>% 
  mutate(row_id = row_number())

#----------- STEP 2: subset points with strong wind ----

windy <- data %>% 
  filter(windSpeed_kmh >= 50) %>% 
  mutate(delta = lead(row_id,1) - row_id ) 

#----------- STEP 3: assign segment IDs to consecutive points retained ----

windy$segID <- c(1,rep(NA,nrow(windy)-1))
  
for(i in 2:nrow(windy)){
  if(windy$delta[i-1] == 1){
    windy$segID[i] <- windy$segID[i-1]
  } else {
    windy$segID[i] <- windy$segID[i-1] + 1
  }
}

#----------- STEP 4: keep segments with at least n points ----

windy %>% 
  group_by(segID) %>% 
  mutate(n_pts = n()) %>% 
  filter(n_pts >= 5) %>% 
  dplyr::select(c(segID, n_pts))
