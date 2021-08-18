#script for looking at the entire nazca booby dataset, annotated directly in Movebank.
#these have not been filtred for flying, and haven't been filtered hourly

#open data
nb <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/annotation/nazca_b_all/Nazca booby Sula granti Isla Espanola, Galapagos.-7401245543257496617/Nazca booby Sula granti Isla Espanola, Galapagos.-7401245543257496617.csv") %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         u10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.,
         v10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component.) %>% 
  mutate(wind_speed_ms = sqrt(u10m^2 + v10m^2),
         wind_speed_kmh = sqrt(u10m^2 + v10m^2) * 3.6,
         month = month(timestamp),
         yday = yday(timestamp))
