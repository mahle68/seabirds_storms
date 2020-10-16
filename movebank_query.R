#playing around with REST API to identify suitable datasets on Movebank
#Oct. 1, 2020. Elham Noruani. Radolfzell, Germany.


library(tidyverse)
library(lubridate)

#"https://github.com/movebank/movebank-api-doc/blob/master/movebank-api.md"

#get a list of all movebank studies
#query https://www.movebank.org/movebank/service/direct-read?entity_type=study

data <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/movebank_query/all_data", stringsAsFactors = F) %>%
  mutate_at(vars(32:33), ~as.POSIXct(strptime(.,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  filter(grepl("GPS",sensor_type_ids)) %>% 
  mutate(yr_start = year(timestamp_first_deployed_location),
         yr_end = year(timestamp_last_deployed_location))

recent <- data %>% #extract studies that have data occurring in 2017 and 2018
  filter(yr_start <= 2017 & yr_end >= 2018) #%>% 
  #select(name)

#idea: download all recent data and filter out spatially in R. But, let's see for which studies I can downlaod data:
can_download <- recent %>% 
  filter(i_have_download_access == "true")

  

#find data that I can either download data or see data
#"query: https://www.movebank.org/movebank/service/direct-read?entity_type=study&i_am_owner=true&i_can_see_data=true&i_have_download_access=true"

data <- read.csv("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/movebank_query/can_see_download", stringsAsFactors = F) %>%
  mutate_at(vars(32:33), ~as.POSIXct(strptime(.,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  filter(grepl("GPS",sensor_type_ids) & 
           year(timestamp_first_deployed_location) %in% c(2017,2018) |
           year(timestamp_last_deployed_location) %in% c(2017,2018)) %>% 
  dplyr::select(name)

