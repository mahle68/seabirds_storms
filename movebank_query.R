#playing around with REST API to identify suitable datasets on Movebank
#Oct. 1, 2020. Elham Noruani. Radolfzell, Germany.


library(tidyverse)
library(lubridate)

#"https://github.com/movebank/movebank-api-doc/blob/master/movebank-api.md"

sample <- read.csv("/home/mahle68/Desktop/direct-read (3)")

#find data that I can either download data or see data
#"query: https://www.movebank.org/movebank/service/direct-read?entity_type=study&i_am_owner=true&i_can_see_data=true&i_have_download_access=true"

data <- read.csv("/home/mahle68/Desktop/all_data_can_see_downlaod", stringsAsFactors = F) %>%
  mutate_at(vars(32:33), ~as.POSIXct(strptime(.,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  filter(grepl("GPS",sensor_type_ids) & 
           year(timestamp_first_deployed_location) %in% c(2017,2018) |
           year(timestamp_last_deployed_location) %in% c(2017,2018)) %>% 
  dplyr::select(name)

  
