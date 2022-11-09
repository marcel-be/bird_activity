library(tidyverse)
library(lubridate)
library(data.table)


df_1min<- fread(paste0(path, "bird_data_storage/tags_1_min_for_analysis.csv"))

df_1min$ydate_f        <- as.factor(df_1min$ydate)
df_1min$species_en     <- as.factor(df_1min$species_en)
df_1min$ring_ID        <- as.factor(df_1min$ring_ID)
df_1min$date_CET       <- date(df_1min$date_CET)
df_1min$date_f         <- as.factor(df_1min$date_CET)
df_1min$ID             <- as.factor(df_1min$ID)
df_1min$brood_patch    <- as.factor(df_1min$brood_patch)
df_1min$timestamp_CET  <- fasttime::fastPOSIXct(df_1min$timestamp_CET, tz="CET")
df_1min$rdate          <- date(df_1min$date_CET)-date(ymd_hms(df_1min$start_datetime))


## create 10-min-intervals for faster plotting:
df_10min <- df_1min %>% 
  mutate(species_en = as.factor(species_en),
         interval   = as_hms(floor_date(timestamp_CET, unit="10minutes"))) %>% 
  group_by(ID, ring_ID, rdate, species_en,year_f,brood_patch,sex, month,week,ydate, date_f,hour,interval) %>% 
  summarise(n_intervals=length(activity),
            n_active=length(activity[activity==1]),
            n_passive=length(activity[activity==0]),
            time_of_day=mean(time_of_day),
            time_to_rise_std = mean(time_to_rise_std)) %>% 
  mutate(active_prop =  n_active/n_intervals)




df_10min.r <- df_10min %>%  # Only bouts > 5 min
  group_by(species_en, rdate, hour) %>% # year, year_f, 
  summarise(mean_activity=mean(active_prop, na.rm = T)) %>% 
  pivot_wider(names_from = species_en, values_from = mean_activity) %>% 
  data.frame()


x<- cor(df_10min.r[,3:10], use = "pairwise.complete.obs")
corrplot::corrplot(x,method = "ellipse", type = "lower")


x<- corrr::correlate(df_10min.r[,3:10])#corrr::shave() %>% 
corrr::network_plot(x, min_cor = .01, repel = T)


