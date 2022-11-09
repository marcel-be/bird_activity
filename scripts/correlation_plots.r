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
#df_1min$timestamp_CET  <- fasttime::fastPOSIXct(df_1min$timestamp_CET, tz="CET")
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
  group_by(species_en, interval) %>% # year, year_f, 
  summarise(mean_activity=mean(active_prop, na.rm = T)) %>% 
  pivot_wider(names_from = species_en, values_from = mean_activity) %>% 
  data.frame()

plot(df_10min.r$Common_Blackbird, df_10min.r$Eurasian_Blackcap)
plot(df_10min.r$Common_Blackbird, df_10min.r$Common_Chaffinch)




x<- cor(df_10min.r[,2:9], use = "pairwise.complete.obs")
corrplot::corrplot(x,method = "ellipse", type = "lower")


x<- corrr::correlate(df_10min.r[,3:10])#corrr::shave() %>% 
corrr::network_plot(x, min_cor = .01, repel = T)







# Hourly activity bouts ----
df_bouts <- df_1min %>% #select()
  group_by(ID, rdate, species_en, ydate, hour) %>% # year, year_f,
  summarise(activity_bout=rle(activity)$lengths,
            AP_value=rle(activity)$values) %>% 
  filter(AP_value==1 & activity_bout >5) %>%  # Only bouts > 5 min
  group_by(ID, rdate, species_en, ydate, hour) %>% # year, year_f, 
  summarise(mean_bout=mean(activity_bout, na.rm = T),
            sd_bout=sd(activity_bout, na.rm = T),
            CV_bout=sd(activity_bout, na.rm = T)/mean(activity_bout, na.rm = T),
            bout_nb=n()) %>% data.frame()
df_bouts$ydate_c <- scale(df_bouts$ydate, scale = F) # center linear date of the year
#df_bouts$hour_sc <- df_bouts$hour-12 # Noon as intercept


df_bouts.r <- df_bouts %>%  # Only bouts > 5 min
  group_by(species_en, rdate, hour) %>% # year, year_f, 
  summarise(mean_bout=mean(mean_bout, na.rm = T)) %>% 
  pivot_wider(names_from = species_en, values_from = mean_bout) %>% 
  data.frame()

df_bouts.r[,2:9] %>% log() %>%  cor(use = "pairwise.complete.obs") %>% 
  corrplot::corrplot(method = "ellipse", type = "lower")

x<- corrr::correlate(df_bouts.r[,3:10])#corrr::shave() %>% 
corrr::network_plot(x, min_cor = .01, repel = T)
