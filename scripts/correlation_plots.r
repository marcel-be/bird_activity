library(tidyr)
library(lubridate)
library(dplyr)
library(data.table)
library(hms)

path<- "J:/rts/rts_activity/"
path<- "G:/rts_activity/"
path<- "F:/Uni_Arbeit/rts_activity/"


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
  group_by(ID, ring_ID, rdate, species_en,year_f,brood_patch,sex, month,week,ydate, date_f,hour,age, weight_g, interval) %>% 
  summarise(n_intervals=length(activity),
            n_active=length(activity[activity==1]),
            n_passive=length(activity[activity==0]),
            time_of_day=mean(time_of_day),
            time_to_rise_std = mean(time_to_rise_std)) %>% 
  mutate(active_prop =  n_active/n_intervals,
         interval_test= round(time_to_rise_std, digits = 1))




df_10min.r <- df_10min %>%  
  group_by(ring_ID, interval_test) %>% # year, year_f, 
  summarise(mean_activity=mean(active_prop, na.rm = T)) %>% 
  pivot_wider(names_from = ring_ID, values_from = mean_activity) %>% 
  data.frame() %>% 
  na.omit()

df_10min.r <- df_10min %>%  
  group_by(species_en, interval_test) %>% # year, year_f, 
  summarise(mean_activity=mean(active_prop, na.rm = T)) %>% 
  pivot_wider(names_from = species_en, values_from = mean_activity) %>% 
  data.frame() %>% 
  na.omit()



x<- cor(df_10min.r[,2:9], use = "pairwise.complete.obs")
na.omit()
corrplot::corrplot(x,method = "ellipse", type = "lower")


x<- corrr::correlate(df_10min.r[,2:66])#corrr::shave() %>% 
corrr::network_plot(x, min_cor = .01, repel = T)


###########################

meta<- read.csv(paste0(path,"bird_data_storage/tags_overview.csv"), stringsAsFactors = T)
names<- levels(df_1min$ring_ID)
meta<- meta %>% 
  filter(ring_ID %in% names) %>% 
  select(-(1:7)) %>% 
  #select( -2, -4, -13,-14) %>% 
  select(1,10,11,12)
meta<- meta[which(duplicated(meta$ring_ID)==FALSE),]

df_1min.r <- df_1min %>%  
  group_by(ring_ID, time_to_rise_std) %>% # year, year_f, 
  summarise(mean_activity=mean(activity, na.rm = T)) %>% 
  pivot_wider(names_from = ring_ID, values_from = mean_activity) %>% 
  data.frame()


library(vegan)
df_10min.r <- df_10min.r %>%
  select(-1)

nmds <- metaMDS(df_10min.r, distance = "bray", k = 2, try = 20)
plot(nmds, display="species") # Eure Plots werden dargestellt
text(nmds, display="species",cex=0.6, pos=4, offset=0.3)

fit<- envfit(nmds, meta, permutations = 1000, na.rm = T) 
fit



decorana(veg = df_10min.r)
pca<- rda(df_10min.r)
ordiplot(pca, display = 'species', type = 'n')
points(pca,display = 'species')
text(pca, display="species",cex=0.6, pos=4, offset=0.3)



## Raphael:
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


df_bouts.r <- df_bouts %>%  # Only bouts > 5 min
  group_by(species_en, rdate, hour) %>% # year, year_f, 
  summarise(mean_bout=mean(mean_bout, na.rm = T)) %>% 
  pivot_wider(names_from = species_en, values_from = mean_bout) %>% 
  data.frame()

df_bouts.r[,3:10] %>% log() %>%  cor(use = "pairwise.complete.obs") %>% 
  corrplot::corrplot(method = "ellipse", type = "lower")

x<- corrr::correlate(df_bouts.r[,3:10])#corrr::shave() %>% 
corrr::network_plot(x, min_cor = .01, repel = T)
