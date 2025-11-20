library(tidyverse)
library(data.table)


path<- "/Users/pandadamda/rts_activity/"
path<- "E:/Uni_Arbeit/rts_activity/"


df_1min<- fread(paste0(path,"bird_data_storage/tags_1min_withmeta.csv"))

##########################################################################
## 1.1 pool woodpecker species (only "woodpecker" form now on)
df_1min$species_en<- gsub("Black_Woodpecker", "woodpecker",
                          gsub("Great_Spotted_Woodpecker", "woodpecker",
                               gsub("Middle_Spotted_Woodpecker", "woodpecker", df_1min$species_en)))

##########################################################################
## 1.2 change data format
df_1min$hour           <- hour(df_1min$timestamp) # Hour of the day
df_1min$minute         <- minute(df_1min$timestamp)
df_1min$month_f        <- factor(month(df_1min$date))
df_1min$year_f         <- factor(df_1min$year)
df_1min$ydate_f        <- as.factor(df_1min$ydate)
df_1min$species_en     <- as.factor(df_1min$species_en)
df_1min$ring_ID        <- as.factor(df_1min$ring_ID)
df_1min$date_CET       <- date(df_1min$date_CET)
df_1min$date_f         <- as.factor(df_1min$date_CET)
df_1min$ID             <- as.factor(df_1min$ID)
df_1min$week           <- week(df_1min$timestamp)
df_1min$brood_patch    <- as.factor(df_1min$brood_patch)
df_1min$time_of_day    <- as.numeric(as_hms(df_1min$timestamp_CET))/3600 # numeric time of day (not used)
df_1min$timestamp_CET  <- fasttime::fastPOSIXct(df_1min$timestamp_CET, tz="CET")
#df_1min$start_datetime<- fasttime::fastPOSIXct(df_1min$start_datetime, tz="CET") # only of needed in further analysis
#df_1min$stop_datetime <- fasttime::fastPOSIXct(df_1min$stop_datetime, tz="CET") # only of needed in further analysis


##########################################################################
## 1.3 set range of "time_to_rise" as -5 to 18. Not all Individuals cover the whole range from -8.7 to 22.1 (all data). This max and min values will be used for cc-smoother. -5 to 18 is the range that all bird individuals fall into (see script "range_time_to_rise.R")
nrow(df_1min[df_1min$time_to_rise_std<=18 & df_1min$time_to_rise_std>=-5,])/nrow(df_1min) # 95.4% of all data
df_1min<- df_1min %>%
  filter(time_to_rise_std >= -5) %>% 
  filter(time_to_rise_std <= 18)


##########################################################################
## 1.4 Subset of coverage > 75%:
nrow(df_1min[df_1min$coverage_daily>=0.75,])/nrow(df_1min) # 94.1 % of all data
df_1min<- df_1min %>% 
  filter(coverage_daily >= 0.75) %>% 
  droplevels()


##########################################################################
## 1.5 Subset of Tags with > 3 days of data (date of capture was removed already)
nrow(df_1min[df_1min$time_total >= 3 & df_1min$ID != "210408_150113_40" & df_1min$ID != "210622_150007_40",])/nrow(df_1min) # 0.99 %
df_1min<- df_1min %>% 
  filter(time_total >= 3) %>% 
  filter(ID != "210408_150113_40") %>% 
  filter(ID != "210622_150007_40")


##########################################################################
## 1.6. Subset of species with at least 4 individuals:
# make a visual inspection based on the following species and remove species with < 4 individuals
df_1min  %>% 
  count(species_en, ring_ID) %>%
  as.data.frame(.) %>% 
  count(species_en) %>%
  mutate(species_en = fct_reorder(species_en, n)) %>% 
  ggplot(. ,aes(y=species_en, x=n, label = round(n,digits=1))) +
  geom_bar(stat="identity", color="steelblue", fill="steelblue")+
  ggtitle("Count of tagged individuals \n(some individuals were tagged two times)")+
  geom_text(size = 5, position = position_stack(vjust = 1.025))+
  xlab("Count")+
  ylab("Species")+
  theme_minimal()+
  theme(text = element_text(size=15),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0, face="bold"),
        legend.title=element_blank(),
        legend.position = c(0.9, 0.1),
        legend.text = element_text(size=15))

df_1min<- df_1min %>% 
  filter(species_en=="Eurasian_Blackcap"|
           species_en=="Great_Tit"|
           species_en=="European_Robin" |
           species_en=="Common_Blackbird" |
           species_en=="Eurasian_Jay" |
           species_en=="Eurasian_Blue_Tit"| 
           species_en=="Common_Chaffinch"|
           species_en=="Wood_Warbler"
  ) %>% 
  droplevels() # subset of species with most individuals

##########################################################################
## 1.7. marking the start of each time series (start-timestamp of each Bird ID) as "TRUE". Will be incorporated as autocorrelation structure. Each ID is treated as an autocorrelated time series.
df_1min<- start_event(df_1min, column="timestamp_CET", event="ID") # Package "itsadug"

## 1.8. remove certain days, based on visual observations of daily activity curves:
df_1min<- df_1min[-which(df_1min$ring_ID=="81948671" & df_1min$date_f=="2019-07-08"),]
df_1min<- df_1min[-which(df_1min$ring_ID=="7974405"  & df_1min$date_f=="2021-05-31"),]
df_1min<- df_1min[-which(df_1min$ring_ID=="82298129" & df_1min$date_f=="2021-05-06"),]
df_1min<- df_1min[-which(df_1min$ring_ID=="90619198" & df_1min$date_f=="2019-06-10"),]
df_1min<- df_1min[-which(df_1min$ring_ID=="90619198" & df_1min$date_f=="2019-06-09"),]
df_1min<- df_1min[-which(df_1min$ring_ID=="90619320" & df_1min$date_f=="2020-04-17"),]
df_1min<- df_1min[-which(df_1min$ring_ID=="6415104"  & df_1min$date_f=="2020-06-03"),]
df_1min<- df_1min[-which(df_1min$ring_ID=="V188907"  & df_1min$date_f=="2020-06-07"),]


## save final dataframe
fwrite(df_1min, paste0(path, "data/bird_data_storage/tags_1_min_for_analysis.csv"))
