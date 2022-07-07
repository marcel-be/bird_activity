library(data.table)
library(fasttime)
library(plyr)
library(dplyr)
library(scales)
library(stringr)
library(ggplot2)
library(hms)
library(lubridate)
library(tRackIT)

rm(list=ls())
path<- "J:/rts/rts_activity/"

#############################################################################################################################
#### 1.1.   Save all 2021 A/P data in one folder (only once) ####

proj<- getProject(projroot = "J:/rts/2021/filter_tRackIT_2021/", plot=F) # filtered data 2021 (change path if needed)

for(id in 1:nrow(proj$tags)){
  print(proj$tags$ID[id])
  fls<-list.files(paste0("J:/rts/2021/bird_data_active_passive/filter_tRackIT_2021/data/individuals/", proj$tags$ID[id] ,"/classification/"), full.names=TRUE)
  fls<- fls[fls %like% "aggregated"]
  
  lapply(fls, function(x){
    data<-fread(x) 
    name<- str_split(basename(x), "_aggregated")[[1]][1]
    fwrite(data,paste0(path, "birds_activity_1min_2021/", name, ".csv"))
  })
}

#############################################################################################################################
#### 1.2.   Save all 2020 A/P data in one folder (only once) ####

tags <-read.csv("J:/rts/2020/birds_enddates_2020_tables/birds_2020_enddate_1104202222.csv") 
tags <- tags %>% 
  filter(remove!="yes")


for(id in 1:nrow(tags)){
  print(tags$ID_20[id])
  fls<-list.files(paste0("J:/rts/2020/filter_tRackIT_2020/data/individuals/", tags$ID_20[id] ,"/classification/"), full.names=TRUE) # filtered data 2020 (change path if needed)
  fls<- fls[fls %like% "aggregated"]
  
  lapply(fls, function(x){
    data<-fread(x) 
    name<- str_split(basename(x), "_aggregated")[[1]][1]
    fwrite(data,paste0(path, "bird_data_active_passive/birds_activity_1min_2020/", name, ".csv"))
  })
}


#############################################################################################################################
#### 1.3.   Save all 2019 A/P data in one folder (only once) ####

proj<- getProject(projroot = "J:/rts/2019/filter_tRackIT_2019_2/", plot=F) 

for(id in 1:nrow(proj$tags)){
  print(proj$tags$ID[id])
  fls<-list.files(paste0("J:/rts/2019/filter_tRackIT_2019_2/data/individuals/", proj$tags$ID[id] ,"/classification/"), full.names=TRUE) # filtered data 2019 (change path if needed)
  fls<- fls[fls %like% "aggregated"]
  
  lapply(fls, function(x){
    data<-fread(x) 
    name<- str_split(basename(x), "_aggregated")[[1]][1]
    fwrite(data,paste0(path, "bird_data_active_passive/birds_activity_1min_2019/", name, ".csv"))
  })
}

############################################################################################################################
#### 1.4. Apply new enddates of tags (only once)

## Visual inspection of A/P Data suggested changing the end date (signals end earlier)

end <-read.csv(paste0(path, "scripts/enddate.csv")) # dataframe with new enddates (then cut from a/p data from steps 1.1 - 1.3)
enddate<-end %>% 
  mutate(enddate_2 = as.POSIXct(enddate, format = "%d.%m.%Y %H:%M", tz="utc"))

for(i in 1:nrow(enddate)){
  
  print(enddate$ID[i])
  
  fls<-paste0(path, "bird_data_active_passive/birds_activity_1min_", enddate$year[i], "/", enddate$ID[i] ,".csv")
  
  lapply(fls, function(x){
    data<-fread(x) 
    data$timestamp <- fasttime::fastPOSIXct(data$timestamp, tz="UTC")
    data<-data[data$timestamp < enddate$enddate_2[i]] 
    fwrite(data,x)
  })
}


## remove certain times from 210518_150077_40 (for now here in A/P data! redo this in raw data ASAP!) only once
df <-read.csv(paste0(path, "bird_data_active_passive/birds_activity_1min_2021/210518_150077_40.csv"))
df<- df %>% 
  filter(timestamp > "2021-06-01 06:00")%>% 
  filter(timestamp < "2021-06-08 06:00" | timestamp > "2021-07-01 00:00")
write.csv(df, paste0(path, "birds_activity_1min_2021/210518_150077_40.csv"))



#############################################################################################################################
#### 2.   Preparation of A/P Data ####

## Create file with all birds; one unique row for each bird --> Overview of all tags used 
## This data will NOT be used for analysis (see Step 3)
## Problem: Each year uses different Bird-IDs and different tables for metadata --> create tables for each year and merge later to one big table ##Get Ring-ID from different metadata-tables. Ring-ID then represents a uniquely assignable ID for all further merging purposes (getting "sex", "age" etc from big-bird-metatable "Dataset_Birds-20220118.csv")

#### 2.1. Prepare data for each year (2019/2020/2021)  ####

#2019:

Data.file.list_2019=list.files(paste0(path,"bird_data_active_passive/birds_activity_1min_2019/"), full.names = TRUE) 

tags_2019_1 <- plyr::ldply(Data.file.list_2019, function(x){
  data<- fread(x)
  data$ID<-str_split(basename(x), ".csv")[[1]][1]
  data$file<-basename(x) # modify
  data$timestamp<-as.POSIXct(data$timestamp, tz="UTC")
  # Calculate start and end dates
  data$start_datetime<-min(data$timestamp)
  data$stop_datetime<-max(data$timestamp)
  # Add year
  data$year <- year(data$timestamp)
  # Compile
  data <- data %>% select(ID,file,year,start_datetime,stop_datetime) %>% unique()
  return(data.table(data))
}) 

# Add species & ring-number information for 2019; Ring number is important for merging the df with metadata (next step):

df <- fread(paste0(path, "bird_metadata/Sendervoegel2019.csv"))  # table of tagged birds 2020

meta_2019 <- df %>% 
  select(ID,Species_DE, Ring_ID, Date) %>% 
  rename(species_de=Species_DE,
         ring_ID=Ring_ID,
         date_capture=Date) %>% 
  mutate(date_capture=as.Date(date_capture, format =  "%d/%m/%Y")) # may need adaption

tags_2019_2 <- merge(tags_2019_1, meta_2019, by="ID", all.x=T)


#2020:

Data.file.list_2020=list.files(paste0(path,"bird_data_active_passive/birds_activity_1min_2020/"), full.names = TRUE) 

tags_2020_1 <- plyr::ldply(Data.file.list_2020, function(x){
  data<- fread(x)
  data$ID<-str_split(basename(x), ".csv")[[1]][1]
  data$file<-basename(x) # modify
  data$timestamp<-as.POSIXct(data$timestamp)
  # Calculate start and end dates
  data$start_datetime<-min(data$timestamp)
  data$stop_datetime<-max(data$timestamp)
  # Add year
  data$year <- year(data$timestamp)
  # Compile
  data <- data %>% select(ID,file,year,start_datetime,stop_datetime) %>% unique()
  return(data.table(data))
}) 

# Add species & ring-number information for 2020; Ring number is important for merching the df with metadata (next step):

df <- fread(paste0(path, "bird_metadata/Aktivitaet_Sender_220420.csv"))  # table of tagged birds 2020

meta_2020 <- df %>% 
  select(ID,Vogelart, Ringnummer, Datum) %>% 
  rename(species_de=Vogelart,
         ring_ID=Ringnummer,
         date_capture=Datum) %>% 
  mutate(date_capture=as.Date(date_capture, format =  "%d/%m/%Y")) # may need adaption

tags_2020_2 <- merge(tags_2020_1, meta_2020, by="ID", all.x=T)


# 2021:

Data.file.list_2021=list.files(paste0(path,"bird_data_active_passive/birds_activity_1min_2021/"), full.names = TRUE) 

tags_2021_1 <- plyr::ldply(Data.file.list_2021, function(x){
  data<-data.table::fread(x)
  data$ID<-str_split(basename(x), ".csv")[[1]][1]
  data$file<-basename(x) 
  #data$timestamp<-as.POSIXct(data$timestamp)
  # Calculate start and end dates
  data$start_datetime<-min(data$timestamp)
  data$stop_datetime<-max(data$timestamp)
  # Add year
  data$year <- year(data$timestamp)
  # Compile
  data <- data %>% select(ID,file,year,start_datetime,stop_datetime) %>% unique()
  return(data.table(data))
}) 

## Add species & ring-number information 2021; Ring number is important for merging the df with metadata (next step):

df <- fread(paste0(path, "bird_metadata/tags_2021_metadata.csv"))  # table of tagged birds 2020

meta_2021 <- df %>% 
  select(ID_21,Species_DE, Ring_ID, Date) %>% 
  rename(species_de=Species_DE,
         ring_ID=Ring_ID,
         date_capture=Date,
         ID=ID_21) %>% 
  mutate(date_capture=as.Date(date_capture, format =  "%d/%m/%Y")) # may need adaption

tags_2021_2 <- merge(tags_2021_1, meta_2021, by="ID", all.x=T)
str(tags_2020_2$date)


#############################################################################
#### 2.2 Combine all years and add metadata

## merge 2019 and 2020 and 2021:
tags_all_1<- rbind(tags_2019_2, tags_2020_2, tags_2021_2)

# Add metadata information for all years: 
df <- fread(paste0(path, "bird_metadata/Dataset_Birds-20220118.csv"))  # big bird-capture-table (data on "sex", "breeding", "size"...)

metadata <- df %>% 
  filter(Tag=="yes") %>% 
  #filter(Gps=="no") %>% 
  select(Ring_ID, Date, Species_engl, Species_lat, sex, age, brood_patch, breeding, family, weight_g, wing_mm, P8_mm) %>%  # add more if needed
  rename(species_en=Species_engl,
         species_lat=Species_lat,
         ring_ID=Ring_ID,
         date_capture=Date)%>% 
  mutate(date_capture=as.Date(date_capture, format =  "%d/%m/%Y")) # may need adaption (data format "/" or "-")

tags_all_2 <- tags_all_1 %>% 
  left_join(metadata, by=c("ring_ID", "date_capture")) %>% 
  mutate(time_total= (stop_datetime - start_datetime))


##########################################################################
#### 2.3. out-filter faulty tags

tags_all_3 <- tags_all_2 %>% 
  filter(ID!="150007_20") %>% # thrush with day-night-shift (predated)
  filter(ID!="150156_4_30_1") %>%  # Jay with missing nights
  filter(ID!="210610_150099_20") %>%    # buzzard 
  droplevels()


write.csv(tags_all_3, paste0(path, "bird_data_storage/tags_overview.csv")) # will be used as "metadata" for rethomics package


###################################################################################
#### 3. Create file with all A/P data combined ####

## one row for each classified minute for all tags
## this data will be used for further analysis
## Problem: Each year uses different Bird-IDs and different tables for metadata --> create tables for each year and merge later to one big table

#### 3.1. Prepare data for each year (2019/2020/2021)  ####

#2019:

df_1min_2019 <- plyr::ldply(Data.file.list_2019, function(x){
  data<- fread(x)
  data$ID<-str_split(basename(x), ".csv")[[1]][1]
  data$file<-basename(x) # modify
  data$timestamp<- as.POSIXct(data$timestamp)
  data$start_datetime<- min(data$timestamp)
  #data$stop_datetime<-max(data$timestamp)
  data$date <- date(ymd_hms(data$timestamp))
  data$year <- year(data$date)
  data$ydate <- yday(data$date) # Day of the year
  data$month <- month(data$date) # Day of the month
  data$t <- data$timestamp-data$start_datetime # time since tagging in seconds
  data$activity <- ifelse(data$prediction=="active",1,0) # Activity classification
  data <- data %>% select(ID,timestamp,date,year,ydate,month,t,activity)
  return(data)
}) # Generate 1 dataframe with all IDs


#2020:

df_1min_2020 <- plyr::ldply(Data.file.list_2020, function(x){
  data<- fread(x)
  data$ID<-str_split(basename(x), ".csv")[[1]][1]
  data$timestamp<- as.POSIXct(data$timestamp)
  data$start_datetime<- min(data$timestamp)
  #data$stop_datetime<-max(data$timestamp)
  data$date <- date(ymd_hms(data$timestamp))
  data$year <- year(data$date)
  data$ydate <- yday(data$date) # Day of the year
  data$month <- month(data$date) # Day of the month
  data$t <- data$timestamp-data$start_datetime # time since tagging in seconds
  data$activity <- ifelse(data$prediction=="active",1,0) # Activity classification
  data <- data %>% select(ID,timestamp,date,year,ydate,month,t,activity)
  return(data)
}) # Generate 1 dataframe with all IDs


# 2021:

df_1min_2021 <- plyr::ldply(Data.file.list_2021, function(x){
  data<-data.table::fread(x)
  data$ID<-str_split(basename(x), ".csv")[[1]][1]
  #data$timestamp<-as.POSIXct(data$timestamp)
  data$start_datetime<-min(data$timestamp)
  #data$stop_datetime<-max(data$timestamp)
  data$date <- date(ymd_hms(data$timestamp))
  data$year <- year(data$date)
  data$ydate <- yday(data$date) # Day of the year
  data$month <- month(data$date) # Day of the month
  data$t <- data$timestamp-data$start_datetime # time since tagging in seconds
  data$activity <- ifelse(data$prediction=="active",1,0) # Activity classification
  data <- data %>% select(ID,timestamp,date,year,ydate,month,t,activity)
  return(data)
}) # Generate 1 dataframe with all IDs


#### 3.2.  Merge all years and save without metadata (used for rhetomic analysis only)

df_1min_all <- rbind(df_1min_2019, df_1min_2020, df_1min_2021) # no metadata!
df_1min_all<- df_1min_all %>% as.data.table()

df_1min_all_2 <- df_1min_all %>% 
  filter(ID!="150007_20")

fwrite(df_1min_all_2, paste0(path,"bird_data_storage/tags_1min_nometa.csv")) # for rhetomics analysis


#### 3.2.  Merge all years and add metadata

## Problem: Each year uses different Bird-IDs and different tables for metadata --> Get Ring-ID from different metadata-tables. Ring-ID then represents a uniquely assignable ID for all further merging purposes (getting "sex", "age" etc from big-bird-metatable "Dataset_Birds-20220118.csv")

# species & ring-number (for merging purposes) for 2019 Data:
df_1min_2019_meta <- df_1min_2019 %>% 
  left_join(meta_2019, by=c("ID"))%>% 
  as.data.table()

# species & ring-number (for merging purposes) for 2020 Data:
df_1min_2020_meta <- df_1min_2020 %>% 
  left_join(meta_2020, by=c("ID"))%>% 
  as.data.table()

# species & ring-number (for merging purposes) for 2021 Data
df_1min_2021_meta <- df_1min_2021 %>% 
  left_join(meta_2021, by=c("ID"))%>% 
  as.data.table()

# merge 2019, 2020 and 2021 data
df_1_min_all_meta<- rbind(df_1min_2019_meta, df_1min_2020_meta, df_1min_2021_meta)

# get metadata from big bird table (now clearly assigned by "Ring ID")
df_1_min_all_meta_2 <- df_1_min_all_meta %>% 
  left_join(metadata, by=c("ring_ID", "date_capture"))



#### 3.3. Filter-out faulty tags 

df_1_min_all_meta_3  <- df_1_min_all_meta_2  %>% 
  filter(ID!="150007_20") %>% # thrush with day-night-shift (predated)
  filter(ID!="150156_4_30_1") %>%  # Jay with missing nights
  filter(ID!="210610_150099_20") %>%    # buzzard 
  droplevels()



#### 3.4.  Add coverage data 
## count datapoints per day (in minutes) and divide by total number of minutes per day

coverage<- df_1_min_all_meta_3 %>% 
  count(ID, date) %>% 
  mutate(min_day = 60*24,
         coverage_daily= n/min_day) %>% 
  select(ID, date, coverage_daily)

nrow(coverage[coverage$coverage_daily<0.8]) / nrow(coverage[coverage$coverage_daily>=0.8]) # % of data if coverage of 80% used as a threshold

df_1_min_all_meta_4 <- df_1_min_all_meta_3 %>% 
  left_join(coverage, by=c("ID", "date"))


#### 3.5. Calculate sunset and sunrise

#df_1_min_all_meta_4$date_suntime<- format(df_1_min_all_meta_4$date, format="%Y/%m/%d")
#df_1_min_all_meta_4$sunrise<- NA
#df_1_min_all_meta_4$sunset<- NA
#for(i in 1:nrow(df_1_min_all_meta_4)){  # takes long! find faster way...
#  sun<- StreamMetabolism::sunrise.set(50.844, 8.661, df_1_min_all_meta_4$date_suntime[i], timezone="UTC")
#  df_1_min_all_meta_4$sunset  <- sun$sunset
#  df_1_min_all_meta_4$sunrise <- sun$sunrise
#}

df_1_min_all_meta_4$data_ID<- seq(1:nrow(df_1_min_all_meta_4)) # get unique ID for each datapoint

Lat<-50.844868
Lon<-8.663649

# apply function "time_of_day" from tRackIT-Package. Based on function StreamMetabolism::sunrise.set (check help for infos):
df <-time_of_day(data = df_1_min_all_meta_4, Lat=Lat, Lon = Lon, tcol = "timestamp", tz = "CET", activity_period = "diurnal") %>% 
  select(data_ID, timestamp, sunrise, sunset, time_to_rise, time_to_set, start_datetime, stop_datetime) %>% 
  rename(timestamp_CET = timestamp) %>% 
  right_join(df_1_min_all_meta_4, by="data_ID")

# For some reason some dates (28/03/2020, 28/03/2021, 29/03/2020, 29/03/2021) are not calculated. Manually set sunset/rise for these dates:

df$timestamp_CET[is.na(df$timestamp_CET)]<- with_tz(df$timestamp[is.na(df$timestamp_CET)], "CET") # timestamp for German timezone (+2h)

df$date_CET <- date(ymd_hms(df$timestamp_CET)) # date for German timezone

# get times for sunrise/set from www.timeanddate.com:
df$sunrise[is.na(df$sunrise) & df$date_CET=="2020-03-28"]<-  "2020-03-28 07:10:35" # replace missing values
df$sunrise[is.na(df$sunrise) & df$date_CET=="2021-03-28"]<-  "2021-03-28 07:10:35"
df$sunrise[is.na(df$sunrise) & df$date_CET=="2020-03-29"]<-  "2020-03-29 07:07:35"
df$sunrise[is.na(df$sunrise) & df$date_CET=="2021-03-29"]<-  "2021-03-29 07:07:35"

df$sunset[is.na(df$sunset) & df$date_CET=="2020-03-28"]<-  "2020-03-28 19:50:20"
df$sunset[is.na(df$sunset) & df$date_CET=="2020-03-29"]<-  "2020-03-29 19:52:20"
df$sunset[is.na(df$sunset) & df$date_CET=="2021-03-28"]<-  "2021-03-28 19:50:20"
df$sunset[is.na(df$sunset) & df$date_CET=="2021-03-29"]<-  "2021-03-29 19:52:20"

# calculate time since/from sunrise and sunset
df$time_to_rise[is.na(df$time_to_rise)]<- (df$timestamp_CET[is.na(df$time_to_rise)] - df$sunrise[is.na(df$time_to_rise)]) / (60*60)
df$time_to_set[is.na(df$time_to_set)]<- (df$timestamp_CET[is.na(df$time_to_set)] - df$sunset[is.na(df$time_to_set)]) / (60*60)


#### 3.6. correction for varying daylengths over the year

df<- df %>% 
  mutate(daylength = sunset - sunrise,
         time_to_rise_std = (as.numeric(mean(daylength)) / as.numeric(daylength)) * time_to_rise) # multiply the time to/since sunset with the quotient of mean daylength and the actual daylength (on shorter days, the time to sunrise becomes a bit longer and vice versa)

#### 3.7 save final dataset

df_1_min_all_meta_5 <- df 
fwrite(df_1_min_all_meta_5, paste0(path,"bird_data_storage/tags_1min_withmeta.csv")) # all data 2019 & 2020 & 2021 with metadata
