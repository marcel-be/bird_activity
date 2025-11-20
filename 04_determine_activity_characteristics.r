########################################################################################################################## 4. Extract timing of activity values (from final daily gam curve) ####

## use data obtained from GAM (predicted values) to determine certain activity characteristics like "activity at sunrise" etc

library(tidyverse)
library(data.table)

path<- "J:/rts/rts_activity/"
path<- "D:/rts_activity/"
path<- "/Users/pandadamda/rts_activity/"

## load data
data_new<- fread(paste0(path,"data/bird_data_storage/models/model_prediction_dates.csv"), stringsAsFactors=T)
data_new$date_f<- as.factor(as.character(data_new$date_f)) 

## activity values from GAM model (Individuals):
# build a new dataframe (from "df1" on) and fill it with the activity characteristics

## 1. activity at time_to_rise_std = 0 (sunrise)
df1<- data_new %>% 
  filter(time_to_rise_std == 0) %>% 
  group_by(species_en, ID, ring_ID, date_f) %>% 
  reframe(act_at_sunrise = mu) %>% 
  distinct()


## 2. area under the curve
data_new_lst<- split(data_new, list(data_new$ring_ID, data_new$date_f),drop=T)

data_new_lst<- lapply(data_new_lst, function(x){
  cbind(x, auc = DescTools::AUC(x=x$time_to_rise_std, y=x$mu))
})

df_auc<- bind_rows(data_new_lst)

df4<- df_auc %>% 
  group_by(ring_ID, ID,  date_f) %>% 
  summarise(auc = mean(auc)) %>% 
  as.data.frame() %>% 
  distinct()%>% 
  left_join(df3,., by=c("ring_ID","ID","date_f"))

## 2.2. area under the curve for subset sunrise+-1h
data_new_sunrise<- data_new %>% 
  filter(time_to_rise_std<1 & time_to_rise_std >(-1))
data_new_sunrise_lst<- split(data_new_sunrise, list(data_new_sunrise$ring_ID, data_new_sunrise$date_f),drop=T)
data_new_sunrise_lst<- lapply(data_new_sunrise_lst, function(x){
  cbind(x, auc = DescTools::AUC(x=x$time_to_rise_std, y=x$mu))
})
df_auc_sunrise<- bind_rows(data_new_sunrise_lst)

df4<- df_auc_sunrise %>% 
  group_by(ring_ID, ID,  date_f) %>% 
  summarise(auc_sunrise = mean(auc)) %>% 
  as.data.frame() %>% 
  distinct()%>% 
  left_join(df4,., by=c("ring_ID","ID","date_f"))


## 3. finding maximum and minimum slope (potential measures for "start" and "end" of activity)

#  calculate differences (increase) in mu for each datapoint for each bird and date (in comparison to the respective, former datapoint) --> highest value for this difference represents the strongest increase in mu (predicted value) ergo the steepest point of the activity slope
# "diff"-function calculates the difference to the former datapoint
# calculating "diff" for the whole dataset results in big diff-values between day-shifts (steepest as/descent at first or last timestamp of the day) --> needs to be calculated for each day seperately
# new dataframe "data_new_diff" will have less datapoints than "data_new", because for the first datapoint of each bird individual and date, no difference in mu can be calculated. This applies only for the very first datapoint of each individual and date (out of thousands) and can be neglected
# no idea if calculating the CI makes any sense here...

data_new_diff<- data.frame()

for(i in 1:nlevels(data_new$ring_ID)){ # find a nicer, vectorized approach! (Problem, that nrow[output df] < nrow[data_new])
  
  diff_1 <- data_new %>% 
    filter(ring_ID==levels(ring_ID)[i]) %>% 
    droplevels()
  
  for(j in 1:nlevels(diff_1$date_f)){
    
    diff_2 <- diff_1 %>% 
      filter(date_f==levels(date_f)[j]) %>% 
      droplevels()
    
    diff_3<- slice_tail(diff_2, n=(nrow(diff_2)-1)) %>% 
      mutate(diff=diff(diff_2$mu))
    
    data_new_diff<- rbind(diff_3,data_new_diff)
  }
}

## in the next step the steepest slope (as- and descent) will be calculated. This slope will be determined within the time: time_to_rise_std >= -1.5 & time_to_rise_std <= 2.25 for the ascent. time_to_rise_std >= 13 & time_to_rise_std <= 17 for the descent
## a visual inspection of the produced values for the steepest ascend and descend showed that the range of daytime (in hours to sunrise) in which the steepest ascend or descend is chosen as activity start or end, differs for some days (the activity starts later or earlier, compared to average), these days must be filtered separately to get the correct values --> The highest "diff" is not always within the range of  [time_to_rise_std >= -1.5 & time_to_rise_std <= 2.25] or [time_to_rise_std >= 13 & time_to_rise_std <= 17]
## all this as actually needed just because the curves are so wobbly...

# steepest_ascend:
df5<- data_new_diff %>% 
  group_by(ring_ID, ID, date_f) %>% 
  filter(case_when(ring_ID=="6415101" & date_f=="2020-05-28"  ~ time_to_rise_std <= -3,# < -3
                   ring_ID=="6415104" & date_f=="2020-05-31"  ~ time_to_rise_std >=  3 & time_to_rise_std <= 5,#>3 <5
                   ring_ID=="6415104" & date_f=="2020-06-07"  ~ time_to_rise_std >=  0 & time_to_rise_std <= 2,#>0 <2
                   ring_ID=="7974327" & date_f=="2019-06-13"  ~ time_to_rise_std <=  1, #<1  
                   ring_ID=="7974426" & date_f=="2021-07-08"  ~ time_to_rise_std >= -2 & time_to_rise_std <= 0, # >-2 <0
                   ring_ID=="90850008" & date_f=="2020-05-11" ~ time_to_rise_std >= -5 & time_to_rise_std <= 3, # <3 >-5
                   ring_ID=="7974426" & date_f=="2021-08-03"  ~ time_to_rise_std >= -2 & time_to_rise_std <= 0, #>-2 <0
                   ring_ID=="90850088" & date_f=="2021-06-11" ~ time_to_rise_std >= -2 & time_to_rise_std <= 0, # <0 >-2
                   ring_ID=="6415104" & date_f=="2020-06-07" ~ time_to_rise_std >= -2 & time_to_rise_std <= 2, # <0 >-2
                   ring_ID=="7974327" & date_f=="2020-04-21" ~ time_to_rise_std >= -2 & time_to_rise_std <= 1, #
                   ring_ID=="7974327" & date_f=="2020-04-19" ~ time_to_rise_std >= -2 & time_to_rise_std <= 1, #
                   ring_ID=="90850086" & date_f=="2021-06-10" ~ time_to_rise_std >= 2 & time_to_rise_std <= 7, 
                   ring_ID=="90850086" & date_f=="2021-06-11" ~ time_to_rise_std >= 2 & time_to_rise_std <= 7, 
                   T ~ time_to_rise_std >= -1.5 & time_to_rise_std <= 2.25) # all other
  ) %>% 
  filter(diff == max(diff,na.rm=TRUE)) %>% 
  reframe(steepest_ascend = time_to_rise_std) %>% 
  as.data.frame() %>% 
  distinct()%>% 
  left_join(df4,., by=c("ring_ID","ID", "date_f"))


# steepest_descend
df6<- data_new_diff %>% 
  group_by(ring_ID,ID, date_f)  %>% 
  filter(case_when(ring_ID=="6412718" & date_f=="2020-05-07"  ~ time_to_rise_std >= 14,  # >14
                   ring_ID=="6415105" & date_f=="2020-06-16"  ~ time_to_rise_std >= 14, # >14
                   ring_ID=="7974327" & date_f=="2020-05-05"  ~ time_to_rise_std >= 15,# >15
                   ring_ID=="7974327" & date_f=="2020-05-20"  ~ time_to_rise_std >= 15, # >15
                   ring_ID=="7974327" & date_f=="2020-06-29"  ~ time_to_rise_std >= 15, # >15
                   ring_ID=="7974426" & date_f=="2021-07-27"  ~ time_to_rise_std >= 15, # >15
                   ring_ID=="81948671" & date_f=="2019-07-28" ~ time_to_rise_std >= 15, #>15
                   ring_ID=="81948728" & date_f=="2019-07-28" ~ time_to_rise_std >= 11 & time_to_rise_std <= 14, # >11 <14
                   ring_ID=="90619323" & date_f=="2020-06-03" ~ time_to_rise_std >= 15, #>15
                   ring_ID=="90786807" & date_f=="2021-08-01" ~ time_to_rise_std >= 15, # >15
                   ring_ID=="90786807" & date_f=="2021-07-30" ~ time_to_rise_std >= 10 & time_to_rise_std <= 14, #>10  <14
                   ring_ID=="90786807" & date_f=="2021-07-31" ~ time_to_rise_std >= 10 & time_to_rise_std <= 14, #>10  <14
                   ring_ID=="82298129" & date_f=="2021-05-05" ~ time_to_rise_std >= 7 & time_to_rise_std <= 10, # >7 <10
                   ring_ID=="90850007" & date_f=="2020-05-02" ~ time_to_rise_std >= 14, # >14
                   ring_ID=="90850019" & date_f=="2021-04-03" ~ time_to_rise_std >= 15, # >15
                   ring_ID=="90850019" & date_f=="2021-04-08" ~ time_to_rise_std >= 15, #> 15
                   ring_ID=="90850094" & date_f=="2021-06-29" ~ time_to_rise_std >= 15, #>15
                   ring_ID=="V188907" & date_f=="2020-06-07" ~ time_to_rise_std >= 10 & time_to_rise_std <= 15, #<15 >10
                   ring_ID=="7974327" & date_f=="2020-07-04" ~ time_to_rise_std >= 14.5 & time_to_rise_std <= 16, 
                   ring_ID=="7974402" & date_f=="2020-04-28" ~ time_to_rise_std >= 14 & time_to_rise_std <= 16, 
                   ring_ID=="7974408" & date_f=="2021-04-21" ~ time_to_rise_std >= 13.5 & time_to_rise_std <= 16,
                   ring_ID=="7974426" & date_f=="2021-07-11" ~ time_to_rise_std >= 14 & time_to_rise_std <= 16,
                   ring_ID=="7974426" & date_f=="2021-07-22" ~ time_to_rise_std >= 14 & time_to_rise_std <= 16,
                   ring_ID=="90850031" & date_f=="2020-08-03" ~ time_to_rise_std >= 14 & time_to_rise_std <= 16,
                   T ~ time_to_rise_std >= 13 & time_to_rise_std <= 17)
  )%>% 
  filter(diff == min(diff,na.rm=TRUE)) %>% 
  reframe(steepest_descend = time_to_rise_std) %>% 
  as.data.frame() %>% 
  distinct()%>%  
  left_join(df5,., by=c("ring_ID","ID","date_f"))

## some descents must be set as NA as thez make no sense (visual inspection)
df6[(df6$ring_ID=="90850031" & df6$date_f=="2020-08-07"),]$steepest_descend <- NA # set offset NA, because it makes no sense for this date
df6[(df6$ring_ID=="90850031" & df6$date_f=="2020-08-08"),]$steepest_descend <- NA # set offset NA, because it makes no sense for this date
df6[(df6$ring_ID=="7974327" & df6$date_f=="2019-06-22"),]$steepest_descend <- NA
df6[(df6$ring_ID=="7974426" & df6$date_f=="2021-07-30"),]$steepest_descend <- NA
df6[(df6$ring_ID=="81948728" & df6$date_f=="2019-07-28"),]$steepest_descend <- NA
df6[(df6$ring_ID=="82298125" & df6$date_f=="2019-04-11"),]$steepest_descend <- NA
df6[(df6$ring_ID=="82298129" & df6$date_f=="2021-05-05"),]$steepest_descend <- NA
df6[(df6$ring_ID=="82298129" & df6$date_f=="2021-05-06"),]$steepest_descend <- NA
df6[(df6$ring_ID=="90619320" & df6$date_f=="2020-04-12"),]$steepest_descend <- NA
df6[(df6$ring_ID=="90786807" & df6$date_f=="2021-07-30"),]$steepest_descend <- NA
df6[(df6$ring_ID=="90850086" & df6$date_f=="2021-06-02"),]$steepest_descend <- NA
df6[(df6$ring_ID=="90850086" & df6$date_f=="2021-06-02"),]$steepest_ascend <- NA
df6[(df6$ring_ID=="90850086" & df6$date_f=="2021-06-06"),]$steepest_ascend <- NA







# 4. Activity at sunset
# the mean daylength that was used to standardize "time to rise" is 15.51328
sunset<- data_new[which(abs(data_new$time_to_rise_std - 15.51328) == min(abs(data_new$time_to_rise_std - 15.51328)))]$time_to_rise_std[1] # value that is closest to sunset 
df7<- data_new %>% 
  filter(time_to_rise_std == sunset) %>% 
  group_by(ring_ID,ID, date_f) %>% 
  reframe(act_at_sunset = mu) %>% 
  as.data.frame()  %>% 
  distinct()%>% 
  left_join(df6,., by=c("ring_ID","ID","date_f"))

# 5. remove certain days, based on visual observations of daily activity curves:

df8<- df7[-which(df7$ring_ID=="81948671" & df7$date_f=="2019-07-08"),]
df8<- df8[-which(df8$ring_ID=="7974405"  & df8$date_f=="2021-05-31"),]
df8<- df8[-which(df8$ring_ID=="82298129" & df8$date_f=="2021-05-06"),]
df8<- df8[-which(df8$ring_ID=="90619198" & df8$date_f=="2019-06-10"),]
df8<- df8[-which(df8$ring_ID=="90619198" & df8$date_f=="2019-06-09"),]
df8<- df8[-which(df8$ring_ID=="90619320" & df8$date_f=="2020-04-17"),]
df8<- df8[-which(df8$ring_ID=="6415104"  & df8$date_f=="2020-06-03"),]
df8<- df8[-which(df8$ring_ID=="V188907"  & df8$date_f=="2020-06-07"),]


# 6. calculate relative activity at sunrise 
# Activity at sunrise is divided by total activity (AUC) to allow for comparability across individuals

df8$act_at_sunset_rel  <- df8$act_at_sunset / df8$auc
df8$act_at_sunrise_rel <- df8$act_at_sunrise/ df8$auc

# 7. calculate time length of activity
df8$length_act<- df8$steepest_descend - df8$steepest_ascend

df_act_charac <- df8

## safe file
fwrite(df_act_charac, paste0(path,"data/bird_data_storage/activity_characteristics/activity_characteristics_individual.csv"))

